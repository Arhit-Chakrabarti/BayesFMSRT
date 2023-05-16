data_list <- c("dark_replicate_1",
               "dark_replicate_2",
               "light_replicate_1",
               "light_replicate_2")

# load gene count matrices
gene_data <-NULL
gene_data$data <- list()
ncells <- rep(NA, length(data_list))
gene_data$location <- list() # To compute the spatial locations

for (i in 1:length(data_list)){
  data <- read.csv(paste("~/BayesFMSRT/STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
  loc <- read.csv(paste("~/BayesFMSRT/STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("~/BayesFMSRT/STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


gene_data$reads_per_cell <-list()
median_expression <- 0

for(i in 1:length(gene_data$data)){
  ind <- which(apply(gene_data$data[[i]], 2, sum) > 100 & apply(gene_data$data[[i]], 2, sum) < 2000)
  gene_data$data[[i]] <- gene_data$data[[i]][ , ind]
  gene_data$location[[i]] <- gene_data$location[[i]][ind, ]
  gene_data$reads_per_cell[[i]] <- apply(gene_data$data[[i]], 2, sum)
  median_expression[i] <- median(gene_data$reads_per_cell[[i]])
  rownames(gene_data$data[[i]]) <- gene_data$name
}


gene_data$N <-  replicate(length(gene_data$data), "list", simplify = FALSE)

for(l in 1:length(gene_data$data)){
  gene_data$N[[l]] <- matrix(NA, nrow = dim(gene_data$data[[l]])[1], ncol = dim(gene_data$data[[l]])[2])
  for (i in 1:ncol(gene_data$N[[l]])){
    # the formulae used in STARmap protocol
    gene_data$N[[l]][,i] <- log(1 + median_expression[l]*((gene_data$data[[l]][,i] + 0.01)/sum(gene_data$data[[l]][,i])))
  }
  rownames(gene_data$N[[l]]) <- gene_data$name
}


head(gene_data$N[[1]], c(6,6))

genes <- list()

for(l in 1:length(gene_data$data)){
  counts <- gene_data$data[[l]]
  if(!require(Matrix))install.packages("Matrix"); suppressPackageStartupMessages(library(Matrix))
  counts <- as(counts, "sparseMatrix")
  colnames(counts) <- 1:ncol(counts)
  if(!require(Seurat))install.packages("Seurat"); suppressPackageStartupMessages(library(Seurat))
  
  meta_data <- data.frame(row = as.numeric(gene_data$location[[l]][,1]), 
                          col = as.numeric(gene_data$location[[l]][,2]),
                          annotation = 1:ncol(counts))
  
  row.names(meta_data) <- 1:ncol(counts)
  ## create Seurat object
  sce <- CreateSeuratObject(counts = counts, meta.data = meta_data)
  # standard log-normalization
  sce <- NormalizeData(sce, verbose = F)
  seu <- FindVariableFeatures(sce, nfeatures = 160, verbose = F)
  if(!require(DR.SC))install.packages("DR.SC"); suppressPackageStartupMessages(library(DR.SC))
  seus <- FindSVGs(seu, nfeatures = 160, verbose = F)
  genes[[l]] = seus@assays$RNA@var.features[1:50]
}

common.genes = intersect(intersect(intersect(genes[[1]], genes[[2]]), genes[[3]]), genes[[4]])


dataa <- list()
locs <- list()
nearest_neighbors <- list()

for(l in 1:length(gene_data$data)){
  dataa[[l]] <- gene_data$N[[l]][common.genes, ]
  if(!require(fields))install.packages("fields"); suppressPackageStartupMessages(library(fields)) 
  # Calculate the distance matrix 
  distance_matrix <- fields::rdist(gene_data$location[[l]])
  distance_matrix <- as.matrix(distance_matrix)
  colnames(distance_matrix) <- NULL; rownames(distance_matrix) = NULL
  distance_matrix_scaled = distance_matrix
  
  if(!require(NPVecchia))install.packages("NPVecchia"); suppressPackageStartupMessages(library(NPVecchia))
  
  order <- order_maximin_dist(distance_matrix_scaled)
  
  # Reorder data and location by maximin ordering
  dataa[[l]] <- dataa[[l]][, order]
  locs[[l]] <- gene_data$location[[l]][order, ]
  # Find the Euclidean neighbors
  distance_matrix_scaled_ordered = rdist(locs[[l]]) # Find the disatnce matrix after ordering
  distance_matrix_scaled_ordered = distance_matrix_scaled_ordered # Scale the matrix so that the distances are between 0 and 1
  n_locs = ncol(dataa[[l]]) # Number of spatial locations
  num_reps = nrow(dataa[[l]]) # Number of replicates
  
  nearest_neighbors[[l]] <- find_nn_dist(distance_matrix_scaled_ordered, n_locs)
}

counts <- list()
for(l in 1:length(gene_data$data)){
  counts[[l]] <- gene_data$data[[l]][common.genes,]  
}
saveRDS(counts, "~/BayesFMSRT/Real Data Application/counts.four.rds")
saveRDS(dataa, "~/BayesFMSRT/Real Data Application/dataa.four.rds")
saveRDS(locs, "~/BayesFMSRT/Real Data Application/locs.four.rds")
saveRDS(nearest_neighbors, "~/BayesFMSRT/Real Data Application/nearest_neighbors.four.rds")
