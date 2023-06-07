################################################################################
# Create the data from the STARMap dataset
################################################################################
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
  data <- read.csv(paste("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/BayesFMSRT/STARmap/", data_list[i],"/",data_list[i],".count.csv", sep = ""), header = TRUE)
  gene_data$data[[i]] <- t(data)
  ncells[i] <- nrow(data)
  loc <- read.csv(paste("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/BayesFMSRT/STARmap/", data_list[i],"/",data_list[i],".loc.csv", sep = ""), header = TRUE)
  gene_data$location[[i]] <- loc
}

gene_data$name <- read.csv(paste("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/BayesFMSRT/STARmap/", data_list[1],"/",data_list[1],".gene_symbol.csv", sep = ""), header = TRUE)[,1]


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

for(l in 1:length(gene_data$data)){
  dataa[[l]] <- gene_data$N[[l]][common.genes, ]
  locs[[l]] <- gene_data$location[[l]]
  
}

################################################################################
# Run NPVecchia
################################################################################
library(NPVecchia)
#############################################
# Sample 1
#############################################
ans <- run_npvecchia(dataa[[1]], locs[[1]], corr_order = FALSE)
Sigma1_MAP = solve(t(ans$u), solve(ans$u))
locs[[1]] = locs[[1]][ans$order, ]

#############################################
# Sample 2
#############################################
ans <- run_npvecchia(dataa[[2]], locs[[2]], corr_order = FALSE)
Sigma2_MAP = solve(t(ans$u), solve(ans$u))
locs[[2]] = locs[[2]][ans$order, ]
# #############################################
# # Sample 2 HAVE TO DO MANUALLY
# #############################################
# #Compute maximin ordering
# order <- orderMaxMinFaster(locs[[2]])
# #Reorder data and location by correlation ordering
# dataa[[2]] <- dataa[[2]][, order]
# locs[[2]] <- locs[[2]][order, ]
# #Euclidean neighbors
# n_locs = nrow(locs[[2]])
# nearest_neighbors <- find_nn_dist(fields::rdist(locs[[2]]), n_locs)
# #initial thetas
# init_theta = c(1,-1,0)
# #Optimization, where one can change limits to match ones desired decay
# thetas_f <- optim(init_theta, minus_loglikeli_c, datum = dataa[[2]],
#                   NNarray = nearest_neighbors, method="L-BFGS-B",
#                   lower=-6, upper=4)
# 
# #Get MAP
# uhat <- get_map(thetas_f$par, dataa[[2]], nearest_neighbors)
# Sigma2_MAP = solve(t(uhat), solve(uhat))
# 

#############################################
# Sample 3
#############################################
ans <- run_npvecchia(dataa[[3]], locs[[3]], corr_order = FALSE)
Sigma3_MAP = solve(t(ans$u), solve(ans$u))
locs[[3]] = locs[[3]][ans$order, ]
#############################################
# Sample 4
#############################################
ans <- run_npvecchia(dataa[[4]], locs[[4]], corr_order = FALSE)
Sigma4_MAP = solve(t(ans$u), solve(ans$u))
locs[[4]] = locs[[4]][ans$order, ]
###############################################################################
################################################################################
### Spectral clustering (SAMPLE 1)
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- as.matrix(cov2cor(as.matrix(Sigma1_MAP)))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)
library(tidyverse)
# Plot of eigen values for Spectral cliustering
data.frame(Index = 1:10, Eigen = sort(eigL$values[(n_locs-10):(n_locs-1)])) %>%
  ggplot(aes(x = Index, y = Eigen)) + geom_point() + geom_line() + labs(y = "Eigen values", x = "") + 
  scale_x_discrete(limits = as.character(1:10)) + theme(axis.text=element_text(size=18),
                                                        axis.title=element_text(size=18))

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-4):(n_locs-1) ]
# To find the number of cluster
WSS <- 0
cent = 2:20
# Plot the Within Sum of Squares for different number of clusters
for(i in 1:length(cent)) WSS[i] <- kmeans(V, centers = cent[i])$tot.withinss
library(tidyverse)
data.frame(x = cent, y = WSS) %>% ggplot(aes(x = x, y = y)) + geom_point() + geom_line() + labs(x = "Number of clusters", y = "Total Within Sum of Squares") + theme(axis.text=element_text(size=18),                                                                           axis.title=element_text(size=18))

################################################################################################
# CLUSTERING PLOTS
################################################################################################
# FROM WSS PLOT SELECT THE NUMBER OF CLUSTERS
# Perform K means clustering
K.means = kmeans(V, centers = 5)


# Plot the Clustering from our method
plot.MAP <- data.frame(x = locs[[1]][, 1], y = locs[[1]][, 2], cluster = factor(K.means$cluster)) %>% ggplot(aes(x = x, y = y, col = cluster)) + geom_point(alpha = 0.7, size = 3)

plot.MAP
################################################################################################
# COMAPRE CLUSTERING PERFORMANCE 
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- data.frame(t(dataa[[1]]), cluster_method_1 = K.means$cluster)
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))
# Function to compare clustering performance by looking at Total Within Sum of Squares
dataa.aug %>%
  gather(method, cluster, cluster_method_1) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)



################################################################################
### Spectral clustering (SAMPLE 2)
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- as.matrix(cov2cor(as.matrix(Sigma2_MAP)))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)
library(tidyverse)
# Plot of eigen values for Spectral cliustering
data.frame(Index = 1:10, Eigen = sort(eigL$values[(n_locs-10):(n_locs-1)])) %>%
  ggplot(aes(x = Index, y = Eigen)) + geom_point() + geom_line() + labs(y = "Eigen values", x = "") + 
  scale_x_discrete(limits = as.character(1:10)) + theme(axis.text=element_text(size=18),
                                                        axis.title=element_text(size=18))

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-4):(n_locs-1) ]
# To find the number of cluster
WSS <- 0
cent = 2:20
# Plot the Within Sum of Squares for different number of clusters
for(i in 1:length(cent)) WSS[i] <- kmeans(V, centers = cent[i])$tot.withinss
library(tidyverse)
data.frame(x = cent, y = WSS) %>% ggplot(aes(x = x, y = y)) + geom_point() + geom_line() + labs(x = "Number of clusters", y = "Total Within Sum of Squares") + theme(axis.text=element_text(size=18),                                                                           axis.title=element_text(size=18))

################################################################################################
# CLUSTERING PLOTS
################################################################################################
# FROM WSS PLOT SELECT THE NUMBER OF CLUSTERS
# Perform K means clustering
K.means = kmeans(V, centers = 5)


# Plot the Clustering from our method
plot.MAP <- data.frame(x = locs[[2]][, 1], y = locs[[2]][, 2], cluster = factor(K.means$cluster)) %>% ggplot(aes(x = x, y = y, col = cluster)) + geom_point(alpha = 0.7, size = 3)

plot.MAP
################################################################################################
# COMAPRE CLUSTERING PERFORMANCE 
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- data.frame(t(dataa[[2]]), cluster_method_1 = K.means$cluster)
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))
# Function to compare clustering performance by looking at Total Within Sum of Squares
dataa.aug %>%
  gather(method, cluster, cluster_method_1) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)



################################################################################
### Spectral clustering (SAMPLE 3)
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- as.matrix(cov2cor(as.matrix(Sigma3_MAP)))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)
library(tidyverse)
# Plot of eigen values for Spectral cliustering
data.frame(Index = 1:10, Eigen = sort(eigL$values[(n_locs-10):(n_locs-1)])) %>%
  ggplot(aes(x = Index, y = Eigen)) + geom_point() + geom_line() + labs(y = "Eigen values", x = "") + 
  scale_x_discrete(limits = as.character(1:10)) + theme(axis.text=element_text(size=18),
                                                        axis.title=element_text(size=18))

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-4):(n_locs-1) ]
# To find the number of cluster
WSS <- 0
cent = 2:20
# Plot the Within Sum of Squares for different number of clusters
for(i in 1:length(cent)) WSS[i] <- kmeans(V, centers = cent[i])$tot.withinss
library(tidyverse)
data.frame(x = cent, y = WSS) %>% ggplot(aes(x = x, y = y)) + geom_point() + geom_line() + labs(x = "Number of clusters", y = "Total Within Sum of Squares") + theme(axis.text=element_text(size=18),                                                                           axis.title=element_text(size=18))

################################################################################################
# CLUSTERING PLOTS
################################################################################################
# FROM WSS PLOT SELECT THE NUMBER OF CLUSTERS
# Perform K means clustering
K.means = kmeans(V, centers = 5)


# Plot the Clustering from our method
plot.MAP <- data.frame(x = locs[[3]][, 1], y = locs[[3]][, 2], cluster = factor(K.means$cluster)) %>% ggplot(aes(x = x, y = y, col = cluster)) + geom_point(alpha = 0.7, size = 3)

plot.MAP
################################################################################################
# COMAPRE CLUSTERING PERFORMANCE
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- data.frame(t(dataa[[3]]), cluster_method_1 = K.means$cluster)
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))
# Function to compare clustering performance by looking at Total Within Sum of Squares
dataa.aug %>%
  gather(method, cluster, cluster_method_1) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)


################################################################################
### Spectral clustering (SAMPLE 4)
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- as.matrix(cov2cor(as.matrix(Sigma4_MAP)))
# Convert correlation matrix to lie between 0 and 1
S <- (S + 1)/2
# Get the Normalized Graph Laplacian matrix
D = diag(rowSums(S))
L = diag(1, nrow = nrow(S)) - solve(D,S)
# Perform eigen analysis on the Normalized Graph Laplacian matrix
eigL = eigen(L)
n_locs <- ncol(S)
library(tidyverse)
# Plot of eigen values for Spectral cliustering
data.frame(Index = 1:10, Eigen = sort(eigL$values[(n_locs-10):(n_locs-1)])) %>%
  ggplot(aes(x = Index, y = Eigen)) + geom_point() + geom_line() + labs(y = "Eigen values", x = "") + 
  scale_x_discrete(limits = as.character(1:10)) + theme(axis.text=element_text(size=18),
                                                        axis.title=element_text(size=18))

# Extract the eigen vectors corresponding to the (k-1) smallest eigen values 
V = eigL$vectors[, (n_locs-4):(n_locs-1) ]
# To find the number of cluster
WSS <- 0
cent = 2:20
# Plot the Within Sum of Squares for different number of clusters
for(i in 1:length(cent)) WSS[i] <- kmeans(V, centers = cent[i])$tot.withinss
library(tidyverse)
data.frame(x = cent, y = WSS) %>% ggplot(aes(x = x, y = y)) + geom_point() + geom_line() + labs(x = "Number of clusters", y = "Total Within Sum of Squares") + theme(axis.text=element_text(size=18),                                                                           axis.title=element_text(size=18))

################################################################################################
# CLUSTERING PLOTS
################################################################################################
# FROM WSS PLOT SELECT THE NUMBER OF CLUSTERS
# Perform K means clustering
K.means = kmeans(V, centers = 5)


# Plot the Clustering from our method
plot.MAP <- data.frame(x = locs[[4]][, 1], y = locs[[4]][, 2], cluster = factor(K.means$cluster)) %>% ggplot(aes(x = x, y = y, col = cluster)) + geom_point(alpha = 0.7, size = 3)

plot.MAP
################################################################################################
# COMAPRE CLUSTERING PERFORMANCE
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- data.frame(t(dataa[[4]]), cluster_method_1 = K.means$cluster)
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))
# Function to compare clustering performance by looking at Total Within Sum of Squares
dataa.aug %>%
  gather(method, cluster, cluster_method_1) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)
