Outputs = readRDS("~/BayesFMSRT/Real Data Application/Outputs.rds")
LL = readRDS("~/BayesFMSRT/Real Data Application/Real Data Application/LL.rds")

burn = 2500; n_iterations = 5000; thin = 1
samples <- seq((burn + 1), n_iterations, by = thin)

library(tidyverse)
log_like <- unlist(LL[samples])
data.frame(x = 1:length(log_like), y = log_like) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood", x = "Iterations post burn-in", y = "") +
  theme(axis.text=element_text(size=12),                                                                          
        axis.title=element_text(size=12),                                                                  
        plot.title = element_text(size=12))


library(forecast)
ggAcf(x = log_like, lag.max = 40) + ggtitle("ACF of log-likelihood") + labs(y = "") + ylim(c(-0.05, 0.25)) + 
  theme(axis.text=element_text(size=12),                                                                           
        axis.title=element_text(size=12), 
        plot.title = element_text(size=12))
################################################################################
### Estimated Gene Network
################################################################################
Lambda.Precision = cov2cor(solve(cov2cor(Outputs$Lambda_post)))
Lambda.Precision.abs = abs(Lambda.Precision)
Lambda.Precision.abs[!lower.tri(Lambda.Precision.abs, diag = FALSE)] <- 0
cell.connection = apply(Lambda.Precision.abs, 1, function(x){which(x > 0.1, arr.ind = TRUE)})

connections <- list()
for(i in 1:length(cell.connection)){
  if(! identical(cell.connection[[i]], integer(0))){
    connections[[i]] <- cbind(cell.connection[[i]], i)
  }
  
}
connections <- Filter(Negate(is.null), connections)
network = NULL
for(i in 1:length(connections)){
  network = rbind(network, connections[[i]])
}

colnames(network) <- c("gene number1", "gene number2")

# Read in the original data to get the names of the genes
dataa <- readRDS("~/BayesFMSRT/Real Data Application/dataa.four.rds")
# Get the genes that form a network and those who are isolated 
gene.numbers = as.vector(t(network))
gene.network = rownames(dataa[[2]])[as.vector(t(network))]
gene.no.network = rownames(dataa[[2]])[setdiff(x = 1:nrow(dataa[[2]]), y = unique(gene.numbers))]

library(igraph)
# Create the Graph
g <- graph(edges = gene.network, isolates = gene.no.network, directed = FALSE)
set.seed(234)
# Plot the Graph
par(mar=c(0,0,0,0) + 0.1)
plot(g, vertex.size = 20, layout = layout.fruchterman.reingold, vertex.label.cex = 1.5, vertex.shape = "none", edge.color="gray8")

################################################################################
# Heatmap of row precision to understand gene-coexpression
################################################################################
# Function to get the lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Function to get the  upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# Rownames and column names of the precision matrix are taken to be names of the genes
rownames(Lambda.Precision) <- rownames(dataa[[1]])
colnames(Lambda.Precision) <- rownames(dataa[[1]])
# Convert Precision matrix to upper triangularized version
upper_tri <- get_upper_tri(Lambda.Precision)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Plot the heatmap
library(corrplot)
corrplot(-Lambda.Precision, method = 'color', order = 'hclust', type = 'lower', diag = FALSE, tl.cex = 1.5)

################################################################################
### Spectral clustering
################################################################################
library(Matrix)
# Change the Spatial Correslation matrices for each of the four samples
S <- as.matrix(cov2cor(Outputs$Sigma1_post))
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
V = eigL$vectors[, (n_locs-3):(n_locs-1) ]
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
# Read the locations
locs <- readRDS("~/BayesFMSRT/Real Data Application/locs.four.rds")
library(DR.SC)
# Read the Counts data
counts <- readRDS("~/BayesFMSRT/Real Data Application/counts.four.rds")
# For every sample number 1 to 4 change the numbers within the list() accordingly
library(Matrix)
counts.c <- as(counts[[1]], "sparseMatrix")
library(Seurat)
colnames(counts.c) = 1:ncol(counts.c)
# Create metadata for plotting using DR.SC package
meta_data <- data.frame(row = as.numeric(locs[[1]][,1]), 
                        col = as.numeric(locs[[1]][,2]),
                        annotation = colnames(counts.c))

row.names(meta_data) <- colnames(counts.c)
## create Seurat object
sce <- CreateSeuratObject(counts=counts.c, meta.data = meta_data)
# standard log-normalization
sce <- NormalizeData(sce, verbose = F)
# standard variable features in DR.SC (This is done for all the features so does not create any changes)
seu <- FindVariableFeatures(sce, nfeatures = nrow(counts.c), verbose = F)
# standard spatially variable features in DR.SC (This is done for all the features so does not create any changes)
seus <- FindSVGs(seu, nfeatures = nrow(counts.c))

# FROM WSS PLOT SELECT THE NUMBER OF CLUSTERS
# Perform K means clustering
K.means = kmeans(V, centers = 5)
# Perform Spatial clustering
seus <- DR.SC(seus, K = 5, platform = 'Visium', verbose=T)
# Plot the DR.SC Clustering
plot.DRSC <- spatialPlotClusters(seus)
# Manual colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 10
cols = gg_color_hue(n)

myvalues = c("1" = "#F8766D", 
             "4" = "#A3A500",
             "5" = "#00BF7D", 
             "2" = "#00B0F6", 
             "3" = "#E76BF3")

# Plot the Clustering from our method
plot.Bayes <- data.frame(x = locs[[1]][, 1], y = locs[[1]][, 2], cluster = factor(K.means$cluster)) %>% ggplot(aes(x = x, y = y, col = cluster)) + geom_point(alpha = 0.7, size = 3)
# + labs(title = "Spectral clustering from Posterior estimate")
plot.Bayesian <- plot.Bayes  +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold", size=18),
        legend.text = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(title="clusters",
                               title.hjust = 0.5)) + 
  scale_color_manual(values = myvalues)

# Add same colors to DR.SC plot for comparison
plot.DR.SC <- plot.DRSC + scale_color_manual(values = myvalues)
library(gridExtra)
plot.Bayesian # Plot from our method
plot.DR.SC # Plot from DR.SC


################################################################################################
# COMAPRE CLUSTERING PERFORMANCE
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- data.frame(t(dataa[[1]]), cluster_method_1 = K.means$cluster, cluster_method_2 = seus$spatial.drsc.cluster)
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))
# Function to compare clustering performance by looking at Total Within Sum of Squares
dataa.aug %>%
  gather(method, cluster, cluster_method_1, cluster_method_2) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)
