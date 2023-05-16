suppressPackageStartupMessages(library(Matrix))
dataa <- readRDS("~/BayesFMSRT/Real Data Application/dataa.four.rds")
Estimates1 <- readRDS("~/BayesFMSRT/Real Data Application/Estimates_Samples1.rds")
Estimates2 <- readRDS("~/BayesFMSRT/Real Data Application/Estimates_Samples2.rds")
U1_samples <- Estimates2$U1_samples
U2_samples <- Estimates2$U2_samples
U3_samples <- Estimates2$U3_samples
U4_samples <- Estimates2$U4_samples
Lambda_samples <- Estimates2$Lambda_samples
n_iterations <-length(Estimates1$log_likelihood_samples)

dmatrix_norm <- function(X, U, V){
  n = nrow(X); p = ncol(X)
  
  first = - 0.5 * sum(diag(solve(V, t(X)) %*% solve(U, X)))
  second = -0.5 * n * p * log(2*pi)
  third = - 0.5 * n * determinant(V, logarithm = TRUE)$modulus[1]
  forth = - 0.5 * p * determinant(U, logarithm = TRUE)$modulus[1]
  return(first + second + third + forth)
} 

log_like_new <- list()

for(i in 1:n_iterations){
  if(i == (1)){
    cat(paste0("Iteration: ", i, "\n"))
  }
  if(i %% floor((5/100)*(n_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", i, "\n"))
  }
  log_like_new[[i]] <- dmatrix_norm(X = dataa[[1]], U = cov2cor(Lambda_samples[[i]]), V = cov2cor(solve(t(U1_samples[[i]]),  solve(U1_samples[[i]])))) +
    dmatrix_norm(X = dataa[[2]], U = cov2cor(Lambda_samples[[i]]), V = cov2cor(solve(t(U2_samples[[i]]),  solve(U2_samples[[i]])))) +
    dmatrix_norm(X = dataa[[3]], U = cov2cor(Lambda_samples[[i]]), V = cov2cor(solve(t(U3_samples[[i]]),  solve(U3_samples[[i]])))) +
    dmatrix_norm(X = dataa[[4]], U = cov2cor(Lambda_samples[[i]]), V = cov2cor(solve(t(U4_samples[[i]]),  solve(U4_samples[[i]]))))
}

saveRDS(log_like_new, "~/BayesFMSRT/Real Data Application/log_like_correlation.rds")
