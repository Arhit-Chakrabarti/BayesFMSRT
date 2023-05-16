# Source all the required libraries
suppressPackageStartupMessages(library(MixMatrix))
suppressPackageStartupMessages(library(MCMCpack))
suppressPackageStartupMessages(library(NPVecchia))
suppressPackageStartupMessages(library(fields)) # Contains rdist function
suppressPackageStartupMessages(library(adaptMCMC))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mvtnorm))


# Quiet Function 
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

suppressPackageStartupMessages(source("functions_replicates.R"))
#######################################################
# READ THE DATA
#######################################################

dataa <- readRDS("~/BayesFMSRT/Real Data Application/dataa.four.rds")

n_locs = c(ncol(dataa[[1]]),
           ncol(dataa[[2]]),
           ncol(dataa[[3]]),
           ncol(dataa[[4]])) # Number of spatial locations

num_reps = nrow(dataa[[1]]) # Number of replicates

nearest_neighbors <- readRDS("~/BayesFMSRT/Real Data Application/nearest_neighbors.four.rds")

m = 10 # Number of nearest neighbors 
N = nrow(dataa[[1]]) # Define the number of replicates

nu = N # Degrees of freedom for the IW prior

Psi = diag(1, num_reps) # Set scaling parameter for IW prior
Lambda = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi))) # Starting value for the IW sampling (taken to be equal to true value)

n_iterations = 5000

# Define lists to store the samples
U_samples <- replicate(length(dataa), list(), simplify = FALSE)
d_samples <- replicate(length(dataa), list(), simplify = FALSE)

post <- list()

theta_samples <- matrix(NA, nrow = n_iterations, ncol = 3)

Lambda_samples <- list()

#A scale matrix similar to this worked well for our application
scale_mat <- matrix(c(0.05, -0.04, 0, -0.04, 0.05, 0, 0, 0, 0.01), nc = 3)

#number of samples
nruns <- 2
log_likelihood <- list()
init_theta = c(1,-1,0)

start_time <- Sys.time()
for(i in 1:n_iterations){
  # Printing the iterations
  if(i == 1){
    cat(paste0("Iteration: ", i, "\n"))
  }
  if(i %% floor((5/100)*(n_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", i, "\n"))
  }
  # Find the optimal thetas as a function of Lambda from the previous iteration
  thetas_new <- quiet(adaptMCMC::MCMC(minus_loglikeli_my_new, datumT = dataa, 
                                      NNarrayT = nearest_neighbors,
                                      Lambda = Lambda,
                                      m = m,
                                      init = init_theta, negativ = FALSE,
                                      scale = scale_mat, adapt = TRUE, 
                                      acc.rate = 0.234, n = nruns, 
                                      showProgressBar = FALSE)$samples[nruns, ])
  
  theta_samples[i, ] = thetas_new
  for(r in 1:length(dataa)){
    # Convert the hyper-parameters thetas to priors which is a function of Lambda from previous iteration
    priors = thetas_to_priors_my(thetas_new, n = nrow(nearest_neighbors[[r]]), m = m)
    # get initial posterior sample
    post[[r]] = get_posts_my(datum = dataa[[r]], priors = priors, Lambda = Lambda, NNarray = nearest_neighbors[[r]])
    
    # Sample posterior sparse matrix U as a function of Lambda from previous iteration
    U_samples[[r]][[i]] = samp_posts_my_new(post[[r]], NNarray = nearest_neighbors[[r]], bayesian = TRUE)
    # Store the diagonal elements of U
    d_samples[[r]][[i]] = diag(U_samples[[r]][[i]])
  }
  
  
  SS = matrix(0, nrow = N, ncol = N) # Define the sample sum of square matrix
  # ALTERNATE
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[r]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
        SS = SS + tcrossprod(dataa[[r]][, ind])/d
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        SS = SS + (tcrossprod(dataa[[r]][, ind] - mean.normal))/d
      }
    }
  }
  
  # Draw sample from IW i.e. from full conditional [Lambda|-]
  Lambda = riwish(v = nu + sum(n_locs), S = as.matrix(forceSymmetric(Psi + SS)))
  # Store Lambda as sample
  Lambda_samples[[i]] = as.matrix(forceSymmetric(Lambda))
  log_like <- replicate(length(dataa), list(), simplify = FALSE)
  
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[[r]]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        
      }
      
      log_like[[r]][[ind]] = dmvnorm(x = dataa[[r]][, ind], mean = mean.normal, sigma = d * Lambda_samples[[i]], log = TRUE)
    }
    log_likelihood[[i]] = sum(unlist(log_like))
  }
}# End of Gibbs Sampling
end_time <- Sys.time()
time.taken = end_time - start_time

U1_post = Matrix(0, nrow = n_locs[1], ncol = n_locs[1], sparse = TRUE)
U2_post = Matrix(0, nrow = n_locs[2], ncol = n_locs[2], sparse = TRUE)
U3_post = Matrix(0, nrow = n_locs[3], ncol = n_locs[3], sparse = TRUE)
U4_post = Matrix(0, nrow = n_locs[4], ncol = n_locs[4], sparse = TRUE)

Lambda_post = matrix(0, nrow = num_reps, ncol = num_reps)
theta_post = 0

burn = 25 # Number of burn-in

thin = 1

samples <- seq(from = (burn + 1), to = n_iterations, by = thin)

for(i in samples){
  U1_post = U1_post + U_samples[[1]][[i]]  # Sum of all matrices from Posterior of U1
  U2_post = U2_post + U_samples[[2]][[i]]  # Sum of all matrices from Posterior of U2
  U3_post = U3_post + U_samples[[3]][[i]]  # Sum of all matrices from Posterior of U3
  U4_post = U4_post + U_samples[[4]][[i]]  # Sum of all matrices from Posterior of U4
  Lambda_post = Lambda_post + Lambda_samples[[i]] # Sum of all matrices from Posterior of Lambda
}

U1_post = U1_post/length(samples) # Posterior Mean as an estimate of U1
U2_post = U2_post/length(samples) # Posterior Mean as an estimate of U2
U3_post = U3_post/length(samples) # Posterior Mean as an estimate of U3
U4_post = U4_post/length(samples) # Posterior Mean as an estimate of U4


Sigma1_post = solve(t(U1_post),  solve(U1_post)) # Posterior estimate of Sigma1
Sigma2_post = solve(t(U2_post),  solve(U2_post)) # Posterior estimate of Sigma2
Sigma3_post = solve(t(U3_post),  solve(U3_post)) # Posterior estimate of Sigma3
Sigma4_post = solve(t(U4_post),  solve(U4_post)) # Posterior estimate of Sigma3

Lambda_post = Lambda_post/length(samples) # Posterior Mean as an estimate of Lambda

Estimates_Samples1 <- list(Sigma1_post = Sigma1_post,
                           Sigma2_post = Sigma2_post,
                           Sigma3_post = Sigma3_post,
                           Sigma4_post = Sigma4_post,
                           Lambda_post = Lambda_post,
                           log_likelihood_samples = log_likelihood,
                           time.taken = time.taken,
                           n_iterations = n_iterations,
                           burn = burn,
                           thin = thin)

Estimates_Samples2 <- list(U1_samples = U_samples[[1]],
                           U2_samples = U_samples[[2]],
                           U3_samples = U_samples[[3]],
                           U4_samples = U_samples[[4]],
                           Lambda_samples = Lambda_samples,
                           theta_samples = theta_samples)

saveRDS(Estimates_Samples1, "~/BayesFMSRT/Real Data Application/Estimates_Samples1.rds")
saveRDS(Estimates_Samples2, "~/BayesFMSRT/Real Data Application/Estimates_Samples2.rds")

