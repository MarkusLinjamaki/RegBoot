library(boot)
library(bootstrap)
library(readxl)
library(ggplot2)
## Bootsrap residuals in regression model 
# @ param model: regressio model lm.function
# @ param iterations: number of iterations
##
# return: dataframe including simulated beta values
resBoot <- function(model,iterations){
  beta_0 <- rep(0,times = iterations)
  beta_1 <- rep(0,times = iterations)
  betas <- as.numeric(coefficients(model))
  residuals <- resid(lm_water)
  n <- length(residuals)
  C <- model.matrix(lm_water) # model in matrix form
  
  
  for(i in 1:iterations){ 
    # sampling w/ replacement from empirical residual dist
    y_residualSample <- sample(residuals,n,replace = T)
    y_star <- C %*% betas + y_residualSample
    beta_star <- solve(t(C) %*% C) %*% t(C) %*% y_star
    beta_0[i] <- beta_star[1]
    beta_1[i] <- beta_star[2]
  }
  betaStars <- data.frame(beta_0,beta_1)
}

## Bootsrap pairs in regression model 
# @ param model: regressio model lm.function
# @ param iterations: number of iterations
##
# return: dataframe including simulated beta values
pairsBoot <- function(model, iterations){
  beta_0 <- rep(0,times = iterations)
  beta_1 <- rep(0,times = iterations)
  betas <- as.numeric(coefficients(model))
  
  # resampled {1,2,...,n}
  n <- length(model$residuals)
  i_set <- seq(1,n,1)
  ## not sampled values
  y_values <- as.matrix(model.frame(lm_water)[1])
  c_values <- as.matrix(model.frame(lm_water)[2])
  # sampled values
  c_vec <- rep(0,times = n)
  y_vec <- rep(0,times = n)
  # for the C matrix
  ones <- rep(1,n)

  for(a in 1:iterations){
    # for every iteration new sample
    i_bootSample <- sample(i_set,n,replace = T)
    for(i in 1:n){
      index <- i_bootSample[i]
      c_vec[i] <- c_values[index]
      y_vec[i] <- y_values[index]
    }
    C <- cbind(ones,c_vec)
    beta_star <- solve(t(C) %*% C) %*% t(C) %*% y_vec
    beta_0[a] <-beta_star[1]
    beta_1[a] <-beta_star[2]
  }
  betaStars <- data.frame(beta_0,beta_1)
}

