## Bootsrap residuals in regression model 
# @ param model: regressio model lm.function
# @ param iterations: number of iterations
##
# return: dataframe including simulated beta values

resBoot <- function(model,iterations){
  betas <- as.numeric(coefficients(model))
  residuals <- resid(model)
  n <- length(residuals)
  C <- model.matrix(model) # model in matrix form
  # simulated beta values
  boot_Betas <- matrix(nrow = iterations, ncol = length(betas))
  for(i in 1:iterations){ 
    # residual sample with replacement
    y_residualSample <- sample(residuals,n,replace = T)
    
    # Bootstrap responses
    y_star <- C %*% betas + y_residualSample
    
    # bootstrap least-square estimates
    beta_star <- solve(t(C) %*% C) %*% t(C) %*% y_star
    
    for(a in 1:length(betas)){
      boot_Betas[i,a] <- beta_star[a]
    }
  }
  return(boot_Betas)
}

values <- resBoot(lm_sepal,100)
values
mean(values[,1])
mean(values[,2])

## Bootsrap pairs in regression model 
# @ param model: regressio model lm.function
# @ param iterations: number of iterations
##
# return: dataframe including simulated beta values

pairsBoot <- function(model, iterations){
  betas <- as.numeric(coefficients(model))
  n <- length(model$residuals)
  i_set <- seq(1,n,1)
  boot_Betas <- matrix(nrow = iterations, ncol = length(betas))
  
  ## not sampled values
  y_values <- as.matrix(model.frame(model)[1])
  c_values <- as.matrix(model.frame(model)[,-1])
  
  # sampled values
  c_vec <- matrix(nrow = length(c_values), ncol = length(c_values[1]))
  y_vec <- matrix(nrow = n)
  # for the C matrix
  ones <- rep(1,n)
  
  for(a in 1:iterations){
    # for every iteration new sample
    i_bootSample <- sample(i_set,n,replace = T)
    for(i in 1:n){
      index <- i_bootSample[i]
      y_vec[i] <- y_values[index]
      c_vec[i,] <- c_values[index,]
    }
    C <- cbind(ones,c_vec)
    beta_star <- solve(t(C) %*% C) %*% t(C) %*% y_vec
    
    for(o in 1:length(betas)){
      boot_Betas[a,o] <- beta_star[o]
    }
    
  }
  return(boot_Betas)
}

