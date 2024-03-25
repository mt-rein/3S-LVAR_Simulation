## this function generates data for a single individual according to a vector autoregressive model
sim_VAR <- function(factors, obs, phi, zeta, mu, burn_in = 0){
  # factors = number of factors
  # obs = number of observations
  # phi = auto-regressive effect (a matrix in case of multiple constructs)
  # zeta = innovation variance (a matrix in case of multiple constructs)
  # mu = the mean of the time series (a vactor in case of multiple constructs)
  # burn_in = length of burn in (i.e., data that are generated to remove influence of initial random draw)
  
  
  # create empty dataframe of length obs + burn_in
  data <- as.data.frame(matrix(NA, nrow = burn_in + obs, ncol = factors))
  names(data) <- paste0("eta", 1:factors)
  intercept <- solve(solve(diag(factors) - phi), mu)
  
  for(i in 1:nrow(data)){
    # simulate the first observation from the person's starting point (= intercept)
    if(i == 1){
      data[i,] <- mu
    }
    
    # then loop through all the rows, predict the current observation from the previous observations, then add random innovation
    if(i > 1){
      predicted <- intercept + phi %*% unlist(data[i-1,])
      data[i, ] <- MASS::mvrnorm(n = 1, mu = predicted, Sigma = zeta)
    }
  }
  
  # remove the first rows, depending on length of burn in
  if(burn_in > 0){
    data <- dplyr::slice(data, burn_in+1:n()) 
  }
  
  data$obs <- 1:nrow(data)
  
  return(data)
}