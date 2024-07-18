# This script defines auxiliary functions to be used in the simulation
# these functions catch errors and warnings:
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))
run_SAM <- quietly(safely(sam))
run_SEM <- quietly(safely(sem))
run_SEcorrection <- quietly(safely(stepwiseSE))

# function to approximate Nickell's bias:
nickells_bias <- function(phi, obs) {
  k <- ncol(phi) # Number of variables in the VAR model
  I_k <- diag(k) # Identity matrix of order k
  
  # Calculate the bias correction term
  bias <- -(I_k + phi) / (obs - 1)
  
  return(bias)
}