# This script defines auxiliary functions to be used in the simulation
# these functions catch errors and warnings:
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))
run_SAM <- quietly(safely(sam))
run_SEM <- quietly(safely(sem))
run_SEcorrection <- quietly(safely(stepwiseSE))

# function to approximate Nickell's bias:
nickells_bias <- function(phi_est, obs) {
  k <- ncol(phi_est) # Number of variables in the VAR model
  I_k <- diag(k) # Identity matrix of order k
  
  # Calculate the bias correction term
  unbiased_phi <- (phi_est*(obs-1)+I_k)/(obs-2)
  
  return(unbiased_phi)
}