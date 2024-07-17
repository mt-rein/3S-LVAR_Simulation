## which (datagenerating) values for theta lead to a (datagenerating) value for rho that is in line with the levels of the simulation?
combinations = expand.grid(phi_size = c("small", "large"),
                           rho_gen = c("small", "medium", "large", "very large"))

output <- data.frame(phi_size = character(0),
                     rho_gen = character(0),
                     theta1 = numeric(0),
                     theta2 = numeric(0),
                     theta3 = numeric(0),
                     theta4 = numeric(0),
                     theta5 = numeric(0),
                     theta6 = numeric(0),
                     theta7 = numeric(0),
                     theta8 = numeric(0),
                     rho1 = numeric(0),
                     rho2 = numeric(0))

for(row in 1:nrow(combinations)){
  phi_size = combinations$phi_size[row] |> as.character()
  rho_gen = combinations$rho_gen[row] |> as.character()
  if(phi_size == "small"){
    phi11_pop <- .25                                                            # effect of f1 on f1
    phi12_pop <- .125                                                           # effect of f2 on f1
    phi21_pop <- .125                                                           # effect of f1 on f2
    phi22_pop <- .25                                                            # effect of f2 on f2
    
  }
  
  if(phi_size == "large"){
    phi11_pop <- .5                                                             # effect of f1 on f1
    phi12_pop <- .25                                                            # effect of f2 on f1
    phi21_pop <- .25                                                            # effect of f1 on f2
    phi22_pop <- .5                                                             # effect of f2 on f2
  }
  phimat <- matrix(c(phi11_pop, phi21_pop, phi12_pop, phi22_pop), ncol = 2)     # combine into matrix
  zetamat <- matrix(c(1.5, .5, .5, 1.5), ncol = 2)
  
  loadings <- rep(1, 4)
  lambda <- matrix(c(loadings, rep(0, 4), rep(0, 4), loadings),
                   nrow = 8, ncol = 2)
  
  # since the total factor variance (psi) differs between conditions (low or large regression
  # parameters), we need to obtain the population psi and then obtain error 
  # variances that lead to a certain rho value
  # get population psi (factor variance)
  psimat <- solve(diag(2*2) - kronecker(phimat, phimat)) %*% c(zetamat) |> matrix(nrow = 2)
  
  # write an objective function that gives the difference between the computed rho
  # (based on the error variances) and desired rho
  rho_difference <- function(errorvars, psi, lambda, target_rho) {
    theta <- matrix(0, nrow = 8, ncol = 8)                                      # error variance matrix with 0 on off-diagonal
    diag(theta) <- errorvars                                                    # add error variances to the diagonal of theta
    computed_rho <- diag(psi %*% t(lambda) %*% solve(lambda %*% psi %*% t(lambda) + theta) %*% lambda)
    diff_rho <- sum((computed_rho - target_rho)^2)                              # squared difference between computed and desired rho as the objective
    return(diff_rho)
  }
  
  # initial values for error variances
  initial_errorvars <- rep(0.5, 8)
  
  # target values for rho
  if(rho_gen == "small"){
    target_rho <- c(0.5, 0.5)
  }
  if(rho_gen == "medium"){
    target_rho <- c(0.7, 0.7)
  }
  if(rho_gen == "large"){
    target_rho <- c(0.9, 0.9)
  }
  if(rho_gen == "very large"){
    target_rho <- c(0.99, 0.99)
  }
  
  
  # run optimization
  result <- optim(
    par = initial_errorvars,
    fn = rho_difference,
    psi = psimat,  # Provide psi, lambda, and target_rho as additional arguments
    lambda = lambda,
    target_rho = target_rho,
    method = "L-BFGS-B",
    lower = rep(0.0001, 8)
  )
  
  # retrieve the optimized theta
  optimized_errorvars <- result$par
  
  
  # compute rho based on these parameters to double check that the values are correct
  psi <- psimat
  theta <- matrix(0, nrow = 8, ncol = 8) # create error covariance matrix
  diag(theta) <- optimized_errorvars
  sigma <- lambda %*% psi %*% t(lambda) + theta
  
  rho <- diag(psi %*% t(lambda) %*% solve(sigma) %*% lambda)
  
  
  output[row, ] = c(phi_size, rho_gen, optimized_errorvars, rho)
}

output
