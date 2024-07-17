## This script defines the simulation function
do_sim <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  phi_size <-  cond$phi_size[pos] |> as.character()
  n <- cond$n[pos]
  obs <- cond$obs[pos]
  rho_gen <- cond$rho_gen[pos] |> as.character()
  variance_means <- cond$variance_means[pos] |> as.character()
  variance_zeta <- cond$variance_zeta[pos] |> as.character()
  variance_phi <- cond$variance_phi[pos] |> as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  #### 1) set data generation parameters ####
  # AR/CR parameters:
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
  
  # innovation (co)variance matrix:
  zeta1_pop <- 1.5                                                              # innovation variance of x1
  zeta2_pop <- 1.5                                                              # innovation variance of x2
  zeta12_pop <- .5                                                              # innovation covariance between x1 and x2
  zetamat <- matrix(c(zeta1_pop, zeta12_pop, zeta12_pop, zeta2_pop), ncol = 2)
  
  # means on the latent constructs
  mu1 <- 5
  mu2 <- 5
  
  # start timer for data generation:
  start_datagen <- Sys.time()
  
  #### 2) generate factor scores ####
  # create empty data frame:
  eta <- tibble(id = numeric(),
                obs = numeric(),
                eta1 = numeric(),
                eta2 = numeric())
  
  # loop over all subjects:
  for(person in 1:n){
    ## generate factor scores for each individual:
    # create data-generating parameter values from population values:
    # regression parameters:
    if(variance_phi == "no"){
      phimat_i <- phimat
    }
    if(variance_phi == "yes"){
      # draw individual-specific effects from truncated normal distribution
      phi11_i <- rtruncnorm(1, a = phi11_pop - .05, b = phi11_pop + .05,
                            mean = phi11_pop, sd = .5)
      phi12_i <- rtruncnorm(1, a = phi12_pop - .05, b = phi12_pop + .05,
                            mean = phi12_pop, sd = .5)
      phi21_i <- rtruncnorm(1, a = phi21_pop - .05, b = phi21_pop + .05,
                            mean = phi21_pop, sd = .5)
      phi22_i <- rtruncnorm(1, a = phi22_pop - .05, b = phi22_pop + .05,
                            mean = phi22_pop, sd = .5)
      phimat_i = matrix(c(phi11_i, phi21_i, phi12_i, phi22_i), ncol = 2)
    }
    
    # innovation (co)variances:
    if(variance_zeta == "no"){
      zetamat_i <- zetamat
    }
    if(variance_zeta == "yes"){
      # draw individual-specific (co)variances from truncated normal distribution
      zeta1_i <- rtruncnorm(1, a = zeta1_pop - .5, b = zeta1_pop + .5,
                            mean = zeta1_pop, sd = .5)
      zeta2_i <- rtruncnorm(1, a = zeta2_pop - .5, b = zeta2_pop + .5,
                            mean = zeta2_pop, sd = .5)
      zeta12_i <- rtruncnorm(1, a = zeta12_pop - .25, b = zeta12_pop + .25,
                             mean = zeta12_pop, sd = .5)
      zetamat_i <- matrix(c(zeta1_i, zeta12_i, zeta12_i, zeta2_i), ncol = 2)
    }
    
    # means:
    if(variance_means == "no"){
      mu1_i <- mu1
      mu2_i <- mu2
    }
    if(variance_means == "yes"){
      # draw individual-specific means from truncated normal distribution
      mu1_i <- rtruncnorm(1, a = mu1 - 2, b = mu1 + 2,
                          mean = mu1, sd = 1)
      mu2_i <- rtruncnorm(1, a = mu2 - 2, b = mu2 + 2,
                          mean = mu2, sd = 1)
    }

    # generate the factor scores
    eta_i <- sim_VAR(factors = 2, obs = obs,
                     phi = phimat_i, zeta = zetamat_i,
                     mu = c(mu1_i, mu2_i),
                     burn_in = 20)
    
    # add id variable:
    eta_i$id <- person
    
    eta <- dplyr::full_join(eta, eta_i, by = join_by(id, obs, eta1, eta2))
    }
  
  
  #### 3) generate indicators ####
  # lambda matrix (loadings), all 1:
  loadings <- rep(1, 4)
  lambda <- matrix(c(loadings, rep(0, 4), rep(0, 4), loadings),
                   nrow = 8, ncol = 2)
  
  # set error variances (see script "rhos.R"):
  theta <- matrix(0, nrow = 8, ncol = 8)                                        # create error covariance matrix
  if(phi_size == "small" & rho_gen == "small"){
    errorvars <- 6.15
  }
  if(phi_size == "large" & rho_gen == "small"){
    errorvars <- 8.832
  }
  if(phi_size == "small" & rho_gen == "medium"){
    errorvars <- 2.546
  }
  if(phi_size == "large" & rho_gen == "medium"){
    errorvars <- 3.398
  }
  if(phi_size == "small" & rho_gen == "large"){
    errorvars <- .639
  }
  if(phi_size == "large" & rho_gen == "large"){
    errorvars <- 0.8
  }
  if(phi_size == "small" & rho_gen == "very large"){
    errorvars <- 0.057
  }
  if(phi_size == "large" & rho_gen == "very large"){
    errorvars <- 0.07
  }
  diag(theta) <- errorvars                                                      # place error variances on theta diagonal
  
  epsilon <- mvrnorm(nrow(eta), mu = rep(0, 8), Sigma = theta, empirical=T)     # generate errors
  
  # generate indicator scores
  # (intercepts left out because we set them to 0):
  data <- as.matrix(eta[, c("eta1", "eta2")]) %*% t(lambda) + epsilon |>
    as.data.frame()
  colnames(data) <- paste0("v", 1:8)
  data$id <- eta$id
  data$obs <- eta$obs
  
  # save duration of data generation
  t_datagen <- difftime(Sys.time(), start_datagen, unit = "s")
  
  #### 4) 3S-LVAR ####
  # within-person mean center the indicators to remove between-person differences
  if(variance_means == "yes"){
    data_cent <- center_within(data, vars = paste0("v", 1:8), id = "id")
  } else {
    data_cent <- data
  }
  
  # write measurement model
  model_step1 <- "
    f1 =~ v1 + v2 + v3 + v4
    f2 =~ v5 + v6 + v7 + v8
    "
  
  ## run step 1
  start_step1 <- Sys.time()                                                     # start timer for step 1
  LVAR_step1 <- run_step1(data = data_cent,
                          measurementmodel = model_step1,
                          id = "id")
  # extract error/warning messages (if applicable):
  LVAR_step1_warning <- ifelse(is_empty(LVAR_step1$warnings),
                                    FALSE, TRUE)
  LVAR_step1_warning_text <- ifelse(is_empty(LVAR_step1$warnings),
                                         "",
                                         paste(c(LVAR_step1$warnings),
                                               collapse = "; ")
                                    )
  LVAR_step1_error <- ifelse(is_empty(LVAR_step1$result$error),
                                  FALSE, TRUE)
  LVAR_step1_error_text <- ifelse(is_empty(LVAR_step1$result$error),
                                       "",
                                       paste(c(LVAR_step1$result$error),
                                             collapse = "; ")
                                  )
  if(is_empty(LVAR_step1$result$error)){                                        # only proceed if there is no error in step 1
    output_LVAR_step1 <- LVAR_step1$result$result
    t_step1 <- difftime(Sys.time(), start_step1, unit = "s")                    # save duration of step 1 (if it was successful)
    } else {
      t_step1 <- NA
    }
  
  ## run step 2 only if step 1 was successful:
  if(!LVAR_step1_error){
    start_step2 <- Sys.time()                                                   # start timer for step 2
    LVAR_step2 <- run_step2(step1output = output_LVAR_step1)
    # extract error/warning messages (if applicable):
    LVAR_step2_warning <- ifelse(is_empty(LVAR_step2$warnings),
                                 FALSE, TRUE)
    LVAR_step2_warning_text <- ifelse(is_empty(LVAR_step2$warnings),
                                      "",
                                      paste(c(LVAR_step2$warnings),
                                            collapse = "; ")
                                      )
    LVAR_step2_error <- ifelse(is_empty(LVAR_step2$result$error),
                               FALSE, TRUE)
    LVAR_step2_error_text <- ifelse(is_empty(LVAR_step2$result$error),
                                    "",
                                    paste(c(LVAR_step2$result$error),
                                          collapse = "; ")
                                    )
    if(is_empty(LVAR_step2$result$error)){                                      # only proceed if there is no error in step
      output_LVAR_step2 <- LVAR_step2$result$result
      t_step2 <- difftime(Sys.time(), start_step2, unit = "s")                  # save duration of step 2 (if it was successful)
      } else {
        t_step2 <- NA
      }
    } else {
      # if step 1 was not successful:
      t_step2 <- NA
      LVAR_step2_warning <- FALSE
      LVAR_step2_warning_text <- "step1 not successful"
      LVAR_step2_error <- FALSE
      LVAR_step2_error_text <- "step1 not successful"
    }
  
  ## run step 3 only if step 1 and step 2 were successful:
  if(!LVAR_step1_error & !LVAR_step2_error){
    start_step3 <- Sys.time()                                                   # start timer of step 3
    LVAR_step3 <- run_step3(step2output = output_LVAR_step2)
    # extract error/warning messages (if applicable):
    LVAR_step3_warning <- ifelse(is_empty(LVAR_step3$warnings),
                                 FALSE, TRUE)
    LVAR_step3_warning_text <- ifelse(is_empty(LVAR_step3$warnings),
                                      "",
                                      paste(c(LVAR_step3$warnings),
                                            collapse = "; ")
                                      )
    LVAR_step3_error <- ifelse(is_empty(LVAR_step3$result$error),
                               FALSE, TRUE)
    LVAR_step3_error_text <- ifelse(is_empty(LVAR_step3$result$error),
                                    "",
                                    paste(c(LVAR_step3$result$error),
                                          collapse = "; ")
                                    )
  } else {
    # if step 1 or step 2 were not successful, set all results to NA:
    t_step3 <- NA
    LVAR_step3_warning <- FALSE
    LVAR_step3_warning_text <- "step1 or step2 not successful"
    LVAR_step3_error <- FALSE
    LVAR_step3_error_text <- "step2 or step2 not successful"
    SEcorr_step3_warning <- FALSE
    SEcorr_step3_warning_text <- "step1 or step2 not successful"
    SEcorr_step3_error <- FALSE
    SEcorr_step3_error_text <- "step1 or step2 not successful"
  }
  
  # if all steps were successful, extract results:
  if(!LVAR_step1_error & !LVAR_step2_error & !LVAR_step3_error){
    output_LVAR_step3 <- LVAR_step3$result$result
    t_step3 <- difftime(Sys.time(), start_step3, unit = "s")                    # save duration of step 3
    
    ## extract results
    # regression parameters and innovation variances
    LVAR_parameters <- parTable(output_LVAR_step3$fit)
    LVAR_phi11 <- LVAR_parameters |> 
      dplyr::filter(lhs == "f1", op == "~", rhs == "f1_lag") |> 
      dplyr::select(est) |> as.numeric()
    LVAR_phi12 <- LVAR_parameters |> 
      dplyr::filter(lhs == "f1", op == "~", rhs == "f2_lag") |> 
      dplyr::select(est) |> as.numeric()
    LVAR_phi21 <- LVAR_parameters |> 
      dplyr::filter(lhs == "f2", op == "~", rhs == "f1_lag") |> 
      dplyr::select(est) |> as.numeric()
    LVAR_phi22 <- LVAR_parameters |> 
      dplyr::filter(lhs == "f2", op == "~", rhs == "f2_lag") |> 
      dplyr::select(est) |> as.numeric()

    LVAR_zeta1 <- LVAR_parameters |>
      dplyr::filter(lhs == "f1", op == "~~", rhs == "f1") |> 
      dplyr::select(est) |> as.numeric()
    LVAR_zeta2 <- LVAR_parameters |>
      dplyr::filter(lhs == "f2", op == "~~", rhs == "f2") |> 
      dplyr::select(est) |> as.numeric()
    LVAR_zeta12 <- LVAR_parameters |>
      dplyr::filter(lhs == "f1", op == "~~", rhs == "f2") |> 
      dplyr::select(est) |> as.numeric()
    
    # (uncorrected) standard errors
    LVAR_phi11_se <- LVAR_parameters |> 
      dplyr::filter(lhs == "f1", op == "~", rhs == "f1_lag") |> 
      dplyr::select(se) |> as.numeric()
    LVAR_phi12_se <- LVAR_parameters |> 
      dplyr::filter(lhs == "f1", op == "~", rhs == "f2_lag") |> 
      dplyr::select(se) |> as.numeric()
    LVAR_phi21_se <- LVAR_parameters |> 
      dplyr::filter(lhs == "f2", op == "~", rhs == "f1_lag") |> 
      dplyr::select(se) |> as.numeric()
    LVAR_phi22_se <- LVAR_parameters |> 
      dplyr::filter(lhs == "f2", op == "~", rhs == "f2_lag") |> 
      dplyr::select(se) |> as.numeric()
    
    LVAR_zeta1_se <- LVAR_parameters |>
      dplyr::filter(lhs == "f1", op == "~~", rhs == "f1") |> 
      dplyr::select(se) |> as.numeric()
    LVAR_zeta2_se <- LVAR_parameters |>
      dplyr::filter(lhs == "f2", op == "~~", rhs == "f2") |> 
      dplyr::select(se) |> as.numeric()
    LVAR_zeta12_se <- LVAR_parameters |>
      dplyr::filter(lhs == "f1", op == "~~", rhs == "f2") |> 
      dplyr::select(se) |> as.numeric()
    
    } else {
      # if any step was not successful, set all results to NA:
      t_step3 <- NA
      
      LVAR_phi11 <- NA
      LVAR_phi12 <- NA
      LVAR_phi21 <- NA
      LVAR_phi22 <- NA
      LVAR_phi11corr <- NA
      LVAR_phi22corr <- NA
      
      LVAR_zeta1 <- NA
      LVAR_zeta2 <- NA
      LVAR_zeta12 <- NA
      
      LVAR_phi11_se <- NA
      LVAR_phi12_se <- NA
      LVAR_phi21_se <- NA
      LVAR_phi22_se <- NA
      
      LVAR_zeta1_se <- NA
      LVAR_zeta2_se <- NA
      LVAR_zeta12_se <- NA
    }
  
  
  # if step3 was successful, also perform the standard error correction
  if(!LVAR_step1_error & !LVAR_step2_error & !LVAR_step3_error){
    start_SEcorr <- Sys.time()                                                  # start timer of SE correction
    SEcorr <- run_SEcorrection(output_LVAR_step2, output_LVAR_step3)
    # extract error/warning messages (if applicable):
    SEcorr_warning <- ifelse(is_empty(SEcorr$warnings),
                             FALSE, TRUE)
    SEcorr_warning_text <- ifelse(is_empty(SEcorr$warnings),
                                  "",
                                  paste(c(SEcorr$warnings),
                                        collapse = "; ")
    )
    SEcorr_error <- ifelse(is_empty(SEcorr$result$error),
                           FALSE, TRUE)
    SEcorr_error_text <- ifelse(is_empty(SEcorr$result$error),
                                "",
                                paste(c(SEcorr$result$error),
                                      collapse = "; ")
    )

    # if correction was successful, extract results
    if(is_empty(SEcorr$result$error)){
      output_SEcorr <- SEcorr$result$result
      t_SEcorr <- difftime(Sys.time(), start_SEcorr, unit = "s")    # save duration of SE correction

      # corrected SEs:
      LVAR_phi11_secorr <- output_SEcorr$SE["f1~f1_lag"] |>
        as.numeric()
      LVAR_phi12_secorr <- output_SEcorr$SE["f1~f2_lag"] |>
        as.numeric()
      LVAR_phi21_secorr <- output_SEcorr$SE["f2~f1_lag"] |>
        as.numeric()
      LVAR_phi22_secorr <- output_SEcorr$SE["f2~f2_lag"] |>
        as.numeric()

      LVAR_zeta1_secorr <- output_SEcorr$SE["f1~~f1"] |>
        as.numeric()
      LVAR_zeta2_secorr <- output_SEcorr$SE["f2~~f2"] |>
        as.numeric()
      LVAR_zeta12_secorr <- output_SEcorr$SE["f1~~f2"] |>
        as.numeric()

    } else {
      # if SEcorrection was not successful, set all corrected SEs to NA:
      t_SEcorr <- NA

      LVAR_phi11_secorr <- NA
      LVAR_phi12_secorr <- NA
      LVAR_phi21_secorr <- NA
      LVAR_phi22_secorr <- NA

      LVAR_zeta1_secorr <- NA
      LVAR_zeta2_secorr <- NA
      LVAR_zeta12_secorr <- NA
    }


  } else {
    # if step3 was not successful, set all corrected SEs to NA:
    t_SEcorr <- NA
    SEcorr_step3_warning <- FALSE
    SEcorr_step3_warning_text <- "step3 not successful"
    SEcorr_step3_error <- FALSE
    SEcorr_step3_error_text <- "step3 not successful"

    LVAR_phi11_secorr <- NA
    LVAR_phi12_secorr <- NA
    LVAR_phi21_secorr <- NA
    LVAR_phi22_secorr <- NA

    LVAR_zeta1_secorr <- NA
    LVAR_zeta2_secorr <- NA
    LVAR_zeta12_secorr <- NA
  }

  #### 5) Naive Factor Scores (NFS) ####
  factorscores <- output_LVAR_step2$data[, c("id", "f1", "f2")]
  # factorscores have already been computed above (step 1 and 2 of 3S-LVAR)
  # so we can just reuse that object
  start_NFS <- Sys.time()                                                       # start timer

  # create additional rows and lagged variables
  factorscores <- factorscores |>
    group_by(id) |>
    do(add_row(.)) |>
    ungroup() |>
    fill(id) |>
    group_by(id) |>
    mutate(f1lag = dplyr::lag(f1),
           f2lag = dplyr::lag(f2)) |>
    ungroup()


  model_NFS <- "
  f1 ~ phi_f1_f1_lag*f1lag + phi_f1_f2_lag*f2lag
  f2 ~ phi_f2_f1_lag*f1lag + phi_f2_f2_lag*f2lag
  f1 ~~ zeta_f1_f1*f1
  f2 ~~ zeta_f2_f2*f2
  f1 ~~ zeta_f1_f2*f2
  "
  NFS <- run_SEM(model_NFS, data = factorscores,
                 missing = "ML",
                 cluster = "id")
  # extract error/warning messages (if applicable):
  NFS_warning <- ifelse(is_empty(NFS$warnings),
                        FALSE, TRUE)
  NFS_warning_text <- ifelse(is_empty(NFS$warnings),
                             "",
                             paste(c(NFS$warnings),
                                   collapse = "; ")
  )
  NFS_error <- ifelse(is_empty(NFS$result$error),
                      FALSE, TRUE)
  NFS_error_text <- ifelse(is_empty(NFS$result$error),
                           "",
                           paste(c(NFS$result$error),
                                 collapse = "; ")
  )
  # if no error during NFS estimation, get results
  if(is_empty(NFS$result$error)){
    output_NFS <- NFS$result$result
    t_NFS <- difftime(Sys.time(), start_NFS, unit = "s")                        # save duration of NFS

    ## extract results:
    # regression parameters and innovation variances
    NFS_parameters <- parTable(output_NFS)
    NFS_phi11 <- NFS_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    NFS_phi12 <- NFS_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(est) |> as.numeric()
    NFS_phi21 <- NFS_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    NFS_phi22 <- NFS_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(est) |> as.numeric()

    NFS_zeta1 <- NFS_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(est) |> as.numeric()
    NFS_zeta2 <- NFS_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(est) |> as.numeric()
    NFS_zeta12 <- NFS_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(est) |> as.numeric()

    # standard errors
    NFS_phi11_se <- NFS_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    NFS_phi12_se <- NFS_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(se) |> as.numeric()
    NFS_phi21_se <- NFS_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    NFS_phi22_se <- NFS_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(se) |> as.numeric()

    NFS_zeta1_se <- NFS_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(se) |> as.numeric()
    NFS_zeta2_se <- NFS_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(se) |> as.numeric()
    NFS_zeta12_se <- NFS_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(se) |> as.numeric()
  } else {
    # if NFS was not successful, set all results to NA
    t_NFS <- NA

    NFS_phi11 <- NA
    NFS_phi12 <- NA
    NFS_phi21 <- NA
    NFS_phi22 <- NA

    NFS_zeta1 <- NA
    NFS_zeta2 <- NA
    NFS_zeta12 <- NA

    NFS_phi11_se <- NA
    NFS_phi12_se <- NA
    NFS_phi21_se <- NA
    NFS_phi22_se <- NA

    NFS_zeta1_se <- NA
    NFS_zeta2_se <- NA
    NFS_zeta12_se <- NA
  }



  #### 6) SAM ####
  start_SAM <- Sys.time()                                                       # start timer for SAM

  data_SEM <- data_cent |>
    group_by(id) |>
    do(add_row(.)) |>
    ungroup() |>
    fill(id) |>
    group_by(id) |>
    mutate(v1lag = dplyr::lag(v1),
           v2lag = dplyr::lag(v2),
           v3lag = dplyr::lag(v3),
           v4lag = dplyr::lag(v4),
           v5lag = dplyr::lag(v5),
           v6lag = dplyr::lag(v6),
           v7lag = dplyr::lag(v7),
           v8lag = dplyr::lag(v8)) |>
    ungroup()


  model_SEM <- "
  f1lag =~ v1lag + v2lag + v3lag + v4lag
  f2lag =~ v5lag + v6lag + v7lag + v8lag
  f1 =~ v1 + v2 + v3 + v4
  f2 =~ v5 + v6 + v7 + v8

  f1 ~ phi_f1_f1_lag*f1lag + phi_f1_f2_lag*f2lag
  f2 ~ phi_f2_f1_lag*f1lag + phi_f2_f2_lag*f2lag
  f1 ~~ zeta_f1_f1*f1
  f2 ~~ zeta_f2_f2*f2
  f1 ~~ zeta_f1_f2*f2
  "

  SAM <- run_SAM(model_SEM, data = data_SEM, missing = "ML",
                 mm.args = list(bounds = "none"),
                 sam.method = "local")
  # extract error/warning messages (if applicable):
  SAM_warning <- ifelse(is_empty(SAM$warnings),
                        FALSE, TRUE)
  SAM_warning_text <- ifelse(is_empty(SAM$warnings),
                             "",
                             paste(c(SAM$warnings),
                                   collapse = "; ")
  )
  SAM_error <- ifelse(is_empty(SAM$result$error),
                      FALSE, TRUE)
  SAM_error_text <- ifelse(is_empty(SAM$result$error),
                           "",
                           paste(c(SAM$result$error),
                                 collapse = "; ")
  )

  # if SAM was succesful, extract results:
  if(is_empty(SAM$result$error)){
    # save duration of SAM:
    output_SAM <- SAM$result$result
    t_SAM <- difftime(Sys.time(), start_SAM, unit = "s")                        # save duration for SAM

    ## extract results:
    # regression parameters and innovation variances
    SAM_parameters <- parTable(output_SAM)
    SAM_phi11 <- SAM_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    SAM_phi12 <- SAM_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(est) |> as.numeric()
    SAM_phi21 <- SAM_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    SAM_phi22 <- SAM_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(est) |> as.numeric()

    SAM_zeta1 <- SAM_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(est) |> as.numeric()
    SAM_zeta2 <- SAM_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(est) |> as.numeric()
    SAM_zeta12 <- SAM_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(est) |> as.numeric()

    # standard errors
    SAM_phi11_se <- SAM_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    SAM_phi12_se <- SAM_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(se) |> as.numeric()
    SAM_phi21_se <- SAM_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    SAM_phi22_se <- SAM_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(se) |> as.numeric()

    SAM_zeta1_se <- SAM_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(se) |> as.numeric()
    SAM_zeta2_se <- SAM_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(se) |> as.numeric()
    SAM_zeta12_se <- SAM_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(se) |> as.numeric()
  } else{
    # if SAM was not succesful, set results to NA
    t_SAM <- NA

    SAM_phi11 <- NA
    SAM_phi12 <- NA
    SAM_phi21 <- NA
    SAM_phi22 <- NA

    SAM_zeta1 <- NA
    SAM_zeta2 <- NA
    SAM_zeta12 <- NA

    SAM_phi11_se <- NA
    SAM_phi12_se <- NA
    SAM_phi21_se <- NA
    SAM_phi22_se <- NA

    SAM_zeta1_se <- NA
    SAM_zeta2_se <- NA
    SAM_zeta12_se <- NA
  }

  #### 7) SEM ####
  start_SEM <- Sys.time()                                                       # start timer for SEM

  SEM <- run_SEM(model_SEM,
                 data = data_SEM,
                 missing = "ML",
                 cluster = "id")
  # extract error/warning messages (if applicable):
  SEM_warning <- ifelse(is_empty(SEM$warnings),
                        FALSE, TRUE)
  SEM_warning_text <- ifelse(is_empty(SEM$warnings),
                             "",
                             paste(c(SEM$warnings),
                                   collapse = "; ")
  )
  SEM_error <- ifelse(is_empty(SEM$result$error),
                      FALSE, TRUE)
  SEM_error_text <- ifelse(is_empty(SEM$result$error),
                           "",
                           paste(c(SEM$result$error),
                                 collapse = "; ")
  )

  # if SEM was succesful, extract results:
  if(is_empty(SEM$result$error)){
    t_SEM <- difftime(Sys.time(), start_SEM, unit = "s")                        # save duration of SEM
    output_SEM <- SEM$result$result

    ## extract results:
    # regression parameters and innovation variances
    SEM_parameters <- parTable(output_SEM)
    SEM_phi11 <- SEM_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    SEM_phi12 <- SEM_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(est) |> as.numeric()
    SEM_phi21 <- SEM_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(est) |> as.numeric()
    SEM_phi22 <- SEM_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(est) |> as.numeric()

    SEM_zeta1 <- SEM_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(est) |> as.numeric()
    SEM_zeta2 <- SEM_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(est) |> as.numeric()
    SEM_zeta12 <- SEM_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(est) |> as.numeric()

    # standard errors
    SEM_phi11_se <- SEM_parameters |>
      dplyr::filter(label == "phi_f1_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    SEM_phi12_se <- SEM_parameters |>
      dplyr::filter(label == "phi_f1_f2_lag") |>
      dplyr::select(se) |> as.numeric()
    SEM_phi21_se <- SEM_parameters |>
      dplyr::filter(label == "phi_f2_f1_lag") |>
      dplyr::select(se) |> as.numeric()
    SEM_phi22_se <- SEM_parameters |>
      dplyr::filter(label == "phi_f2_f2_lag") |>
      dplyr::select(se) |> as.numeric()

    SEM_zeta1_se <- SEM_parameters |>
      dplyr::filter(label == "zeta_f1_f1") |>
      dplyr::select(se) |> as.numeric()
    SEM_zeta2_se <- SEM_parameters |>
      dplyr::filter(label == "zeta_f2_f2") |>
      dplyr::select(se) |> as.numeric()
    SEM_zeta12_se <- SEM_parameters |>
      dplyr::filter(label == "zeta_f1_f2") |>
      dplyr::select(se) |> as.numeric()
  } else{
    # if SEM was not succesful, set results to NA
    t_SEM <- NA

    SEM_phi11 <- NA
    SEM_phi12 <- NA
    SEM_phi21 <- NA
    SEM_phi22 <- NA

    SEM_zeta1 <- NA
    SEM_zeta2 <- NA
    SEM_zeta12 <- NA

    SEM_phi11_se <- NA
    SEM_phi12_se <- NA
    SEM_phi21_se <- NA
    SEM_phi22_se <- NA

    SEM_zeta1_se <- NA
    SEM_zeta2_se <- NA
    SEM_zeta12_se <- NA
  }


  # compile output:
  output <- c("iteration" = iteration, "replication" = replication,
              "phi_size" = phi_size, "n" = n, "obs" = obs, "rho_gen" = rho_gen,
              "variance_means" = variance_means, "variance_zeta" = variance_zeta,
              "variance_phi" = variance_phi,
              "t_datagen" = t_datagen, "t_step1" = t_step1, "t_step2" = t_step2, 
              "t_step3" = t_step3, "t_SEcorr" = t_SEcorr,
              "t_NFS" = t_NFS, "t_SAM" = t_SAM, "t_SEM" = t_SEM,
              "phi11_pop" = phi11_pop, "phi12_pop" = phi12_pop, "phi21_pop" = phi21_pop, "phi22_pop" = phi22_pop,
              "zeta1_pop" = zeta1_pop, "zeta2_pop" = zeta2_pop, "zeta12_pop" = zeta12_pop,
              "LVAR_phi11" = LVAR_phi11, "LVAR_phi12" = LVAR_phi12,
              "LVAR_phi21" = LVAR_phi21, "LVAR_phi22" = LVAR_phi22,
              "LVAR_zeta1" = LVAR_zeta1, "LVAR_zeta2" = LVAR_zeta2, "LVAR_zeta12" = LVAR_zeta12,
              "LVAR_phi11_se" = LVAR_phi11_se, "LVAR_phi12_se" = LVAR_phi12_se,
              "LVAR_phi21_se" = LVAR_phi21_se, "LVAR_phi22_se" = LVAR_phi22_se,
              "LVAR_zeta1_se" = LVAR_zeta1_se, "LVAR_zeta2_se" = LVAR_zeta2_se, "LVAR_zeta12_se" = LVAR_zeta12_se,
              "LVAR_phi11_secorr" = LVAR_phi11_secorr, "LVAR_phi12_secorr" = LVAR_phi12_secorr,
              "LVAR_phi21_secorr" = LVAR_phi21_secorr, "LVAR_phi22_secorr" = LVAR_phi22_secorr,
              "LVAR_zeta1_secorr" = LVAR_zeta1_secorr, "LVAR_zeta2_secorr" = LVAR_zeta2_secorr, "LVAR_zeta12_secorr" = LVAR_zeta12_secorr,
              "NFS_phi11" = NFS_phi11, "NFS_phi12" = NFS_phi12,
              "NFS_phi21" = NFS_phi21, "NFS_phi22" = NFS_phi22,
              "NFS_zeta1" = NFS_zeta1, "NFS_zeta2" = NFS_zeta2, "NFS_zeta12" = NFS_zeta12,
              "NFS_phi11_se" = NFS_phi11_se, "NFS_phi12_se" = NFS_phi12_se,
              "NFS_phi21_se" = NFS_phi21_se, "NFS_phi22_se" = NFS_phi22_se,
              "NFS_zeta1_se" = NFS_zeta1_se, "NFS_zeta2_se" = NFS_zeta2_se, "NFS_zeta12_se" = NFS_zeta12_se,
              "SAM_phi11" = SAM_phi11, "SAM_phi12" = SAM_phi12,
              "SAM_phi21" = SAM_phi21, "SAM_phi22" = SAM_phi22,
              "SAM_zeta1" = SAM_zeta1, "SAM_zeta2" = SAM_zeta2, "SAM_zeta12" = SAM_zeta12,
              "SAM_phi11_se" = SAM_phi11_se, "SAM_phi12_se" = SAM_phi12_se,
              "SAM_phi21_se" = SAM_phi21_se, "SAM_phi22_se" = SAM_phi22_se,
              "SAM_zeta1_se" = SAM_zeta1_se, "SAM_zeta2_se" = SAM_zeta2_se, "SAM_zeta12_se" = SAM_zeta12_se,
              "SEM_phi11" = SEM_phi11, "SEM_phi12" = SEM_phi12,
              "SEM_phi21" = SEM_phi21, "SEM_phi22" = SEM_phi22,
              "SEM_zeta1" = SEM_zeta1, "SEM_zeta2" = SEM_zeta2, "SEM_zeta12" = SEM_zeta12,
              "SEM_phi11_se" = SEM_phi11_se, "SEM_phi12_se" = SEM_phi12_se,
              "SEM_phi21_se" = SEM_phi21_se, "SEM_phi22_se" = SEM_phi22_se,
              "SEM_zeta1_se" = SEM_zeta1_se, "SEM_zeta2_se" = SEM_zeta2_se, "SEM_zeta12_se" = SEM_zeta12_se,
              "LVAR_step1_warning" = LVAR_step1_warning, "LVAR_step2_warning" = LVAR_step2_warning,
              "LVAR_step3_warning" = LVAR_step3_warning, "SEcorr_warning" = SEcorr_warning, 
              "NFS_warning" = NFS_warning, "SAM_warning" = SAM_warning, "SEM_warning" = SEM_warning,
              "LVAR_step1_error" = LVAR_step1_error, "LVAR_step2_error" = LVAR_step2_error,
              "LVAR_step3_error" = LVAR_step3_error, "SEcorr_error" = SEcorr_error, 
              "NFS_error" = NFS_error, "SAM_error" = SAM_error, "SEM_error" = SEM_error,
              "seed" = seed_cond, "pos" = pos,
              "LVAR_step1_warning_text" = LVAR_step1_warning_text, 
              "LVAR_step2_warning_text" = LVAR_step2_warning_text,
              "LVAR_step3_warning_text" = LVAR_step3_warning_text,
              "SEcorr_warning_text" = SEcorr_warning_text, 
              "NFS_warning_text" = NFS_warning_text, 
              "SAM_warning_text" = SAM_warning_text, 
              "SEM_warning_text" = SEM_warning_text,
              "LVAR_step1_error_text" = LVAR_step1_error_text, 
              "LVAR_step2_error_text" = LVAR_step2_error_text, 
              "LVAR_step3_error_text" = LVAR_step3_error_text,
              "SEcorr_error_text" = SEcorr_error_text, 
              "NFS_error_text" = NFS_error_text, 
              "SAM_error_text" = SAM_error_text, 
              "SEM_error_text" = SEM_error_text)

  for(i in 104:117){
    output[i] <- str_squish(output[i])                                          # removes all whitespace and linebreaks from the error and warning strings
    output[i] <- gsub(",", "", output[i])                                       # removes all commata from error and warning strings (to prevent messing up the CSV file)
  }
  
  
  # check if file exists
  if(!file.exists(outputfile)){
    # if file does not yet exist
    write.table(t(output), file = outputfile, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    # lock the file to prevent multiple processes accessing it simultaneously
    lock <- flock::lock(outputfile)
    write.table(t(output), file = outputfile, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
    # unlock the file
    flock::unlock(lock)
  }
  
  if(verbose == TRUE){
    print(paste("Simulation", pos, "completed at", Sys.time()))                 # prints a message when a replication is done (as a sign that R did not crash)
  }
  
  return(output)
}