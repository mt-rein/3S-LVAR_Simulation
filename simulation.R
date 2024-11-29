library(conflicted)
library(dplyr)
library(flock)
library(lavaan)
library(MASS)
library(numDeriv)
library(purrr)
library(RcppAlgos)
library(tidyr)
library(truncnorm)
library(stringr)

library(parallel)
library(parabar)


source("sim_VAR.R")
source("center_within.R")
source("step1.R")
source("step2.R")
source("step3.R")
source("SEcorrection.R")
source("auxilliary functions.R")
source("do_sim.R")

# condition grid:
cond <- expand.grid(replication = 1:500,
                    phi_size = c("small", "large"),
                    n = c(25, 50),
                    obs = c(25, 50),
                    rho_gen = c("small", "medium", "large", "very large"),
                    variance_means = c("no", "yes"),
                    variance_zeta = c("no", "yes"),
                    variance_phi = c("no", "yes"))

# add seeds:
set.seed(123)
cond$seed <- sample(1:nrow(cond)*5, size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration


#### set up parallel computing ####
## open cluster
numCores <- detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


## load libraries and additional functions into cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(numDeriv)
  library(purrr)
  library(RcppAlgos)
  library(tidyr)
  library(truncnorm)
  library(stringr)
  
  source("sim_VAR.R")
  source("center_within.R")
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("SEcorrection.R")
  source("auxilliary functions.R")
  source("do_sim.R")
})

## load condition grid into cluster
export(backend, "cond")

# simulation with parallel processing:
outputfile = "output_sim.csv"
if(file.exists(outputfile)){
  print("Careful! File already exists!")
} else {
  output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = outputfile, verbose = FALSE)
}

## close cluster
stop_backend(backend)


#### Re-estimation ####
library(readr)
results <- read_csv("Data/output_sim.csv",
                    col_types = cols(LVAR_step1_warning_text = col_character(),
                                     LVAR_step2_warning_text = col_character(),
                                     LVAR_step3_warning_text = col_character(),
                                     SEcorr_warning_text = col_character(),
                                     NFS_warning_text = col_character(),
                                     SAM_warning_text = col_character(),
                                     SEM_warning_text = col_character(),
                                     LVAR_step1_error_text = col_character(),
                                     LVAR_step2_error_text = col_character(),
                                     LVAR_step3_error_text = col_character(),
                                     SEcorr_error_text = col_character(),
                                     NFS_error_text = col_character(),
                                     SAM_error_text = col_character(),
                                     SEM_error_text = col_character())) |> 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large")),
         variance_means = factor(variance_means, levels = c("yes", "no")),
         variance_zeta = factor(variance_zeta, levels = c("yes", "no")),
         variance_phi = factor(variance_phi, levels = c("yes", "no"))
  ) |>  
  arrange(iteration)

results |> 
  dplyr::select(ends_with("error")) |> 
  colSums()
# 58 errors for SAM, 0 for rest

failed_iterations <- results$iteration[results$SAM_error]

cond_reestimation <- cond[cond$iteration %in% failed_iterations, ]

# replace do_sim() function with modified function (new seed before SAM estimation)
source("do_sim_reestimation.R")

## open cluster
numCores <- detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


## load libraries and additional functions into cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(numDeriv)
  library(purrr)
  library(RcppAlgos)
  library(tidyr)
  library(truncnorm)
  library(stringr)
  
  source("sim_VAR.R")
  source("center_within.R")
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("SEcorrection.R")
  source("auxilliary functions.R")
  source("do_sim_reestimation.R")
})

## load condition grid into cluster
export(backend, "cond_reestimation")

# simulation with parallel processing:
outputfile = "output_sim_reestimation.csv"
if(file.exists(outputfile)){
  print("Careful! File already exists!")
} else {
  output <- par_lapply(backend, 1:nrow(cond_reestimation), do_sim, cond = cond_reestimation, outputfile = outputfile, verbose = FALSE)
}

## close cluster
stop_backend(backend)


## check if re-estimation solved the errors
results_reestimation <- read_csv("Data/output_sim_reestimation.csv",
                                 col_types = cols(LVAR_step1_warning_text = col_character(),
                                                  LVAR_step2_warning_text = col_character(),
                                                  LVAR_step3_warning_text = col_character(),
                                                  SEcorr_warning_text = col_character(),
                                                  NFS_warning_text = col_character(),
                                                  SAM_warning_text = col_character(),
                                                  SEM_warning_text = col_character(),
                                                  LVAR_step1_error_text = col_character(),
                                                  LVAR_step2_error_text = col_character(),
                                                  LVAR_step3_error_text = col_character(),
                                                  SEcorr_error_text = col_character(),
                                                  NFS_error_text = col_character(),
                                                  SAM_error_text = col_character(),
                                                  SEM_error_text = col_character())) |> 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large")),
         variance_means = factor(variance_means, levels = c("yes", "no")),
         variance_zeta = factor(variance_zeta, levels = c("yes", "no")),
         variance_phi = factor(variance_phi, levels = c("yes", "no"))
  ) |> 
  arrange(iteration)                                                            # sort by iteration

# check for any errors in re-estimated iterations:
results_reestimation |> 
  dplyr::select(ends_with("error")) |> 
  colSums()
# Still errors in all 58 data sets