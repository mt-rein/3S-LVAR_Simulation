library(conflicted)
library(dplyr)
library(lavaan)
library(MASS)
library(numDeriv)
library(purrr)
library(RcppAlgos)
library(truncnorm)
library(stringr)

library(parallel)
library(parabar)


source("sim_VAR.R")
source("center_within.R")
source("step1.R")
source("step2.R")
source("step3.R")
source("auxilliary functions.R")
source("do_sim.R")

# condition grid:
cond <- expand.grid(replication = 1:250,
                    phi_size = c("small", "large"),
                    n = c(25, 50),
                    obs = c(25, 50, 100),
                    rho_gen = c("small", "medium", "large", "very large"),
                    meanvar = c("small", "large"))

# add seeds:
set.seed(123)
cond$seed <- sample(1:nrow(cond)*5, size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration


#### set up parallel computing ####
## open cluster
numCores <- 4
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


## load libraries and additional functions into cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(lavaan)
  library(MASS)
  library(numDeriv)
  library(purrr)
  library(RcppAlgos)
  library(truncnorm)
  library(stringr)
  
  source("sim_VAR.R")
  source("center_within.R")
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("auxilliary functions.R")
  source("do_sim.R")
})

## load condition grid into cluster
export(backend, "cond")

# simulation with parallel processing:
outputfile = "test.csv"
if(file.exists(outputfile)){
  print("Careful! File already exists!")
} else {
  output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = outputfile, verbose = FALSE)
}


## close cluster
stop_backend(backend)


results <- readr::read_csv("test.csv") |> 
  arrange(iteration) |> 
  mutate(phi_size = factor(phi_size, levels = c("small", "large")),
         rho_gen = factor(rho_gen, levels = c("small", "medium", "large", "very large")))

sum(results$LVAR_step3_warning)

blubb <- results |> 
  group_by(phi_size, n, obs, rho_gen) |> 
  summarise(
    AB_phi11 = mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11 = (mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi12 = mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12 = (mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21 = mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21 = (mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_phi22 = mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22 = (mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_zeta1 = mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1 = (mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2 = mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2 = (mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12 = mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12 = (mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    nickellsbias_phi11 = mean(nickellsbias_phi11, na.rm = TRUE),
    nickellsbias_phi22 = mean(nickellsbias_phi22, na.rm = TRUE),
    rhof1 = mean(rho_f1),
    rhof2 = mean(rho_f2),
    .groups = "drop")


blubb_overall <- blubb |> 
  summarise(across(AB_phi11:rhof2, ~ mean(.x, na.rm = TRUE)))

#### bias by phi_size ####
blubb_phi_size <- blubb |> 
  group_by(phi_size) |> 
  summarise(across(AB_phi11:rhof2, ~ mean(.x, na.rm = TRUE)))

#### bias by obs ####
blubb_obs <- blubb |> 
  group_by(obs) |> 
  summarise(across(AB_phi11:rhof2, ~ mean(.x, na.rm = TRUE)))

#### bias by n ####
blubb_n <- blubb |> 
  group_by(n) |> 
  summarise(across(AB_phi11:rhof2, ~ mean(.x, na.rm = TRUE)))

#### bias by rho_gen ####
blubb_rho_gen <- blubb |> 
  group_by(rho_gen) |> 
  summarise(across(AB_phi11:rhof2, ~ mean(.x, na.rm = TRUE)))















## perform simulation:
output <- lapply(1:nrow(cond), do_sim, cond = cond, outputfile = "output_sim.csv", verbose = TRUE)
# note: the output object is irrelevant, the results are written into the CSV file

# we checked the simulation for errors and warnings, see file "check_warnings.R".
# 516 iterations failed to convergence when using SAM. We re-estimated these 
# using the default bounds when estimating the MM

#### re-estimation with SAM ####
## find iterations where SAM failed to converge
library(tidyverse)
results <- read_csv("Data/output_sim.csv")
failed_iterations <- results$iteration[results$SAM_error]
cond_reestimation <- cond[cond$iteration %in% failed_iterations, ]

## load reestimation function
source("do_sim_reestimation.R")

## run re-estimation
output <- lapply(1:nrow(cond_reestimation), do_sim_reestimation, cond = cond_reestimation, outputfile = "output_sim.csv", verbose = TRUE)