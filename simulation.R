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
