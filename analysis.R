library(tidyverse)

#### preparations ####
## load data and merge original and re-estimation data set
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
  arrange(iteration)                                                            # sort by iteration

# remove data sets where SAM had an error:
results_clean <- results[!results$SAM_error, ]

# check if any parameters were not estimated
results_clean |> 
  summarise(across(25:91, ~ sum(is.na(.)))) |> 
  print(width = Inf)

# how many replications remain in each condition?
results_clean |> 
  group_by(phi_size, n, obs, rho_gen, variance_means, variance_phi, variance_zeta) |> 
  summarise(n_rep = n()) |> 
  arrange(n_rep)
# largest number of removed replications for a condition is 5 (1%)

# split into two data frames (one with mean differences, one without)
results_nocent <- results_clean |> 
  dplyr::filter(variance_means == "no")
results_cent <- results_clean |> 
  dplyr::filter(variance_means == "yes")

#### apply correction for Nickell's bias ####
results_cent <- results_cent |> 
  rowwise() |> 
  mutate(LVAR_phi11corr = (LVAR_phi11*(obs-1) + 1)/(obs-2),
         LVAR_phi22corr = (LVAR_phi22*(obs-1) + 1)/(obs-2),
         NFS_phi11corr = (NFS_phi11*(obs-1) + 1)/(obs-2),
         NFS_phi22corr = (NFS_phi22*(obs-1) + 1)/(obs-2),
         SAM_phi11corr = (SAM_phi11*(obs-1) + 1)/(obs-2),
         SAM_phi22corr = (SAM_phi22*(obs-1) + 1)/(obs-2),
         SEM_phi11corr = (SEM_phi11*(obs-1) + 1)/(obs-2),
         SEM_phi22corr = (SEM_phi22*(obs-1) + 1)/(obs-2))

#### performance no centering ####
performance_nocent <- results_nocent |>
  group_by(phi_size, n, obs, variance_phi, variance_zeta, rho_gen) |> 
  ## compute the evaluation criteria:
  summarise(
    ## bias:
    AB_phi11_LVAR = mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_LVAR = (mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_LVAR = mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_LVAR = (mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_LVAR = mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_LVAR = (mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_LVAR = mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_LVAR = (mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_LVAR = mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_LVAR = (mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_LVAR = mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_LVAR = (mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_LVAR = mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_LVAR = (mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # LVAR
    AB_phi11_NFS = mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_NFS = (mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_NFS = mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_NFS = (mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_NFS = mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_NFS = (mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_NFS = mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_NFS = (mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_NFS = mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_NFS = (mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_NFS = mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_NFS = (mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_NFS = mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_NFS = (mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # NFS
    AB_phi11_SAM = mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SAM = (mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SAM = mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SAM = (mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SAM = mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SAM = (mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SAM = mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SAM = (mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_SAM = mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SAM = (mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SAM = mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SAM = (mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SAM = mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SAM = (mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SAM
    AB_phi11_SEM = mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SEM = (mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SEM = mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SEM = (mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SEM = mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SEM = (mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SEM = mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SEM = (mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_zeta1_SEM = mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SEM = (mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SEM = mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SEM = (mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SEM = mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SEM = (mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SEM
    ## RMSE:
    RMSE_phi11_LVAR = mean(sqrt((LVAR_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_LVAR = mean(sqrt((LVAR_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_LVAR = mean(sqrt((LVAR_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_LVAR = mean(sqrt((LVAR_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_LVAR = mean(sqrt((LVAR_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_LVAR = mean(sqrt((LVAR_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_LVAR = mean(sqrt((LVAR_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # LVAR
    RMSE_phi11_NFS = mean(sqrt((NFS_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_NFS = mean(sqrt((NFS_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_NFS = mean(sqrt((NFS_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_NFS = mean(sqrt((NFS_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_NFS = mean(sqrt((NFS_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_NFS = mean(sqrt((NFS_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_NFS = mean(sqrt((NFS_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # NFS
    RMSE_phi11_SAM = mean(sqrt((SAM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SAM = mean(sqrt((SAM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SAM = mean(sqrt((SAM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SAM = mean(sqrt((SAM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SAM = mean(sqrt((SAM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SAM = mean(sqrt((SAM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SAM = mean(sqrt((SAM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SAM
    RMSE_phi11_SEM = mean(sqrt((SEM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SEM = mean(sqrt((SEM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SEM = mean(sqrt((SEM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SEM = mean(sqrt((SEM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SEM = mean(sqrt((SEM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SEM = mean(sqrt((SEM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SEM = mean(sqrt((SEM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SEM
    ## Standard error recovery:
    SER_phi11_LVAR = mean(LVAR_phi11_se, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22_LVAR = mean(LVAR_phi22_se, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12_LVAR = mean(LVAR_phi12_se, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21_LVAR = mean(LVAR_phi21_se, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1_LVAR = mean(LVAR_zeta1_se, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2_LVAR = mean(LVAR_zeta2_se, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12_LVAR = mean(LVAR_zeta12_se, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    SER_phi11_LVARcorr = mean(LVAR_phi11_secorr, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22_LVARcorr = mean(LVAR_phi22_secorr, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12_LVARcorr = mean(LVAR_phi12_secorr, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21_LVARcorr = mean(LVAR_phi21_secorr, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1_LVARcorr = mean(LVAR_zeta1_secorr, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2_LVARcorr = mean(LVAR_zeta2_secorr, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12_LVARcorr = mean(LVAR_zeta12_secorr, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    # LVAR
    SER_phi11_NFS = mean(NFS_phi11_se, na.rm = TRUE) / sd(NFS_phi11, na.rm = TRUE),
    SER_phi22_NFS = mean(NFS_phi22_se, na.rm = TRUE) / sd(NFS_phi22, na.rm = TRUE),
    SER_phi12_NFS = mean(NFS_phi12_se, na.rm = TRUE) / sd(NFS_phi12, na.rm = TRUE),
    SER_phi21_NFS = mean(NFS_phi21_se, na.rm = TRUE) / sd(NFS_phi21, na.rm = TRUE),
    SER_zeta1_NFS = mean(NFS_zeta1_se, na.rm = TRUE) / sd(NFS_zeta1, na.rm = TRUE),
    SER_zeta2_NFS = mean(NFS_zeta2_se, na.rm = TRUE) / sd(NFS_zeta2, na.rm = TRUE),
    SER_zeta12_NFS = mean(NFS_zeta12_se, na.rm = TRUE) / sd(NFS_zeta12, na.rm = TRUE),
    # NFS
    SER_phi11_SAM = mean(SAM_phi11_se, na.rm = TRUE) / sd(SAM_phi11, na.rm = TRUE),
    SER_phi22_SAM = mean(SAM_phi22_se, na.rm = TRUE) / sd(SAM_phi22, na.rm = TRUE),
    SER_phi12_SAM = mean(SAM_phi12_se, na.rm = TRUE) / sd(SAM_phi12, na.rm = TRUE),
    SER_phi21_SAM = mean(SAM_phi21_se, na.rm = TRUE) / sd(SAM_phi21, na.rm = TRUE),
    SER_zeta1_SAM = mean(SAM_zeta1_se, na.rm = TRUE) / sd(SAM_zeta1, na.rm = TRUE),
    SER_zeta2_SAM = mean(SAM_zeta2_se, na.rm = TRUE) / sd(SAM_zeta2, na.rm = TRUE),
    SER_zeta12_SAM = mean(SAM_zeta12_se, na.rm = TRUE) / sd(SAM_zeta12, na.rm = TRUE),
    # SAM
    SER_phi11_SEM = mean(SEM_phi11_se, na.rm = TRUE) / sd(SEM_phi11, na.rm = TRUE),
    SER_phi22_SEM = mean(SEM_phi22_se, na.rm = TRUE) / sd(SEM_phi22, na.rm = TRUE),
    SER_phi12_SEM = mean(SEM_phi12_se, na.rm = TRUE) / sd(SEM_phi12, na.rm = TRUE),
    SER_phi21_SEM = mean(SEM_phi21_se, na.rm = TRUE) / sd(SEM_phi21, na.rm = TRUE),
    SER_zeta1_SEM = mean(SEM_zeta1_se, na.rm = TRUE) / sd(SEM_zeta1, na.rm = TRUE),
    SER_zeta2_SEM = mean(SEM_zeta2_se, na.rm = TRUE) / sd(SEM_zeta2, na.rm = TRUE),
    SER_zeta12_SEM = mean(SEM_zeta12_se, na.rm = TRUE) / sd(SEM_zeta12, na.rm = TRUE),
    # # SEM
    .groups = "drop") |> 
  # bring into wide format and then back into long format (bring the methods into distinct rows)
  pivot_longer(cols = 7:125,
               names_to = c("type", "variable", "method"),
               names_pattern = "^(.*)_(phi\\d{2}|zeta\\d{1,2})_(.*)$") |> 
  mutate(name = paste(type, variable, sep = "_")) |>
  select(-type, -variable) |> 
  pivot_wider(names_from = name,
              values_from = value)


## performance overall
performance_nocent_overall <- performance_nocent |> 
  group_by(method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE))) |> 
  mutate(manipulated_aspect = "overall",
         level = as.character(NA),
         .before = 1)
## performance by effect size
performance_nocent_phi_size <- performance_nocent |> 
  group_by(phi_size, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "phi_size", .before = 1) |> 
  rename(level = phi_size) |> 
  mutate(level = as.character(level))

## performance by n
performance_nocent_n <- performance_nocent |> 
  group_by(n, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "n", .before = 1) |> 
  rename(level = n) |> 
  mutate(level = as.character(level))

## performance by obs
performance_nocent_obs <- performance_nocent |> 
  group_by(obs, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "obs", .before = 1) |> 
  rename(level = obs) |> 
  mutate(level = as.character(level))

## performance by rho
performance_nocent_rho_gen <- performance_nocent |> 
  group_by(rho_gen, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "rho_gen", .before = 1) |> 
  rename(level = rho_gen) |> 
  mutate(level = as.character(level))

## performance by variance_phi
performance_nocent_variance_phi <- performance_nocent |> 
  group_by(variance_phi, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "variance_phi", .before = 1) |> 
  rename(level = variance_phi) |> 
  mutate(level = as.character(level))

## performance by variance_zeta
performance_nocent_variance_zeta <- performance_nocent |> 
  group_by(variance_zeta, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "variance_zeta", .before = 1) |> 
  rename(level = variance_zeta) |> 
  mutate(level = as.character(level))



## combine into three different tables (for bias, RMSE, SE recovery)
performance_nocent_bias <- bind_rows(performance_nocent_overall,
                                     performance_nocent_phi_size,
                                     performance_nocent_n,
                                     performance_nocent_obs,
                                     performance_nocent_rho_gen,
                                     performance_nocent_variance_phi,
                                     performance_nocent_variance_zeta) |> 
  select(1:17) |>  # only select bias columns
  dplyr::filter(method != "LVARcorr") # remove these rows because they're only relevant for SE correction

performance_nocent_RMSE <- bind_rows(performance_nocent_overall,
                                     performance_nocent_phi_size,
                                     performance_nocent_n,
                                     performance_nocent_obs,
                                     performance_nocent_rho_gen,
                                     performance_nocent_variance_phi,
                                     performance_nocent_variance_zeta) |> 
  select(1:3, 18:24) |>  # only select RMSE columns
  dplyr::filter(method != "LVARcorr") # remove these rows because they're only relevant for SE correction

performance_nocent_SER <- bind_rows(performance_nocent_overall,
                                    performance_nocent_phi_size,
                                    performance_nocent_n,
                                    performance_nocent_obs,
                                    performance_nocent_rho_gen,
                                    performance_nocent_variance_phi,
                                    performance_nocent_variance_zeta) |> 
  select(1:3, 25:31)  # only select SE recovery columns

# write_csv(performance_nocent_bias, file = "Data/table_nocent_bias.csv")
# write_csv(performance_nocent_RMSE, file = "Data/table_nocent_RMSE.csv")
# write_csv(performance_nocent_SER, file = "Data/table_nocent_SER.csv")


#### performance with centering ####
performance_cent <- results_cent |>
  group_by(phi_size, n, obs, variance_phi, variance_zeta, rho_gen) |> 
  ## compute the evaluation criteria:
  summarise(
    ## bias:
    AB_phi11_LVAR = mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_LVAR = (mean(LVAR_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_LVAR = mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_LVAR = (mean(LVAR_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_LVAR = mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_LVAR = (mean(LVAR_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_LVAR = mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_LVAR = (mean(LVAR_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_phi11_LVARcorr = mean(LVAR_phi11corr, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_LVARcorr = (mean(LVAR_phi11corr, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_LVARcorr = mean(LVAR_phi22corr, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_LVARcorr = (mean(LVAR_phi22corr, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_zeta1_LVAR = mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_LVAR = (mean(LVAR_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_LVAR = mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_LVAR = (mean(LVAR_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_LVAR = mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_LVAR = (mean(LVAR_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # LVAR
    AB_phi11_NFS = mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_NFS = (mean(NFS_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_NFS = mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_NFS = (mean(NFS_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_NFS = mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_NFS = (mean(NFS_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_NFS = mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_NFS = (mean(NFS_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_phi11_NFScorr = mean(NFS_phi11corr, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_NFScorr = (mean(NFS_phi11corr, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_NFScorr = mean(NFS_phi22corr, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_NFScorr = (mean(NFS_phi22corr, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_zeta1_NFS = mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_NFS = (mean(NFS_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_NFS = mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_NFS = (mean(NFS_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_NFS = mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_NFS = (mean(NFS_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # NFS
    AB_phi11_SAM = mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SAM = (mean(SAM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SAM = mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SAM = (mean(SAM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SAM = mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SAM = (mean(SAM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SAM = mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SAM = (mean(SAM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_phi11_SAMcorr = mean(SAM_phi11corr, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SAMcorr = (mean(SAM_phi11corr, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SAMcorr = mean(SAM_phi22corr, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SAMcorr = (mean(SAM_phi22corr, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_zeta1_SAM = mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SAM = (mean(SAM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SAM = mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SAM = (mean(SAM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SAM = mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SAM = (mean(SAM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SAM
    AB_phi11_SEM = mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SEM = (mean(SEM_phi11, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SEM = mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SEM = (mean(SEM_phi22, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_phi12_SEM = mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop),
    RB_phi12_SEM = (mean(SEM_phi12, na.rm = TRUE) - unique(phi12_pop)) / unique(phi12_pop),
    AB_phi21_SEM = mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop),
    RB_phi21_SEM = (mean(SEM_phi21, na.rm = TRUE) - unique(phi21_pop)) / unique(phi21_pop),
    AB_phi11_SEMcorr = mean(SEM_phi11corr, na.rm = TRUE) - unique(phi11_pop),
    RB_phi11_SEMcorr = (mean(SEM_phi11corr, na.rm = TRUE) - unique(phi11_pop)) / unique(phi11_pop),
    AB_phi22_SEMcorr = mean(SEM_phi22corr, na.rm = TRUE) - unique(phi22_pop),
    RB_phi22_SEMcorr = (mean(SEM_phi22corr, na.rm = TRUE) - unique(phi22_pop)) / unique(phi22_pop),
    AB_zeta1_SEM = mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop),
    RB_zeta1_SEM = (mean(SEM_zeta1, na.rm = TRUE) - unique(zeta1_pop)) / unique(zeta1_pop),
    AB_zeta2_SEM = mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop),
    RB_zeta2_SEM = (mean(SEM_zeta2, na.rm = TRUE) - unique(zeta2_pop)) / unique(zeta2_pop),
    AB_zeta12_SEM = mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop),
    RB_zeta12_SEM = (mean(SEM_zeta12, na.rm = TRUE) - unique(zeta12_pop)) / unique(zeta12_pop),
    # SEM
    ## RMSE:
    RMSE_phi11_LVAR = mean(sqrt((LVAR_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_LVAR = mean(sqrt((LVAR_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_LVAR = mean(sqrt((LVAR_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_LVAR = mean(sqrt((LVAR_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_phi11_LVARcorr = mean(sqrt((LVAR_phi11corr - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_LVARcorr = mean(sqrt((LVAR_phi22corr - phi22_pop)^2), na.rm = TRUE),
    RMSE_zeta1_LVAR = mean(sqrt((LVAR_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_LVAR = mean(sqrt((LVAR_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_LVAR = mean(sqrt((LVAR_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # LVAR
    RMSE_phi11_NFS = mean(sqrt((NFS_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_NFS = mean(sqrt((NFS_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_NFS = mean(sqrt((NFS_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_NFS = mean(sqrt((NFS_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_phi11_NFScorr = mean(sqrt((NFS_phi11corr - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_NFScorr = mean(sqrt((NFS_phi22corr - phi22_pop)^2), na.rm = TRUE),
    RMSE_zeta1_NFS = mean(sqrt((NFS_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_NFS = mean(sqrt((NFS_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_NFS = mean(sqrt((NFS_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # NFS
    RMSE_phi11_SAM = mean(sqrt((SAM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SAM = mean(sqrt((SAM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SAM = mean(sqrt((SAM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SAM = mean(sqrt((SAM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_phi11_SAMcorr = mean(sqrt((SAM_phi11corr - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SAMcorr = mean(sqrt((SAM_phi22corr - phi22_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SAM = mean(sqrt((SAM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SAM = mean(sqrt((SAM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SAM = mean(sqrt((SAM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SAM
    RMSE_phi11_SEM = mean(sqrt((SEM_phi11 - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SEM = mean(sqrt((SEM_phi22 - phi22_pop)^2), na.rm = TRUE),
    RMSE_phi12_SEM = mean(sqrt((SEM_phi12 - phi12_pop)^2), na.rm = TRUE),
    RMSE_phi21_SEM = mean(sqrt((SEM_phi21 - phi21_pop)^2), na.rm = TRUE),
    RMSE_phi11_SEMcorr = mean(sqrt((SEM_phi11corr - phi11_pop)^2), na.rm = TRUE),
    RMSE_phi22_SEMcorr = mean(sqrt((SEM_phi22corr - phi22_pop)^2), na.rm = TRUE),
    RMSE_zeta1_SEM = mean(sqrt((SEM_zeta1 - zeta1_pop)^2), na.rm = TRUE),
    RMSE_zeta2_SEM = mean(sqrt((SEM_zeta2 - zeta2_pop)^2), na.rm = TRUE),
    RMSE_zeta12_SEM = mean(sqrt((SEM_zeta12 - zeta12_pop)^2), na.rm = TRUE),
    # SEM
    ## Standard error recovery:
    SER_phi11_LVAR = mean(LVAR_phi11_se, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22_LVAR = mean(LVAR_phi22_se, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12_LVAR = mean(LVAR_phi12_se, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21_LVAR = mean(LVAR_phi21_se, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1_LVAR = mean(LVAR_zeta1_se, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2_LVAR = mean(LVAR_zeta2_se, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12_LVAR = mean(LVAR_zeta12_se, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    SER_phi11_LVARcorr = mean(LVAR_phi11_secorr, na.rm = TRUE) / sd(LVAR_phi11, na.rm = TRUE),
    SER_phi22_LVARcorr = mean(LVAR_phi22_secorr, na.rm = TRUE) / sd(LVAR_phi22, na.rm = TRUE),
    SER_phi12_LVARcorr = mean(LVAR_phi12_secorr, na.rm = TRUE) / sd(LVAR_phi12, na.rm = TRUE),
    SER_phi21_LVARcorr = mean(LVAR_phi21_secorr, na.rm = TRUE) / sd(LVAR_phi21, na.rm = TRUE),
    SER_zeta1_LVARcorr = mean(LVAR_zeta1_secorr, na.rm = TRUE) / sd(LVAR_zeta1, na.rm = TRUE),
    SER_zeta2_LVARcorr = mean(LVAR_zeta2_secorr, na.rm = TRUE) / sd(LVAR_zeta2, na.rm = TRUE),
    SER_zeta12_LVARcorr = mean(LVAR_zeta12_secorr, na.rm = TRUE) / sd(LVAR_zeta12, na.rm = TRUE),
    # LVAR
    SER_phi11_NFS = mean(NFS_phi11_se, na.rm = TRUE) / sd(NFS_phi11, na.rm = TRUE),
    SER_phi22_NFS = mean(NFS_phi22_se, na.rm = TRUE) / sd(NFS_phi22, na.rm = TRUE),
    SER_phi12_NFS = mean(NFS_phi12_se, na.rm = TRUE) / sd(NFS_phi12, na.rm = TRUE),
    SER_phi21_NFS = mean(NFS_phi21_se, na.rm = TRUE) / sd(NFS_phi21, na.rm = TRUE),
    SER_zeta1_NFS = mean(NFS_zeta1_se, na.rm = TRUE) / sd(NFS_zeta1, na.rm = TRUE),
    SER_zeta2_NFS = mean(NFS_zeta2_se, na.rm = TRUE) / sd(NFS_zeta2, na.rm = TRUE),
    SER_zeta12_NFS = mean(NFS_zeta12_se, na.rm = TRUE) / sd(NFS_zeta12, na.rm = TRUE),
    # NFS
    SER_phi11_SAM = mean(SAM_phi11_se, na.rm = TRUE) / sd(SAM_phi11, na.rm = TRUE),
    SER_phi22_SAM = mean(SAM_phi22_se, na.rm = TRUE) / sd(SAM_phi22, na.rm = TRUE),
    SER_phi12_SAM = mean(SAM_phi12_se, na.rm = TRUE) / sd(SAM_phi12, na.rm = TRUE),
    SER_phi21_SAM = mean(SAM_phi21_se, na.rm = TRUE) / sd(SAM_phi21, na.rm = TRUE),
    SER_zeta1_SAM = mean(SAM_zeta1_se, na.rm = TRUE) / sd(SAM_zeta1, na.rm = TRUE),
    SER_zeta2_SAM = mean(SAM_zeta2_se, na.rm = TRUE) / sd(SAM_zeta2, na.rm = TRUE),
    SER_zeta12_SAM = mean(SAM_zeta12_se, na.rm = TRUE) / sd(SAM_zeta12, na.rm = TRUE),
    # SAM
    SER_phi11_SEM = mean(SEM_phi11_se, na.rm = TRUE) / sd(SEM_phi11, na.rm = TRUE),
    SER_phi22_SEM = mean(SEM_phi22_se, na.rm = TRUE) / sd(SEM_phi22, na.rm = TRUE),
    SER_phi12_SEM = mean(SEM_phi12_se, na.rm = TRUE) / sd(SEM_phi12, na.rm = TRUE),
    SER_phi21_SEM = mean(SEM_phi21_se, na.rm = TRUE) / sd(SEM_phi21, na.rm = TRUE),
    SER_zeta1_SEM = mean(SEM_zeta1_se, na.rm = TRUE) / sd(SEM_zeta1, na.rm = TRUE),
    SER_zeta2_SEM = mean(SEM_zeta2_se, na.rm = TRUE) / sd(SEM_zeta2, na.rm = TRUE),
    SER_zeta12_SEM = mean(SEM_zeta12_se, na.rm = TRUE) / sd(SEM_zeta12, na.rm = TRUE),
    # # SEM
    .groups = "drop") |> 
  # bring into wide format and then back into long format (bring the methods into distinct rows)
  pivot_longer(cols = 7:149,
               names_to = c("type", "variable", "method"),
               names_pattern = "^(.*)_(phi\\d{2}|zeta\\d{1,2})_(.*)$") |> 
  mutate(name = paste(type, variable, sep = "_")) |>
  select(-type, -variable) |> 
  pivot_wider(names_from = name,
              values_from = value)


## performance overall
performance_cent_overall <- performance_cent |> 
  group_by(method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE))) |> 
  mutate(manipulated_aspect = "overall",
         level = as.character(NA),
         .before = 1)
## performance by effect size
performance_cent_phi_size <- performance_cent |> 
  group_by(phi_size, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "phi_size", .before = 1) |> 
  rename(level = phi_size) |> 
  mutate(level = as.character(level))

## performance by n
performance_cent_n <- performance_cent |> 
  group_by(n, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "n", .before = 1) |> 
  rename(level = n) |> 
  mutate(level = as.character(level))

## performance by obs
performance_cent_obs <- performance_cent |> 
  group_by(obs, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "obs", .before = 1) |> 
  rename(level = obs) |> 
  mutate(level = as.character(level))

## performance by rho
performance_cent_rho_gen <- performance_cent |> 
  group_by(rho_gen, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "rho_gen", .before = 1) |> 
  rename(level = rho_gen) |> 
  mutate(level = as.character(level))

## performance by variance_phi
performance_cent_variance_phi <- performance_cent |> 
  group_by(variance_phi, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "variance_phi", .before = 1) |> 
  rename(level = variance_phi) |> 
  mutate(level = as.character(level))

## performance by variance_zeta
performance_cent_variance_zeta <- performance_cent |> 
  group_by(variance_zeta, method) |> 
  summarise(across(AB_phi11:SER_zeta12, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "variance_zeta", .before = 1) |> 
  rename(level = variance_zeta) |> 
  mutate(level = as.character(level))



## combine into three different tables (for bias, RMSE, SE recovery)
performance_cent_bias <- bind_rows(performance_cent_overall,
                                   performance_cent_phi_size,
                                   performance_cent_n,
                                   performance_cent_obs,
                                   performance_cent_rho_gen,
                                   performance_cent_variance_phi,
                                   performance_cent_variance_zeta) |> 
  select(1:17) # only select bias columns

performance_cent_RMSE <- bind_rows(performance_cent_overall,
                                   performance_cent_phi_size,
                                   performance_cent_n,
                                   performance_cent_obs,
                                   performance_cent_rho_gen,
                                   performance_cent_variance_phi,
                                   performance_cent_variance_zeta) |> 
  select(1:3, 18:24)  # only select RMSE columns
  
performance_cent_SER <- bind_rows(performance_cent_overall,
                                  performance_cent_phi_size,
                                  performance_cent_n,
                                  performance_cent_obs,
                                  performance_cent_rho_gen,
                                  performance_cent_variance_phi,
                                  performance_cent_variance_zeta) |> 
  select(1:3, 25:31) |>  # only select SE recovery columns
  filter(method %in% c("LVAR", "LVARcorr", "NFS", "SAM", "SEM"))

# write_csv(performance_cent_bias, file = "Data/table_cent_bias.csv")
# write_csv(performance_cent_RMSE, file = "Data/table_cent_RMSE.csv")
# write_csv(performance_cent_SER, file = "Data/table_cent_SER.csv")