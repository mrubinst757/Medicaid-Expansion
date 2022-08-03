# program: 03-jackknife.r
# purpose: leave-one-out-state jackknife estimates
# author: max rubinstein
# date modified: december 14, 2020

source("03_Analysis/03_jackknife.R")

variables <- read_csv("../02_Specs/tol_specs.csv") %>%
  mutate(`Reference Variable` = case_when(grepl("child", Variable) ~ 0, TRUE ~ `Reference Variable`)) %>%
  filter(`Reference Variable` == 0) %>%
  .$Variable %>%
  sort()

tol_list <- read_csv("../02_Specs/tol_specs.csv") %>%
  filter(Variable %in% variables) %>%
  arrange(Variable) %>%
  .$`Base Tol` %>%
  set_names(variables)

tol_list <- list(tol_list)

cov_models <- c("sigma_uu_i_modeled", "sigma_uu_avg", "sigma_zero")

cov_groups <- c("None")

model_names <- c("SBW", "H-SBW", "BC-SBW", "BC-HSBW")

variable_list <- map(c("", 1:5), ~read_varnames(.x))
names(variable_list) <- cov_groups

# run program and save output -------------------------------------------------------------------------
c1_data <-readRDS("../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-true.rds")[[1]] %>%
  filter(treatment < 2) %>%
  arrange(state, cpuma)

c2_data <- c1_data %>%
  filter(!state %in% c("CA", "CT", "MN", "NJ", "WA"))

c1_tdat <- subset(c1_data, treatment == 1)
c2_tdat <- c1_tdat %>%
  filter(!state %in% c("CA", "CT", "MN", "NJ", "WA"))
c1_cdat <- subset(c1_data, treatment == 0)
c2_cdat <- c1_cdat

c2_indices <- c1_tdat %>%
  mutate(c2_remove = ifelse(state %in% c("CA", "CT", "MN", "NJ", "WA"), FALSE, TRUE)) %>%
  .$c2_remove

sigma_uu_tx_c1 <- readRDS("../01_ProcessedData/sigma-uu-i-txonly-true.rds")
sigma_uu_tx_c2 <- sigma_uu_tx_c1[c2_indices]
sigma_uu_ct_c1 <- readRDS("../01_ProcessedData/sigma-uu-i-ctonly-true.rds")
sigma_uu_ct_c2 <- sigma_uu_ct_c1

calibrated_data_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.RDS")
cdat_imputed_c1 <- filter(calibrated_data_c1, treatment == 0, set == "true")
tdat_imputed_c1 <- filter(calibrated_data_c1, treatment == 1, set == "true")

calibrated_data_c2 <- readRDS("../01_ProcessedData/calibrated-data-c2.RDS")
cdat_imputed_c2 <- filter(calibrated_data_c2, treatment == 0)
tdat_imputed_c2 <- filter(calibrated_data_c2, treatment == 1)

count_xwalk <- read_csv("../02_Specs/tol_specs.csv") %>%
  select(Variable, sample_var = `Sample Size`) %>%
  filter(Variable %in% variables)

targets <- colMeans(c1_cdat[variables])

# create c1 jackknife data
etc_jackknife_dat_c1 <- cluster_jackknife_etc(c1_tdat, sigma_uu_tx_c1, variables, targets, count_xwalk) 
etc_jackknife_dat_c1 <- etc_jackknife_dat_c1[-grep("NH", names(etc_jackknife_dat_c1))]
etc_jackknife_dat_c1 <- map(etc_jackknife_dat_c1, ~mutate(.x, data = map(data, ~filter(.x, state != "NH"))))

# create c2 jackknife data
etc_jackknife_dat_c2 <- cluster_jackknife_etc(c2_tdat, sigma_uu_tx_c2,
                                              variables, targets, count_xwalk) 
etc_jackknife_dat_c2 <- etc_jackknife_dat_c2[-grep("NH", names(etc_jackknife_dat_c2))]
etc_jackknife_dat_c2 <- map(etc_jackknife_dat_c2, ~mutate(.x, data = map(data, ~filter(.x, state != "NH"))))

write_rds(etc_jackknife_dat_c1, "../01_ProcessedData/etc-jackknife-data-c1.RDS")
write_rds(etc_jackknife_dat_c2, "../01_ProcessedData/etc-jackknife-data-c2.RDS")

etc_jackknife_dat_c1 <- readRDS("../01_ProcessedData/etc-jackknife-data-c1.RDS")
etc_jackknife_dat_c2 <- readRDS("../01_ProcessedData/etc-jackknife-data-c2.RDS")

# calculate weights on jackknife data
etc_results_c1 <- etc_weights(etc_jackknife_dat_c1, tol_list, cov_groups, targets)
etc_results_c2 <- etc_weights(etc_jackknife_dat_c2, tol_list, cov_groups, targets)

saveRDS(etc_results_c1, "../04_Output/etc-jackknife-c1.rds")
saveRDS(etc_results_c2, "../04_Output/etc-jackknife-c2.rds")
