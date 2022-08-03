# program: 01-adjust-data.R
# purpose: run data calibration programs and output files
# author: max rubinstein
# date modified: december 14, 2020

# load libraries and data ------------------------------------------------------------------------------------
#remotes::install_github("mrubinst757/optweight") # install modified version of optweight
source("03_Analysis/01_calibrate-data.R")
library(tidyverse)
library(assertthat)

modify_variables <- function(filename, variables) {
  suffix <- gsub("\\.rds", "", stringr::str_extract(filename, "[a-z]{4,5}\\.rds"))
  
  count_xwalk <- read_csv('../02_Specs/tol_specs.csv') %>%
    select(Variable, sample_var = `Sample Size`) %>%
    filter(Variable %in% variables)
  
  adjustment <- case_when(
    grepl("true", filename) ~ 0,
    grepl("test", filename) ~ 1,
    grepl("valid", filename) ~ 2
  )
  
  variables <- ifelse(grepl("[0-9]{4}$", variables), 
                      paste0(gsub("[0-9]{4}", "", variables), 
                             as.numeric(stringr::str_extract(variables, "[0-9]{4}")) - adjustment), 
                      variables)
  
  count_xwalk <- count_xwalk %>%
    mutate(year = as.numeric(stringr::str_extract(sample_var, "[0-9]{4}")) - adjustment) %>%
    mutate_at(c("Variable", "sample_var"), ~stringr::str_replace(., "[0-9]{4}$", as.character(year))) %>%
    select(-year) 

  list(variables, count_xwalk)
}

read_data <- function(filename, state_names = "ZZ") {
  full_data <- read_rds(filename)[[1]] %>%
    ungroup() %>%
    arrange(state, cpuma) %>%
    filter(!grepl(state_names, state))
  
  replicate_data <- read_rds(filename)[-1] %>%
    map(~ungroup(.x) %>% arrange(state, cpuma) %>% filter(!grepl(state_names, state)))
  
  treatment_data <- subset(full_data, treatment == 1) %>%
    arrange(state, cpuma) %>%
    filter(!grepl(state_names, state))
  
  control_data <- subset(full_data, treatment == 0) %>%
    arrange(state, cpuma) %>%
    filter(!grepl(state_names, state))
  
  total_data <- full_data %>%
    filter(treatment < 2) %>%  
    arrange(treatment, state, cpuma) %>%
    filter(!grepl(state_names, state))
  
  replicate_treatment <- map(replicate_data, ~subset(.x, treatment == 1) %>% arrange(state, cpuma))
  replicate_control <- map(replicate_data, ~subset(.x, treatment == 0) %>% arrange(state, cpuma))
  
  list(treatment_data = treatment_data,
       control_data = control_data,
       replicate_treatment = replicate_treatment,
       replicate_control = replicate_control,
       total_data = total_data)
}

# generate individual level covariate estimates ----------------------------------------------------
sigma_uu_i_estimates <- function(data_list, variable_list, filename, suffix2 = "") {
  suffix <- stringr::str_extract(filename, "[a-z]{1,5}\\.rds")
  
  treatment_data <- data_list$treatment_data
  control_data <- data_list$control_data
  replicate_treatment <- data_list$replicate_treatment
  replicate_control <- data_list$replicate_control
  
  variables <- variable_list[[1]]
  
  treat_uu_i <- map(1:nrow(treatment_data), ~generate_full_cov_matrix(treatment_data, 
                                                                      replicate_treatment, 
                                                                      variables, .x))
  
  control_uu_i <- map(1:nrow(control_data), ~generate_full_cov_matrix(control_data, 
                                                                      replicate_control, 
                                                                      variables, .x))
  
  saveRDS(treat_uu_i, paste0("../01_ProcessedData/sigma-uu-i-txonly-", suffix2, suffix))
  saveRDS(control_uu_i, paste0("../01_ProcessedData/sigma-uu-i-ctonly-", suffix2, suffix))
  
  list(treat_uu_i = treat_uu_i, control_uu_i = control_uu_i)
}

finalize_data <- function(data_list, sigma_uu_i, variable_list, output_name) {
  treatment_data <- data_list$treatment_data
  control_data <- data_list$control_data
  treat_uu_i <- sigma_uu_i$treat_uu_i
  control_uu_i <- sigma_uu_i$control_uu_i
  
  variables <- variable_list[[1]]; count_xwalk <- variable_list[[2]]
  
  # generate error covariance matrix estimates -------------------------------------------------------
  treat_uu <- generate_all_covariance_estimates(treatment_data, treat_uu_i, count_xwalk, variables)
  control_uu <- generate_all_covariance_estimates(control_data, control_uu_i, count_xwalk, variables)
  
  # calculate kappa for each covariance estimator ----------------------------------------------------
  kappa_treat <- calculate_kappa_all(treatment_data, treat_uu, variables)
  kappa_control <- calculate_kappa_all(control_data, control_uu, variables)
  
  # impute data using kappa estimates ----------------------------------------------------------------
  tdat_trans_tx <- transform_data(treatment_data, variables, kappa_treat)
  cdat_trans_ct <- transform_data(control_data, variables, kappa_control)
  
  # final files ---------------------------------------------------------------------------------------
  all_data <- tdat_trans_tx %>%
    mutate(treatment = 1) %>% 
    rbind( 
      cdat_trans_ct %>%
        mutate(treatment = 0) 
    ) %>%
    mutate(data = map(data, ~arrange(.x, state, cpuma))) 
  
  all_data
}

variables <- read_csv('../02_Specs/tol_specs.csv') %>%
  mutate(`Reference Variable` = case_when(grepl('child', Variable) ~ 0, TRUE ~ `Reference Variable`)) %>%
  filter(`Reference Variable` == 0) %>%
  .$Variable %>%
  sort()

c2_indices <-readRDS("../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-true.rds")[[1]] %>%
  filter(treatment == 1) %>%
  arrange(state, cpuma) %>%
  mutate(c2_remove = ifelse(state %in% c("CA", "CT", "MN", "NJ", "WA"), FALSE, TRUE)) %>%
  .$c2_remove

true_file  <- "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-true.rds"
test_file  <- "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-test.rds"
valid_file <- "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-valid.rds"

varlist_true  <- modify_variables(true_file, variables)
varlist_test  <- modify_variables(test_file, variables)
varlist_valid <- modify_variables(valid_file, variables)

data_true <- read_data(true_file) 
data_true_c2 <- read_data(true_file, c("CA|CT|MN|NJ|WA"))
data_test <- read_data(test_file)
data_test_c2 <- read_data(test_file, c("CA|CT|MN|NJ|WA"))
data_valid <- read_data(valid_file)
data_valid_c2 <- read_data(valid_file, c("CA|CT|MN|NJ|WA"))

sigma_uu_i_true <- sigma_uu_i_estimates(data_true, varlist_true, true_file)
sigma_uu_i_true_c2 <- sigma_uu_i_estimates(data_true_c2, varlist_true, true_file, "-c2")
sigma_uu_i_test <- sigma_uu_i_estimates(data_test, varlist_test, test_file)
sigma_uu_i_test_c2 <- sigma_uu_i_estimates(data_true_c2, varlist_test, test_file, "-c2")
sigma_uu_i_valid <- sigma_uu_i_estimates(data_valid, varlist_valid, valid_file)
sigma_uu_i_valid_c2 <- sigma_uu_i_estimates(data_true_c2, varlist_valid, valid_file, "-c2")

sigma_uu_i_true <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly-true.rds"),
                        control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly-true.rds"))

sigma_uu_i_true_c2 <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly--c2true.rds"),
                           control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly--c2true.rds"))

sigma_uu_i_test <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly-test.rds"),
                        control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly-test.rds"))

sigma_uu_i_test_c2 <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly--c2test.rds"),
                           control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly--c2test.rds"))

sigma_uu_i_valid <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly-valid.rds"),
                         control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly-valid.rds"))

sigma_uu_i_valid_c2 <- list(treat_uu_i   = readRDS("../01_ProcessedData/sigma-uu-i-txonly--c2valid.rds"),
                            control_uu_i = readRDS("../01_ProcessedData/sigma-uu-i-ctonly--c2valid.rds"))

all_true <- finalize_data(data_true, sigma_uu_i_true, varlist_true)
all_true_c2 <- finalize_data(data_true_c2, sigma_uu_i_true_c2, varlist_true)
all_test <- finalize_data(data_test, sigma_uu_i_test, varlist_test)
all_test_c2 <- finalize_data(data_test_c2, sigma_uu_i_test_c2, varlist_test)
all_valid <- finalize_data(data_valid, sigma_uu_i_valid, varlist_valid)
all_valid_c2 <- finalize_data(data_valid_c2, sigma_uu_i_valid_c2, varlist_valid)

all_data <- bind_rows(
  all_true %>% mutate(set = "true"),
  all_test %>% mutate(set = "test"),
  all_valid %>% mutate(set = "valid")
) %>%
  mutate(data = map(data, ~filter(.x, state != "NH")))

saveRDS(all_data, "../01_ProcessedData/calibrated-data-all.rds")

all_data_c2 <- bind_rows(
  all_true_c2 %>% mutate(set = "true"),
  all_test_c2 %>% mutate(set = "test"),
  all_valid_c2 %>% mutate(set = "valid")
) %>%
  mutate(data = map(data, ~filter(.x, state != "NH")))

saveRDS(all_data_c2, "../01_ProcessedData/calibrated-data-c2.rds")

all_data    <- readRDS("../01_ProcessedData/calibrated-data-all.rds")
all_data_c2 <- readRDS("../01_ProcessedData/calibrated-data-c2.rds")

######################################################################################################
##################################### LEAVE ONE OUT DATASETS #########################################
######################################################################################################
state_list <- append(list(c("")), as.list(unique(all_data$data[[1]]$state))) 
names(state_list) <- c("Preferred", unique(all_data$data[[1]]$state))
data_list <- map(state_list, ~data_subset(all_data$data, .x))
saveRDS(data_list, "../01_ProcessedData/calibrated-data-all-subsets-c1-all.rds")

state_list <- append(list(c("")), as.list(unique(all_data_c2$data[[1]]$state)))
names(state_list) <- c("Early expansion", unique(all_data_c2$data[[1]]$state))
data_list <- map(state_list, ~data_subset(all_data_c2$data, .x))
saveRDS(data_list, "../01_ProcessedData/calibrated-data-all-subsets-c2-all.rds")

######################################################################################################
##################################### CORRELATED ADJUSTMENT ##########################################
######################################################################################################

# correlated adjustment test: model correlations among the data?
xhat.cor <- transform_correlated_data(data_true$treatment_data, 
                                      sigma_uu_i_true$treat_uu_i,
                                      variables)

diag(cov(xhat.cor[, variables])) - diag(cov(data_true$treatment_data[, variables])) #unsatisfactory
