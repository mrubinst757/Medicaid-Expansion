# program: 02-estimate-models.R
# purpose: generate weighting estimators using various data subsets
# author: max rubinstein
# date modified: december 14, 2020

# test setup ------------------------------------------------------------------------------------
source("03_Analysis/02_model-estimation.R")

all_data_c1 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c1-all.rds") %>%
  map(~.x[1:3])

all_data_c2 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c2-all.rds") %>%
  map(~.x[1:3])

cdata <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c1-all.rds")[[1]][[6]]

variables <- read_csv("../02_Specs/tol_specs.csv") %>%
  filter(`Reference Variable` == 0) %>%
  arrange(Variable) %>%
  .$Variable

targets <- colMeans(cdata[variables])

tol_list <- read_csv("../02_Specs/tol_specs.csv") 

tol_list <- map(0:5, ~tol_list %>% mutate(`Base Tol` = if_else(grepl(.x, Group), 100, `Base Tol`))) %>%
  map(~filter(.x, Variable %in% variables)) %>%
  map(~arrange(.x, Variable)) %>%
  map(~.x$`Base Tol`) %>%
  map(~set_names(.x, variables))

cov_groups <- c("None", "Republican", "Unins & Unemp", "Urb-Age-Educ-Cit-Mar-Stu-Dis-F",
                "Race-Eth-For-Inc-Pov", "Child-PGrowth-HHRatio")

model_names <- c("SBW", "H-SBW", "BC-SBW", "BC-HSBW")

cov_models <- c("sigma_uu_i_modeled", "sigma_uu_avg", "sigma_zero")

# helper fun to set names of list of list of lists
set_cov_mod_name <- function(result_list) {
  map(result_list, ~set_names(.x, cov_models))
}

data_subset <- function(data_list, state_name) {
  map(data_list, ~filter(.x, !state %in% state_name))
}

# calculate results --------------------------------------------------------------------------------------
bw_models <- function(c1_data, c2_data, tol_list, targets, stop, file_extension) {
  c1_results <- map(c1_data, ~iterate_covariate_subsets(.x, tol_list, targets,
                                                                    distance = linf_imbalance, stop_criterion = stop))
  c1_results <- map(c1_results, ~set_names(.x, cov_groups))
  c1_results <- map(1:length(c1_data), ~map(c1_results[[.x]], ~set_names(.x, model_names)))
  c1_results <- map(1:length(c1_data), ~map(c1_results[[.x]], set_cov_mod_name))
  names(c1_results) <- names(c1_data)
  saveRDS(c1_results, paste0("../04_Output/c1-results", file_extension, ".rds"))
  
  c2_results <- map(c2_data, ~iterate_covariate_subsets(.x, tol_list, targets,
                                                        distance = linf_imbalance, stop_criterion = stop))
  c2_results <- map(c2_results, ~set_names(.x, cov_groups))
  c2_results <- map(1:length(c2_data), ~map(c2_results[[.x]], ~set_names(.x, model_names)))
  c2_results <- map(1:length(c2_data), ~map(c2_results[[.x]], set_cov_mod_name))
  names(c2_results) <- names(c2_data)
  saveRDS(c2_results, paste0("../04_Output/c2-results", file_extension, ".rds"))
}

bw_models(all_data_c1, all_data_c2, tol_list, targets, stop = 0.5, "")
