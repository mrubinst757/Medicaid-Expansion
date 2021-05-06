# program: oate-analysis
# purpose: analyze overlap weights
# author: max rubinstein
# date modified: december 14, 2020

# load libraries and read data --------------------------------------------------------------
library(tidyverse)
library(assertthat)
library(scales)

# analyze results ------------------------------------------------------------------------------------
sigma_estimator <- c("sigma_uu_i_modeled", "sigma_uu_avg", "sigma_zero")

process_results <- function(results) {
  map(results, ~invoke(rbind, .x)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    mutate(txfx = unlist(txfx)) %>%
    mutate(sigma_estimate = rep(sigma_estimator, 6)) %>%
    mutate(variables_subset = rep(covariate_group, each = 5)) 
}

loo_covariates_plot <- function(results) {
  process_results(results) %>%
    filter(sigma_estimate != "sigma_avg") %>%
    mutate_at("variables_subset", 
              ~stringr::str_replace_all(., c("Urb-Age-Educ-Cit-Mar-Stu-Dis-F" = "Urb-Age-Educ-Cit\n-Mar-Stu-Dis-F",
                                             "Child-PGrowth-HHRatio" = "Child-PGrowth-\nHHRatio",
                                             "Race-Eth-For-Inc-Pov" = "Race-Eth-\nFor-Inc-Pov"))) %>%
    mutate_at("variables_subset", ~factor(.x, levels = c(
      "None", "Republican", "Unins & Unemp", "Urb-Age-Educ-Cit\n-Mar-Stu-Dis-F",
      "Race-Eth-\nFor-Inc-Pov", "Child-PGrowth-\nHHRatio"
    ))) %>%
    ggplot(aes(x = variables_subset, y = txfx, fill = sigma_estimate)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1", name = "Sigma Estimator") +
    xlab("Variables Removed") +
    ylab("Estimated Contrast") 
}

variables = read_csv('../02_Specs/tol_specs.csv') %>%
  mutate(`Reference Variable` = case_when(grepl('child', Variable) ~ 0, TRUE ~ `Reference Variable`)) %>%
  filter(`Reference Variable` == 0) %>%
  .$Variable %>%
  sort()

orig_data <- readRDS("../01_ProcessedData/cpuma-analytic-file-2011-2014-r0-r80.rds")[[1]] %>%
  ungroup() %>%
  select(state, treatment) %>%
  distinct() %>%
  filter(treatment < 2)

oate_results_c1 <- readRDS("../04_Output/oate-results-c1.rds")
oate_results_c2 <- readRDS("../04_Output/oate-results-c2.rds")

state_names_c1 <- names(oate_results_c1)[-1]
state_names_c2 <- names(oate_results_c2)[-1]

jackknife_c1 <- readRDS("../04_Output/oate-jackknife-c1.rds")
jackknife_c2 <- readRDS("../04_Output/oate-jackknife-c2.rds") %>%
  map(~set_names(.x, state_names_c2))

covariate_group <- c('None', 'Republican', "Unins & Unemp", "Urb-Age-Educ-Cit-Mar-Stu-Dis-F",
                     "Race-Eth-For-Inc-Pov", "Child-PGrowth-HHRatio")

assert_that(identical(names(jackknife_c1$weights), names(oate_results_c1)[-c(1)]))
assert_that(identical(names(jackknife_c2$weights), names(oate_results_c2)[-c(1)]))

# calculate standard errors ---------------------------------------------------------------------
num_states_c1 <- length(jackknife_c1$data)
num_states_c2 <- length(jackknife_c2$data)

c1_txfx <- process_results(oate_results_c1$Preferred)
c2_txfx <- process_results(oate_results_c2$`Early Expansion`)

c1_jackknife_states <- map(oate_results_c1[-1], ~process_results(.x)$txfx)
c2_jackknife_states <- map(oate_results_c2[-1], ~process_results(.x)$txfx)
c1_jackknife_proc <- map(jackknife_c1$weights, ~process_results(.x)$txfx)
c2_jackknife_proc <- map(jackknife_c2$weights, ~process_results(.x)$txfx)

c1_mean_states <- Reduce(`+`, c1_jackknife_states)/length(jackknife_c1$data)
c2_mean_states <- Reduce(`+`, c2_jackknife_states)/length(jackknife_c2$data)
c1_mean_proc   <- Reduce(`+`, c1_jackknife_proc)/length(jackknife_c1$data)
c2_mean_proc   <- Reduce(`+`, c2_jackknife_proc)/length(jackknife_c2$data)

seest_c1_states <- sqrt((num_states_c1 - 1)/(num_states_c1)*Reduce(`+`, map(c1_jackknife_states, ~(.x - c1_mean_states)^2)))
seest_c2_states <- sqrt((num_states_c2 - 1)/(num_states_c2)*Reduce(`+`, map(c2_jackknife_states, ~(.x - c2_mean_states)^2)))
seest_c1_proc <- sqrt((num_states_c1 - 1)/(num_states_c1)*Reduce(`+`, map(c1_jackknife_proc, ~(.x - c1_mean_proc)^2)))
seest_c2_proc <- sqrt((num_states_c2 - 1)/(num_states_c2)*Reduce(`+`, map(c2_jackknife_proc, ~(.x - c2_mean_proc)^2)))

# create tables ---------------------------------------------------------------------------------------
create_result_table <- function(txfx, sehat_states, sehat_proc) {
  txfx %>%
    mutate(sehat1 = sehat_states, 
           l95ci1 = txfx - 1.96*sehat1, 
           u95ci1 = txfx + 1.96*sehat1,
           sehat2 = sehat_proc, 
           l95ci2 = txfx - 1.96*sehat2, 
           u95ci2 = txfx + 1.96*sehat2)
}

confint_table <- function(result_table) {
  result_table %>%
    select(-weights, -subset) %>%
    mutate(ci_data = paste0("(", round(l95ci1, 2), ", ", round(u95ci1, 2), ")"),
           ci_boot = paste0("(", round(l95ci2, 2), ", ", round(u95ci2, 2), ")")) %>%
    select(-contains('sehat'), -contains("95")) %>%
    mutate_at("variables_subset", ~factor(., levels = covariate_group)) %>%
    mutate_at("variables_subset", ~as.numeric(.) - 1) %>%
    mutate_at("sigma_estimate", ~factor(., levels = sigma_estimator)) %>%
    arrange(variables_subset, sigma_estimate) 
}

# are loo estimates centered at original value?
process_loo_data <- function(primary_results, loo_results, covariate_group) {
  loo_results %>%
    invoke(cbind, .) %>%
    as.data.frame() %>%
    mutate(main_estimate = primary_results$txfx) %>%
    mutate(sigma_estimate = rep(sigma_estimator, 6)) %>%
    mutate(variables_subset = rep(covariate_group, each = length(sigma_estimator))) %>%
    gather(key, value, -main_estimate, -sigma_estimate, -variables_subset) 
}

centered <- function(primary_results, loo_results, covariate_group) {
  process_loo_data(primary_results, loo_results, covariate_group) %>%
    group_by(variables_subset, sigma_estimate) %>%
    mutate(mean_jackknife = mean(value),
           diff = mean_jackknife - main_estimate) %>%
    select(-value, -key) %>%
    ungroup() %>%
    distinct() %>%
    arrange(-abs(diff))
}


