# program: hte-oate-analysis.R
# purpose: calculate HTEs
# author: max rubinstein
# date modified: january 25, 2021

library(tidyverse)
library(clubSandwich)

# read variable names
read_varnames <- function(group = "") {
  variable_names <- read_csv('../02_Specs/tol_specs.csv') %>%
    filter(`Reference Variable` == 0) %>%
    filter(!Group == group) %>%
    arrange(Variable) %>%
    .$Variable
  variable_names
}

data_subset <- function(data_list, state_name) {
  map(data_list, ~filter(.x, !state %in% state_name))
}

calc_hte <- function(contrast, data, model, alpha) {
  lincom <- rep(0, length(coef(model)[!is.na(coef(model))]))
  ivars <- grep("treatment:repub", names(coef(model)[!is.na(coef(model))]))
  lincom[ivars] <- rep(contrast, length(ivars))
  dhat <- lincom %*% coef(model)[!is.na(coef(model))]
  dhat
}

# run analysis --------------------------------------------------------------------------------
merged_data_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.rds") %>%
  filter(set == "true") %>%
  select(-treatment) %>%
  unnest(cols = c(data)) %>%
  nest(-key) %>%
  mutate(data = map(data, ~arrange(.x, state, cpuma)))

merged_data_c2 <- readRDS("../01_ProcessedData/calibrated-data-c2.rds") %>%
  select(-treatment) %>%
  unnest(cols = c(data)) %>%
  nest(-key) %>%
  mutate(data = map(data, ~arrange(.x, state, cpuma)))

variables <- map(c("", 1:5), ~read_varnames(.x))[[1]]

full_jackknife_dat_c1 <- readRDS("../01_ProcessedData/full-jackknife-data-c1.rds")
full_jackknife_dat_c2 <- readRDS("../01_ProcessedData/full-jackknife-data-c2.rds")

treatvar <- paste0("treatment:", variables)

formula.all <- paste0("hins_unins_pct_2014 ~ treatment + ", paste0(c(variables, treatvar), collapse = "+"), collapse = "")

gen_ests <- function(data_list) {
  d1 <- data_list[[1]]
  d2 <- data_list[[2]]
  d3 <- data_list[[3]]

  m1a <- lm(as.formula(formula.all), data = d1)
  m2a <- lm(as.formula(formula.all), data = d2)
  m3a <- lm(as.formula(formula.all), data = d3)

  table95 <- map2(list(m1a, m2a, m3a), 
                  list(d1, d1, d1), ~calc_hte(-50, .y, .x, 0.05)) %>%
    invoke(rbind, .) %>%
    as_tibble()  
  table95
}

ests_c1 <- gen_ests(merged_data_c1$data) %>%
  rename(ests = V1) %>%
  mutate(sigma_estimate = c("Heterogeneous", "Homogeneous", "None"))

ests_c2 <- gen_ests(merged_data_c2$data) %>%
  rename(ests = V1) %>%
  mutate(sigma_estimate = c("Heterogeneous", "Homogeneous", "None"))

se_c1 <- map(full_jackknife_dat_c1, gen_ests) %>%
  map(~mutate(.x, sigma_estimate = c("Heterogeneous", "Homogeneous", "None"))) %>%
  invoke(rbind, .) 

se_c2 <- map(full_jackknife_dat_c2, gen_ests) %>%
  map(~mutate(.x, sigma_estimate = c("Heterogeneous", "Homogeneous", "None"))) %>%
  invoke(rbind, .) 

proc_rslt <- function(ests, se_ests) {
  ests %>%
    group_by(sigma_estimate) %>%
    summarize_all(list(mean = mean)) %>%
    right_join(se_ests, by = "sigma_estimate") %>%
    mutate(se = (V1 - mean)^2) %>%
    group_by(sigma_estimate) %>%
    summarize_at(vars(contains("se")), ~sqrt(((n(.) - 1)/n(.))sum(.))) %>%
    left_join(ests) %>%
    mutate(l95ci = ests - 1.96*se, u95ci = ests + 1.96*se)
}

final <- proc_rslt(ests_c1, se_c1) %>%
  mutate(Dataset = "Primary") %>%
  bind_rows(
    proc_rslt(ests_c2, se_c2) %>%
      mutate(Dataset = "Early expansion")
  ) %>%
  rename(Adjustment = sigma_estimate, Estimate = ests) %>%
  mutate(`CI (95%)` = paste0("(", round(l95ci, 2), ", ", round(u95ci, 2), ")")) %>%
  select(Estimate, Adjustment, `CI (95%)`, Dataset) 

print(xtable::xtable(final), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)
