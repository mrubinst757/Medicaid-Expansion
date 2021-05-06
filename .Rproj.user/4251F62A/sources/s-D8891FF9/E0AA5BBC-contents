# program: balance-plots.R
# purpose: analyze covariate balance for preferred weighting estimators
# author: max rubinstein
# date modified: december 22, 2020

# load libraries 
library(tidyverse)

# calculate oate region means
calc_oate_region <- function(merged_data, weights_oate, variables) {
  merged_data %>%
    ungroup() %>%
    arrange(state, cpuma) %>%
    mutate(weights_oate = weights_oate) %>%
    select(all_of(variables), weights_oate, treatment) %>%
    gather(key, value, -weights_oate, -treatment) %>%
    group_by(key, treatment) %>%
    dplyr::summarize(value = sum(value*weights_oate)/sum(weights_oate)) %>%
    filter(treatment == 1) %>%
    select(key, value)
}

# calculate weighted balance for etc
calc_hsbw_balance <- function(merged_data, weights_etu, variables) {
  merged_data %>%
    ungroup() %>%
    filter(treatment == 1) %>%
    arrange(state, cpuma) %>%
    mutate(weights_etu = weights_etu) %>%
    summarize_at(all_of(variables), ~sum(.*weights_etu)/sum(weights_etu)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column()
#    select(all_of(variables), weights_etu, treatment) %>%
#    gather(key, value, -weights_etu, -treatment) %>%
#    mutate_at("key", ~stringr::str_replace_all(., var_names)) %>%
#    group_by(key, treatment) %>%
#    dplyr::summarize(value = sum(value*weights_etu)/sum(weights_etu)) %>%  
#    select(key, value)
}

# calculate covariate means for treatment and control data
calc_balance_data <- function(merged_data, weights_oate, weights_etu, variables) {
  oate_region <- calc_oate_region(merged_data, weights_oate, variables)
  hsbw_balance <- calc_hsbw_balance(merged_data, weights_etu, variables)
  
  control_data <- subset(merged_data, treatment == 0)
  treatment_data <- subset(merged_data, treatment == 1)
  tx_mean <- colMeans(treatment_data[variables])
  ct_mean <- colMeans(control_data[variables])
  
  tibble(
    ct_mean = ct_mean,
    tx_mean = tx_mean,
    wtx_mean = hsbw_balance$V1,
    oate_mean = oate_region$value,
    variable = variables
  ) %>%
    gather(key, value, -variable) 
}

# OATE region comparison --------------------------------------------------------------------------------
calc_oate_distance <- function(merged_data, weights_oate, weights_etu, variables) {
  merged_data %>%
    arrange(state, cpuma) %>%
    calc_balance_data(weights_oate, weights_etu, variables) %>%
    filter(key != "wtx_mean") %>%
    mutate_at("variable", ~stringr::str_replace_all(., var_names)) %>%
    group_by(variable) %>%
    mutate(dist_ct = value[key == "oate_mean"] - value[key == "ct_mean"],
           dist_tx = value[key == "oate_mean"] - value[key == "tx_mean"]) %>%
    distinct(variable, dist_ct, dist_tx) 
}

calc_mean_oate_dist <- function(oate_distance) {
  oate_distance %>%
    ungroup() %>%
    dplyr::summarize(dist_ct = mean(abs(dist_ct)), dist_tx = mean(abs(dist_tx)))
}

oate_distance_plot <- function(oate_distance, file_extension) {
  oate_distance %>%
    gather(key, value, -variable) %>%
    filter(abs(value[key == "dist_ct"]) > 1 & abs(value[key == "dist_tx"]) > 1) %>%
    mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
    ggplot(aes(x = reorder(variable, -value), y = value, fill = key)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    theme_minimal() +
    ylab("Difference between overlap region and treated, untreated regions") +
    xlab("Variables") +
    scale_fill_brewer(palette = "Set1", name = "Region") +
    ggsave(paste0("../../02_Paper/01_Plots/oate-imbalances-", file_extension, ".png"))
}

# HSBW imbalances ---------------------------------------------------------------------------------------
hsbw_baltab <- function(data, weights_etu, variables) {
  hsbw_balance <- calc_hsbw_balance(data, weights_etu, variables)
  
  control_data <- subset(data, treatment == 0) 
  treatment_data <- subset(data, treatment == 1) 
  tx_mean <- colMeans(treatment_data[variables])
  targets <- colMeans(control_data[variables])
  
  tibble(
    target = targets,
    wtx_mean = hsbw_balance$V1,
    variable = variables,
    tx_mean = tx_mean
  ) %>%
    mutate(`Unweighted Diff` = tx_mean - target,
           `Weighted Diff` = wtx_mean - target) %>%
    select(Variables = variable, contains("Diff")) 
  
}

hsbw_balplot <- function(hsbw_balance, file_extension) {
  hsbw_balance %>%
    filter(abs(`Unweighted Diff`) > 1 | abs(`Weighted Diff`) > 1) %>%
    gather(key, value, -Variables) %>%
    mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
    mutate_at("key", ~factor(., levels = c("Weighted Diff", "Unweighted Diff"))) %>%
    ggplot(aes(x = reorder(Variables, value), y = value, fill = key)) +
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_brewer(palette = "Set1", name = "Difference") +
    coord_flip() +
    theme_minimal() +
    ylab("Treated, weighted treated means, \n minus untreated means") +
    xlab("Variables") +
    ggsave(paste0("../../02_Paper/01_Plots/balance-plot-etu", file_extension, ".png"))
}

# other OATE plots --------------------------------------------------------------------------
oate_region_table <- function(data, w1, w2, w3) {
  map2(data, list(w1, w2, w3), ~.x %>%
         arrange(state, cpuma) %>%
         mutate(weights = .y) %>%
         group_by(treatment, state) %>%
         dplyr::summarize(weights = sum(weights))
  ) %>%
    map2(c("sigma_uu_i", "sigma_uu_avg", "sigma_zero"),
         ~mutate(.x, `Sigma estimate` = .y)) %>%
    invoke(rbind, .) %>%
    group_by(`Sigma estimate`, treatment) %>%
    mutate(weights = 100*weights/sum(weights)) %>%
    spread(`Sigma estimate`, weights) %>%
    arrange(-sigma_uu_i) %>%
    select(State = state, Treatment = treatment, 
           `Heterogeneous adjustment` = sigma_uu_i, 
           `Homogeneous adjustment` = sigma_uu_avg, 
           `No adjustment` = sigma_zero)
}

oate_region_plot <- function(data, weights, file_extension) {
  data %>%
    arrange(state, cpuma) %>%
    mutate(weights = weights) %>%
    group_by(treatment, state) %>%
    dplyr::summarize(weights = sum(weights)) %>%
    ungroup() %>%
    group_by(treatment) %>%
    mutate(weights = 100*weights/sum(weights)) %>%
    ggplot(aes(x = reorder(state, -weights), y = weights, fill = factor(treatment))) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    xlab("State") +
    scale_fill_brewer(palette = "Set1", name = "Expansion") +
    ylab("Sum of Weights within State") +
    xlab("State") +
    scale_x_discrete(guide = guide_axis(n.dodge=2))+
    ggsave(paste0("../../02_Paper/01_Plots/oate-region-", file_extension, ".png"))
}
