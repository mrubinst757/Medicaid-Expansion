# program: balance-plots.R
# purpose: analyze covariate balance for preferred weighting estimators
# author: max rubinstein
# date modified: december 22, 2020

# load libraries 
library(tidyverse)

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
}

# HSBW imbalances ---------------------------------------------------------------------------------------
hsbw_baltab <- function(data, weights_etu, variables) {
  hsbw_balance <- calc_hsbw_balance(data, weights_etu, variables)
  
  control_data <- subset(data, treatment == 0) 
  treatment_data <- subset(data, treatment == 1) 
  tx_mean <- colMeans(treatment_data[variables])
  targets <- colMeans(control_data[variables])
  var1 <- diag(cov(treatment_data[variables]))#/nrow(treatment_data)
  var0 <- diag(cov(control_data[variables]))#/nrow(control_data)
  
  tibble(
    target = targets,
    wtx_mean = hsbw_balance$V1,
    variable = variables,
    tx_mean = tx_mean
  ) %>%
    mutate(`Unweighted Diff` = tx_mean - target,
           `Unweighted SMD` = (tx_mean - target)/sqrt(var1 + var0),
           `Weighted Diff` = wtx_mean - target,
           `Weighted SMD` = `Weighted Diff`/sqrt(var1 + var0)) %>%
    select(Variables = variable, contains("Diff"), contains("SMD")) 
}

hsbw_balplot <- function(hsbw_balance, file_extension, diff_criterion,
                         SMD = "Pct") {
  if (SMD != "SMD") {
    plot <- hsbw_balance %>%
      filter(abs(`Unweighted Diff`) > diff_criterion | abs(`Weighted Diff`) > diff_criterion) %>%
      rename(`Unweighted` = `Unweighted Diff`,
             `Weighted` = `Weighted Diff`) %>%
      select(-contains("SMD")) %>%
      gather(key, value, -Variables) %>%
      mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
      mutate_at("key", ~factor(., levels = c("Weighted", "Unweighted"))) %>%
      ggplot(aes(x = reorder(Variables, value), y = value, fill = rev(key))) +
      geom_bar(stat = "identity", position = "dodge")+
      #scale_fill_brewer(palette = "Set1", name = "Difference") +
      scale_fill_grey(name = "Difference") +
      coord_flip() +
      theme_minimal() +
      ylab("Treated, weighted treated means, \n minus untreated means") +
      xlab("Variables") 
  }
  
  if (SMD == "SMD") {
    plot <- hsbw_balance %>%
      filter(abs(`Unweighted Diff`) > diff_criterion | abs(`Weighted Diff`) > diff_criterion) %>%
      rename(`Unweighted` = `Unweighted SMD`,
             `Weighted` = `Weighted SMD`) %>%
      select(-contains("Diff")) %>%
      gather(key, value, -Variables) %>%
      mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
      mutate_at("key", ~factor(., levels = c("Weighted", "Unweighted"))) %>%
      ggplot(aes(x = reorder(Variables, value), y = value, fill = rev(key))) +
      geom_bar(stat = "identity", position = "dodge")+
      #scale_fill_brewer(palette = "Set1", name = "Difference") +
      scale_fill_grey(name = "Difference") +
      coord_flip() +
      theme_minimal() +
      ylab("Standardized Mean Difference (treatment minus control)") +
      xlab("Variables") 
  }
  ggsave(paste0("../../02_Paper/01_Plots/balance-plot-etu", file_extension, "-", SMD, ".png"), plot)
  return(plot)
}

hsbw_balplot_smd <- function(results, file_extension) {
  plot <- results %>%
    nest(`Standardized mean difference` = c("Unweighted SMD", "Weighted SMD"), 
         `Percentage point difference` = c("Unweighted Diff", "Weighted Diff"))  %>%
    gather(key, value, -Variables) %>%
    mutate_at("value", ~map(., ~set_names(.x, c("Unweighted", "Weighted")))) %>%
    unnest(cols = c(value)) %>%
    mutate_at("Variables", ~stringr::str_replace_all(., var_names)) %>%
    gather(weighted, difference, Unweighted, Weighted) %>%
    mutate_at("weighted", ~factor(., levels = c("Weighted", "Unweighted"))) %>%
    group_by(Variables) %>%
    mutate(pctdiff = difference[key == "Standardized mean difference" & weighted == "Unweighted"]) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(Variables, abs(pctdiff)), y = difference, fill = weighted)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~key, scales = "free_x") +
    theme_bw() +
    scale_fill_grey(name = "Difference", breaks = c("Unweighted", "Weighted")) +
    coord_flip() +
    ylab("Mean Difference (treatment mean minus control mean)") +
    xlab("Variables \n (ordered by absolute magnitude \n of unweighted percentage point imbalance)")
  
  ggsave(paste0("../../02_Paper/01_Plots/balance-plot-all-etu", file_extension, ".png"), plot,
         width = 12, height = 7)
  return(plot)
}


