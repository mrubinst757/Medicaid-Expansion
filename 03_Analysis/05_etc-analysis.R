# program: etc-analysis
# purpose: calculate treatment effects and variance estimates for ETC
# author: max rubinstein
# date modified: december 28, 2020

# load libraries and initialize functions ---------------------------------------------------------
library(tidyverse)
library(scales)
library(xtable)
library(clubSandwich)
library(sandwich)
library(lmtest)

# calculate treatment effects 
calculate_txfx <- function(data_list, weight_list, y0est) {
  
  outer_product <- function(vector) {
    vector %*% t(vector)
  }
  
  estimates_per_weight_type <- function(data_list, y0est, variables, weights) {
    
    estimates <- map2(data_list, weights, ~sum(.x$hins_unins_pct_2014*.y$weights)/sum(.y$weights) - y0est) %>%
      unlist() %>%
      set_names(names(weights))
    
    estimates
  }
  
  estimates_per_cov_subset <- function(data_list, y0est, variables, weight_list) {
    
    map(weight_list, ~estimates_per_weight_type(data_list, y0est, variables, .x)) %>%
      set_names(names(weight_list)) %>%
      invoke(rbind, .) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      gather(sigma_estimate, psihat, -rowname) 
  }
  
  variables <- map(weight_list, ~.x[[1]][[1]]$tols) %>% 
    map(~names(.x)[.x != 100])
  
  map2(weight_list, variables, ~estimates_per_cov_subset(data_list, y0est, .y, .x)) %>%
    set_names(names(weight_list)) %>%
    invoke(rbind, .) %>%
    rename(weight_type = rowname) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(variable_subset = rowname) %>%
    mutate_at("variable_subset", ~gsub("\\.[0-9]?[0-9]", "", .)) 
}

# read variable names
read_varnames <- function(group = "") {
  variable_names <- read_csv("../02_Specs/tol_specs.csv") %>%
    filter(`Reference Variable` == 0) %>%
    filter(!Group == group) %>%
    arrange(Variable) %>%
    .$Variable
  variable_names
}

# state by state sensitivity analysis
sensitivity <- function(primary_result, loo_states, number)  {
  primary_result %>%
    mutate(
      loostate = loo_states[[number]]$psihat,
      diff = loostate - psihat
    ) %>%
    select(variable_subset, weight_type, sigma_estimate, psihat, loostate, diff) 
}

# full result tables -------------------------------------------------------------------------------
calc_all_txfx <- function(c1_data, c1_results, c2_data, c2_results, 
                          c1_jackknife, c2_jackknife, y0est, y0varest) {
  
  # generate estimate for observed outcomes among controls
  variable_list <- map(c("", 1:5), read_varnames) 
  
  # calculate treatment effects ----------------------------------------------------------------------
  txfx_c1 <- calculate_txfx(c1_data[[1]], c1_results[[1]], y0est) 
  txfx_c2 <- calculate_txfx(c2_data[[1]], c2_results[[1]], y0est) 
  
  loo_states_c1 <- map2(c1_data[-1], c1_results[-1], ~calculate_txfx(.x, .y, y0est))
  loo_states_c2 <- map2(c2_data[-1], c2_results[-1], ~calculate_txfx(.x, .y, y0est))
  
  loo_proc_c1 <- map2(c1_jackknife$data, c1_jackknife$weights, ~calculate_txfx(.x, .y, y0est))
  loo_proc_c2 <- map2(c2_jackknife$data, c2_jackknife$weights, ~calculate_txfx(.x, .y, y0est))
  
  # variance estimation -----------------------------------------------------------------------------
  c1_sensitivity_state <- map(1:21, ~sensitivity(txfx_c1, loo_states_c1, .x))
  c1_sensitivity_proc  <- map(1:21, ~sensitivity(txfx_c1, loo_proc_c1, .x))
  
  c2_sensitivity_state <- map(1:16, ~sensitivity(txfx_c2, loo_states_c2, .x))
  c2_sensitivity_proc  <- map(1:16, ~sensitivity(txfx_c2, loo_proc_c2, .x))
  
  num_states_c1 <- 21
  num_states_c2 <- 16
  
  loo_c1_mean_state <- Reduce(`+`, map(c1_sensitivity_state, ~.x$loostate))/length(c1_sensitivity_state)
  varhat_c1_state <- ((num_states_c1 - 1)/num_states_c1)*Reduce(`+`, map(c1_sensitivity_state, ~(.x$loostate - loo_c1_mean_state)^2))
  sehat_c1_state <- sqrt(varhat_c1_state + rep(y0varest, length(varhat_c1_state)))
  
  loo_c1_mean_proc <- Reduce(`+`, map(c1_sensitivity_proc, ~.x$loostate))/length(c1_sensitivity_proc)
  varhat_c1_proc <- ((num_states_c1 - 1)/num_states_c1)*Reduce(`+`, map(c1_sensitivity_proc, ~(.x$loostate - loo_c1_mean_proc)^2))
  sehat_c1_proc <- sqrt(varhat_c1_proc + rep(y0varest, length(varhat_c1_proc)))
  
  loo_c2_mean_state <- Reduce(`+`, map(c2_sensitivity_state, ~.x$loostate))/length(c2_sensitivity_state)
  varhat_c2_state <- ((num_states_c2 - 1)/num_states_c2)*Reduce(`+`, map(c2_sensitivity_state, ~(.x$loostate - loo_c2_mean_state)^2))
  sehat_c2_state <- sqrt(varhat_c2_state + rep(y0varest, length(varhat_c2_state)))
  
  loo_c2_mean_proc <- Reduce(`+`, map(c2_sensitivity_proc, ~.x$loostate))/length(c2_sensitivity_proc)
  varhat_c2_proc <- ((num_states_c2 - 1)/num_states_c2)*Reduce(`+`, map(c2_sensitivity_proc, ~(.x$loostate - loo_c2_mean_proc)^2))
  sehat_c2_proc <- sqrt(varhat_c2_proc + rep(y0varest, length(varhat_c2_proc)))
  
  list(txfx_c1 = txfx_c1, txfx_c2 = txfx_c2, loo_states_c1 = loo_states_c1, loo_states_c2 = loo_states_c2,
       loo_proc_c1 = loo_proc_c1, loo_proc_c2 = loo_proc_c2, sehat_c1_state = sehat_c1_state,
       sehat_c1_proc = sehat_c1_proc, sehat_c2_state = sehat_c2_state, sehat_c2_proc = sehat_c2_proc,
       c1_sensitivity_state = c1_sensitivity_state, c1_sensitivity_proc = c1_sensitivity_proc,
       c2_sensitivity_state = c2_sensitivity_state, c2_sensitivity_proc = c2_sensitivity_proc)
}

create_full_table <- function(txfx_results, subset) {
  se_state <- txfx_results[paste0("sehat_", subset, "_state")][[1]]
  se_proc  <- txfx_results[paste0("sehat_", subset, "_proc")][[1]]
  txfx <- txfx_results[paste0("txfx_", subset)][[1]]
  dof <- length(txfx_results[[paste0("loo_states_", subset)]]) - 1
  quant <- qt(0.975, dof)
  
  txfx %>%
    mutate(sehat = se_state, l95ci = psihat - quant*sehat, u95ci = psihat + quant*sehat,
           data = "Preferred", jackknife = "Data") %>%
    bind_rows(
      txfx %>%
        mutate(sehat = se_proc, l95ci = psihat - quant*sehat, u95ci = psihat + quant*sehat,
               data = "Preferred", jackknife = "All")
    ) 
}

point_est_table <- function(result_table) {
  result_table %>%
    mutate_at("sigma_estimate", ~gsub("_modeled", "", .)) %>%
    mutate_at("variable_subset", ~factor(., levels = covariate_group)) %>%
    mutate_at("variable_subset", ~as.numeric(.) - 1) %>%
    mutate_at("sigma_estimate", ~factor(., levels = sigma_estimator)) %>%
    mutate_at("weight_type", ~factor(., levels = estimator)) %>%
    filter(jackknife == "Data") %>%
    select(-sehat, -l95ci, -u95ci, -jackknife, -data) %>%
    spread(weight_type, psihat) %>%
    rename(`Variable subset` = variable_subset, Adjustment = sigma_estimate) %>%
    mutate_at("Adjustment", ~stringr::str_replace_all(., adjust_replace)) %>%
    mutate_at("Variable subset", as.character)
}

confint_table <- function(result_table) {
  result_table %>%
    mutate(ci = paste0("(", round(l95ci, 2), ", ", round(u95ci, 2), ")")) %>%
    select(-l95ci, -u95ci, -sehat) %>%
    spread(jackknife, ci) %>%
    mutate_at("sigma_estimate", ~gsub("_modeled", "", .)) %>%
    mutate_at("variable_subset", ~factor(., levels = covariate_group)) %>%
    mutate_at("variable_subset", ~as.numeric(.) - 1) %>%
    mutate_at("sigma_estimate", ~factor(., levels = sigma_estimator)) %>%
    mutate_at("weight_type", ~factor(., levels = estimator)) %>%
    select(-data) %>%
    ungroup() %>%
    arrange(weight_type, sigma_estimate, variable_subset) %>%
    select(`Weight type` = weight_type, Adjustment = sigma_estimate, 
           `Variable subset` = variable_subset,
           Psihat = psihat,
           `CI (states)` = Data, 
           `CI (proc)` = All,
    ) %>%
    mutate_at("Adjustment", ~stringr::str_replace_all(., adjust_replace))
}

paper_tables <- function(txfx_results) {
  c1_results_table <- create_full_table(txfx_results, "c1") 
  c2_results_table <- create_full_table(txfx_results, "c2") 
  
  c1_point_estimates <- point_est_table(c1_results_table)
  c2_point_estimates <- point_est_table(c2_results_table)
  
  c1_confint_table <- confint_table(c1_results_table)
  c2_confint_table <- confint_table(c2_results_table)
  
  c1_confint_filtered <- c1_confint_table %>%
    filter(`Variable subset` == 0) %>%
    select(-`Variable subset`)
  
  c2_confint_filtered <- c2_confint_table %>%
    filter(`Variable subset` == 0) %>%
    select(-`Variable subset`)
  
  list(c1_results_table = c1_results_table, c2_results_table = c2_results_table,
       c1_point_estimate = c1_point_estimates, c1_confint_table = c1_confint_table,
       c2_point_estimate = c2_point_estimates, c2_confint_table = c2_confint_table,
       c1_confint_filtered = c1_confint_filtered, c2_confint_filtered = c2_confint_filtered)
}

GenMainTable <- function(tables, sigma_exclude = "") {
  tab_new <- tables$c1_confint_filtered %>%
    filter(Adjustment != sigma_exclude) %>%
    group_by(`Weight type`) %>%
    mutate(Difference = Psihat - Psihat[Adjustment == "Unadjusted"]) %>%
    select(-`CI (states)`) %>%
    mutate(`Estimate (95% CI)` = paste0(round(Psihat, 2), " ", `CI (proc)`)) %>%
    select(`Weight type`, "Adjustment", `Estimate (95% CI)`, Difference)
  
  tab_new1 <- tables$c2_confint_filtered %>%
    filter(Adjustment != sigma_exclude) %>%
    group_by(`Weight type`) %>%
    mutate(Difference1 = Psihat - Psihat[Adjustment == "Unadjusted"]) %>%
    select(-`CI (states)`) %>%
    mutate(`Estimate (95% CI) \n (Early excluded)` = paste0(round(Psihat, 2), " ", `CI (proc)`)) %>%
    select(`Weight type`, "Adjustment", `Estimate (95% CI) \n (Early excluded)`, Difference1)
  
  tabnew2 <- left_join(tab_new, tab_new1) 
  return(tabnew2)
}

# Weight diagnostic plots -------------------------------------------------------------------------- 

# add weights to dataset
attach_weights <- function(data_subset, weights) {
  w1 = weights$`H-SBW`$sigma_uu_avg$weights
  w2 = weights$SBW$sigma_uu_avg$weights
  w3 = weights$`BC-SBW`$sigma_uu_avg$weights
  w4 = weights$`BC-HSBW`$sigma_uu_avg$weights
  
  data_subset %>%
    mutate(`H-SBW` = w1, SBW = w2, `BC-SBW` = w3,
           `BC-HSBW` = w4) 
}

# weights summed by state
aggregate_plot <- function(data_subset, weights, weight_type, file_extension) {
  myplot <- data_subset %>%
    attach_weights(weights) %>%
    select(state, matches("SBW")) %>%
    mutate_at(vars(matches("SBW")), 
              funs(pos = if_else(. > 0, ., 0),
                   neg = if_else(. <= 0, -1*., 0))) %>%
    select(state, contains("pos"), contains("neg")) %>%
    gather(key, value, -state) %>%
    mutate(pos = if_else(grepl("pos", key), "pos", "neg")) %>% 
    mutate_at("key", ~gsub("_pos|_neg", "", .)) %>% 
    mutate_at("key", ~factor(., levels = estimator)) %>% 
    filter(key %in% weight_type) %>%
    group_by(key, state, pos) %>% 
    dplyr::summarize(value = sum(value)) %>%
    mutate(value = if_else(pos == "pos", value, -value)) %>%
    mutate_at("state", factor) %>%
    ggplot(aes(x = state, y = value, fill = pos)) +
    guides(fill = FALSE) +
    geom_bar(stat = "identity") +
    facet_wrap(~key) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels("state"))) +
    xlab("") + ylab("Sum of Weights within State") 
    
  ggsave(paste0("../../02_Paper/01_Plots/weights-by-state-", file_extension ,".png"), myplot)
}

# Treatment effect plots ---------------------------------------------------------------------------

# plots main results -------------------------------------------------------------------------------
output_all_plots <- function(c1_data, c2_data, c1_results, c2_reuslts, txfx_results) {
  weight_type1 <- c("BC-HSBW", "H-SBW")
  weight_type2 <- c("BC-SBW", "SBW")
  weight_type4 <- c("SBW", "H-SBW")
  
  d1 <- c1_data$Preferred
  d2 <- c2_data$`Early expansion`
  w1 <- c1_results[[1]]
  w2 <- c2_results[[1]]
  txfx_c1 <- txfx_results$txfx_c1
  txfx_c2 <- txfx_results$txfx_c2
  
  aggregate_plot(d1[[1]], w1[[1]], weight_type1, "hsbw-c1") 
  aggregate_plot(d2[[1]], w2[[1]], weight_type1, "hsbw-c2") 
  aggregate_plot(d1[[1]], w1[[1]], weight_type2, "sbw-c1") 
  aggregate_plot(d1[[1]], w1[[1]], weight_type4, "sbw-hsbw-c1") 
}

# leave-one-out states analysis -----------------------------------------------------------------------------
centered <- function(loo_results, loo_state_names) {
  loo_results  %>%
    map2(loo_state_names, ~mutate(.x, state = .y)) %>%
    invoke(rbind, .) %>%
    group_by(variable_subset, weight_type, sigma_estimate) %>%
    mutate(mean = mean(loostate)) %>%
    select(-diff, -loostate, -state) %>%
    distinct() %>%
    mutate(diff = mean - psihat) %>%
    arrange(-abs(diff))
}

loo_table <- function(loo_results_proc, loo_results_state, loo_state_names) {
  loo_results_proc %>%
    map2(loo_state_names, ~mutate(.x, state = .y)) %>%
    invoke(rbind, .) %>%
    mutate(type = "proc") %>%
    rbind(
      loo_results_state %>%
        map2(loo_state_names, ~mutate(.x, state = .y)) %>%
        invoke(rbind, .) %>%
        mutate(type = "data")
    )
}

loo_heatmap <- function(result_table, sigma_estimator, type, file_extension) {
  replacements <- c("Homogeneous", "Heterogeneous", "Unadjusted")
  names(replacements) <- c("sigma_uu_avg", "sigma_uu_i_modeled", "sigma_zero")
  
  heatplot <- result_table %>%
    select(-psihat, -loostate) %>%
    filter(sigma_estimate %in% c("sigma_zero", sigma_estimator)) %>%
    filter(type == !!type) %>%
    filter(variable_subset == "None") %>%
    mutate_at("weight_type", ~factor(., estimator)) %>%
    mutate_at("sigma_estimate", ~stringr::str_replace_all(., replacements)) %>%
    ggplot(aes(x = weight_type, y = state, fill = diff)) +
    geom_tile() +
    facet_wrap(~sigma_estimate) +
    scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"),
                         name = "LOO estimate minus \n original estimate") +
    ylab("State Excluded") +
    xlab("Weight Type") 
  
  ggsave(paste0("../../02_Paper/01_Plots/loostate-sensitivity", file_extension, ".png"),
         heatplot)
}


produce_heatmaps <- function(txfx_results, file_extension) {
  c1_states <- names(txfx_results$loo_states_c1)
  c2_states <- names(txfx_results$loo_states_c2)
  
  c1_loo_table <- loo_table(txfx_results$loo_proc_c1, txfx_results$loo_states_c1, c1_states)
  c2_loo_table <- loo_table(txfx_results$loo_proc_c2, txfx_results$loo_states_c2, c2_states)
  
  main <- c("c1-state-main", "c1-proc-main", "c2-state-main", "c2-proc-main")
  main <- paste0(main, file_extension)
  full <- c("c1-state-full", "c1-proc-full", "c2-state-full", "c2-proc-full")
  full <- paste0(full, file_extension)
  othr <- c("c1-state-uu-avg", "c1-proc-uu-avg", "c2-state-uu-avg", "c2-proc-uu-avg")
  othr <- paste0(othr, file_extension)

  c1_diff_proc <- txfx_results$c1_sensitivity_proc %>%
    map2(names(txfx_results$loo_states_c1), ~mutate(.x, state = .y)) %>%
    invoke(rbind, .)
  
  c2_diff_proc <- txfx_results$c2_sensitivity_proc %>%
    map2(names(txfx_results$loo_states_c2), ~mutate(.x, state = .y)) %>%
    invoke(rbind, .)
  
  loo_heatmap(c1_diff_proc, c("sigma_uu_avg"), "data", othr[2])
  loo_heatmap(c2_diff_proc, c("sigma_uu_avg"),"data", othr[4])
}

centering_check <- function(txfx_results) {
  c1_sensitivity_state <- txfx_results$c1_sensitivity_state
  c2_sensitivity_state <- txfx_results$c2_sensitivity_state
  c1_sensitivity_proc  <- txfx_results$c1_sensitivity_proc
  c2_sensitivity_proc  <- txfx_results$c2_sensitivity_proc
  c1_names <- names(all_txfx$loo_states_c1)
  c2_names <- names(all_txfx$loo_states_c2)
  
  t1 <- centered(c1_sensitivity_state, c1_names)
  t2 <- centered(c1_sensitivity_proc, c1_names) 
  t3 <- centered(c2_sensitivity_state, c2_names) 
  t4 <- centered(c2_sensitivity_proc, c2_names) 
  
  list(c1_state = t1, c1_proc = t2, 
       c2_state = t3, c2_proc = t4)
}

