# program: 06_jackknife.R
# purpose: calculate point estimates when recomputing entire imputation procedure for ETC, OATE
# author: max rubinstein
# date modified: december 22, 2020

source('03_Analysis/01_calibrate-data.R')
source('03_Analysis/02_model-estimation.R')

# leave one cluster out jackknife functions -------------------------------------------------------------------
generate_covariance_matrix <- function(Sigma_UU_i, original_data, variables, count_xwalk) {
  Sigma_SS_i <- map(1:nrow(original_data), ~generate_count_covariance(original_data, variables,
                                                                      count_xwalk, .x))
  
  Sigma_UU_modeled <- Reduce(`+`, map2(Sigma_SS_i, Sigma_UU_i, ~.x*.y))/length(Sigma_SS_i)
  
  Sigma_UU_i_modeled <- map(Sigma_SS_i, ~Sigma_UU_modeled/.x)
  
  Sigma_UU_Avg <- rep(list(Reduce(`+`, Sigma_UU_i)/length(Sigma_UU_i)), nrow(original_data))
  
  all_vars <- sort(variables)
  
  Sigma_Zero <- rep(list(matrix(rep(0, (length(all_vars))^2), length(all_vars), length(all_vars))),
                    nrow(original_data))
  
  list("sigma_uu_i_modeled" = Sigma_UU_i_modeled, 
       "sigma_uu_avg" = Sigma_UU_Avg, 
       "sigma_zero" = Sigma_Zero)
}

# jackknife for etc
cluster_jackknife_etc <- function(treatment_data, sigma_uu_i, variables, targets, count_xwalk) {

  data <- treatment_data %>%
    nest(-state) 
  
  state_names <- data$state
  all_data <- list()
  
  for (i in 1:nrow(data)) {
    jackknife <- data[-i,] %>%
      unnest()
    
    variables_new <- variables
    
    Sigma_X_hat <- cov(as.matrix(jackknife[,variables]))
    rank_X1 <- Matrix::rankMatrix(Sigma_X_hat)
    
    if (rank_X1 != ncol(Sigma_X_hat)) {
      qr.x <- qr(Sigma_X_hat)
      Sigma_X_hat <- Sigma_X_hat[, qr.x$pivot[1:rank_X1]]
      targets <- targets[names(targets) %in% colnames(Sigma_X_hat)]
      variables_new <- variables[variables %in% names(targets)]
    }
    
    full_ids <- paste0(treatment_data$state, "_", treatment_data$cpuma)
    sub_ids <- paste0(jackknife$state, "_", jackknife$cpuma)
    observation_indices <- unlist(map(full_ids, ~any(grepl(.x, sub_ids))))
    observation_indices <- grep(TRUE, observation_indices)
    
    variable_indices <- variables %in% variables_new

    Sigma_UU_i_f <- sigma_uu_i[observation_indices]
    Sigma_UU_i_f <- map(Sigma_UU_i_f, ~.x[variable_indices, variable_indices])
    
    sigma_uu <- generate_covariance_matrix(Sigma_UU_i_f, jackknife, variables_new, count_xwalk)
    
    kappa_list <- calculate_kappa_all(jackknife, sigma_uu, variables_new) 
    
    data_list <- transform_data(jackknife, variables_new, kappa_list) %>%
      mutate(data = map(data, ~arrange(.x, state, cpuma))) 

    all_data <- append(all_data, list(data_list))
  }
  names(all_data) <- state_names
  all_data
}

etc_weights <- function(jackknife_data, tol_list, cov_groups, targets, max_iter = 1e6) {
  data <- map(jackknife_data, ~.x$data)
  weights <- map(data, ~iterate_covariate_subsets(.x, tol_list, targets,
                                                      linf_imbalance, stop_criterion = 0.5, 
                                                      max_iter = max_iter))
  weights <- map(weights, ~set_names(.x, cov_groups))
  weights <- map(1:length(jackknife_data), ~map(weights[[.x]], ~set_names(.x, model_names)))
  weights <- map(1:length(jackknife_data), ~map(weights[[.x]], set_cov_mod_name))
  names(weights) <- names(jackknife_data)
  
  list(weights = weights, data = jackknife_data)
}

# jackknife for oate
cluster_jackknife_oate <- function(full_data, variables_list, cdat_imputed, tdat_imputed,
                                   sigma_uu_i_tx, sigma_uu_i_ct, count_xwalk) {
  
  nest_data <- full_data %>%
    arrange(state, cpuma) %>%
    nest(-state, -treatment) 
  
  state_list <- nest_data$state

  variables <- variables_list[[1]]
  
  remove_collinear <- function(data, variables, targets, treatment) {
    variables_new <- variables
    X <- as.matrix(data[variables])
    
    rank_X <- Matrix::rankMatrix(X)

    if (rank_X != ncol(X)) {
      qr.x <- qr(X)
      X <- X[, qr.x$pivot[1:rank_X]]
      targets <- targets[names(targets) %in% colnames(X)]
      variables_new <- variables[variables %in% names(targets)]
    }
    
    for (variable in variables_new) {
      if (var(X[,variable]) == 0) {
        variables_new <- variables_new[-grep(variable, variables_new)]
      }
    }
    
    list(variables = variables_new)
  }
  
  all_data <- list()

  for (i in 1:nrow(nest_data)) {
    state <- nest_data$state[i]
    treatment <- nest_data$treatment[i]
    jackknife <- nest_data[-i,] 
    variables_tx <- variables
    
    X_new <- filter(jackknife, treatment == !!treatment) %>%
      unnest() %>%
      arrange(state, cpuma)
    
    X_full <- filter(nest_data, treatment == !!treatment) %>%
      unnest() %>%
      arrange(state, cpuma)
    
    full_ids <- paste0(X_full$state, "_", X_full$cpuma)
    sub_ids <- paste0(X_new$state, "_", X_new$cpuma)
    observation_indices <- unlist(map(full_ids, ~any(grepl(.x, sub_ids))))
    observation_indices <- grep(TRUE, observation_indices)
    
    X_old <- switch(as.character(treatment),
                    "1" = cdat_imputed,
                    "0" = tdat_imputed)
    
    sigma_uu_i <- switch(as.character(treatment),
                    "1" = sigma_uu_i_tx,
                    "0" = sigma_uu_i_ct)
    
    
    if (treatment == 0) {
      index <- grep("repub_total_control", variables_tx)
      variables_tx <- variables_tx[-index]
      sigma_uu_i <- map(sigma_uu_i, ~.x[-index, -index])
    }
    
    aligned <- remove_collinear(X_new, variables_tx, targets, treatment)
    
    variables_new <- aligned$variables

    variable_indices <- variables_tx %in% variables_new

    sigma_uu_i <- sigma_uu_i[observation_indices]
    sigma_uu_i <- map(sigma_uu_i, ~.x[variable_indices, variable_indices])
    
    sigma_uu <- generate_covariance_matrix(sigma_uu_i, X_new, variables_new, count_xwalk)

    kappa_list <- calculate_kappa_all(X_new, sigma_uu, variables_new)

    trans_data <- transform_data(X_new, variables_new, kappa_list) %>%
      mutate(data = map(data, ~arrange(.x, state, cpuma)))

    final_data <- trans_data %>%
      mutate(treatment = treatment) %>%
      bind_rows(X_old) %>%
      unnest() %>% 
      nest(-key) %>%
      mutate(data = map(data, ~arrange(.x, state, cpuma)))
    
    all_data <- append(all_data, list(final_data))
  }
  
  data_list <- all_data %>% 
    map(~.x$data) %>%
    map(~set_names(.x, cov_models))
}

oate_weights <- function(jackknife_data, variables_list, state_list) {
  state_list <- names(jackknife_data)
  
  weights <- map(jackknife_data, ~iterate_covariates_oate(.x, variables_list))
  
  final <- list(weights = weights, data = jackknife_data) %>%
    map(~set_names(.x, state_list))
  final
}

set_cov_mod_name <- function(result_list) {
  map(result_list, ~set_names(.x, cov_models))
}

