# program: 01_error-variance-estimation.R
# purpose: estimate the CPUMA-level estimation error and calibrate covariates/outcomes
# author: max rubinstein
# date modified: december 14, 2020

# estimate covariance matrix for each row of data using replicates
generate_full_cov_matrix <-function(original, replicates, variables, row_num) {
  selections <- sort(variables)
  
  data_point <- map(replicates, ~.x %>% 
                     ungroup() %>% 
                     arrange(state, cpuma) %>%
                     select(all_of(selections)) %>% 
                     slice(row_num) %>% 
                     as.matrix())
  
  original_point <- original %>% 
    ungroup() %>%
    arrange(state, cpuma) %>%
    select(all_of(selections)) %>% 
    slice(row_num) %>% 
    as.matrix()
    
  cov_matrix <- 4*Reduce(`+`, map(data_point, ~t(.x - original_point) %*% (.x - original_point)))/length(data_point)
  
  cov_matrix <- cov_matrix %>%
    as.data.frame() %>%
    rownames_to_column() %>% 
    gather(key, value, -rowname) %>%
    mutate(year1 = stringr::str_extract_all(rowname, '201[1-4]', simplify = TRUE),
           year2 = stringr::str_extract_all(key, '201[1-4]', simplify = TRUE)) %>% 
    mutate(value = if_else(year1 != year2 & !(year1 == '' | year2 == ''), 0, value),
           value = if_else((year1 == '2014' | year2 == '2014') & (year1 != '2014' | year2 != '2014'), 0, value)) %>%
    select(-year1, -year2) %>%
    spread(key, value) %>%
    select(-rowname) %>%
    as.matrix()
  
  assertthat::assert_that(identical(colnames(cov_matrix), selections))
  assertthat::assert_that(length(grep(FALSE, diag(cov_matrix) >= 0)) == 0)
  
  cov_matrix
}

# generate counts covariance_matrix
generate_count_covariance <- function(original_data, variables, count_xwalk, row_num) {
  counts <- original_data %>% 
    ungroup() %>% 
    arrange(state, cpuma) %>%
    select(contains('count')) %>% 
    slice(row_num) %>%
    gather(sample_var, counts)
  
  row_counts = count_xwalk %>%
    left_join(counts, by = 'sample_var') %>%
    replace_na(list(counts = 1e7)) %>%
    mutate(root_counts = sqrt(counts)) %>%
    filter(Variable %in% variables) %>%
    arrange(Variable) 

  row_count_covs <- row_counts$root_counts %*% t(row_counts$root_counts)
  colnames(row_count_covs) <- sort(row_counts$Variable)
  row_count_covs
}

# generate different estimators of measurement error covariance matrices
generate_all_covariance_estimates <- function(original_data, Sigma_UU_i, count_xwalk, variables) {
  original_data <- arrange(ungroup(original_data), treatment, state, cpuma)
  
  # 1. Rescaled version of one covariance matrix
  Sigma_SS_i <- map(1:nrow(original_data), ~generate_count_covariance(original_data, variables,
                                                                      count_xwalk, .x))
  
  Sigma_UU_modeled <- Reduce(`+`, map2(Sigma_SS_i, Sigma_UU_i, ~.x*.y))/length(Sigma_SS_i)
  
  # verify that the order of the column names in the same order
  assert_that(identical(colnames(Sigma_UU_modeled), colnames(Sigma_SS_i[[1]])))
  
  Sigma_UU_i_modeled <- map(Sigma_SS_i, ~Sigma_UU_modeled/.x)
  
  # 2. Averaged full covariance matrix
  Sigma_UU_Avg <- rep(list(Reduce(`+`, Sigma_UU_i)/length(Sigma_UU_i)), nrow(original_data))
  
  # 3. No measurement error
  all_vars <- sort(variables)
  
  Sigma_Zero <- rep(list(matrix(rep(0, (length(all_vars))^2), length(all_vars), length(all_vars))),
                   nrow(original_data))
  
  list("sigma_uu_i_modeled" = Sigma_UU_i_modeled, 
       "sigma_uu_avg" = Sigma_UU_Avg, 
       "sigma_zero" = Sigma_Zero)
}

# calculate kappa ------------------------------------------------------------------------------------------
calculate_kappa_all <- function(original_data, Sigma_UU_estimates, variables) {
  
  calculate_avg_cov <- function(cov_matrix_list) {
    Reduce(`+`, cov_matrix_list)/length(cov_matrix_list)
  }
  
  calculate_kappa_list <- function(Sigma_X_hat, Sigma_U_hat_list) {
    map(Sigma_U_hat_list, ~solve(Sigma_X_hat + .x) %*% Sigma_X_hat)
  }
  
  set_colnames <- function(mymatrix, mynames) {
    colnames(mymatrix) <- mynames
    mymatrix
  }
  
  # calculate covariance matrix on original dataset
  selections <- sort(variables)
  
  Sigma_W <- cov(original_data[selections])

  # test that order of covariance matrices identical
  test_case <- Sigma_UU_estimates[[1]][[1]]
  assert_that(identical(colnames(Sigma_W), colnames(test_case)))
  
  # calculate list of average covariance matrices
  Sigma_U_avg <- map(Sigma_UU_estimates, ~calculate_avg_cov(.x))

  # calculate list of sigma_x estimates based on averages
  Sigma_X_hat <- map(Sigma_U_avg, ~Sigma_W - .x)

  # generate unit-level kappa values
  kappa_all  <- map2(Sigma_X_hat, Sigma_UU_estimates, ~calculate_kappa_list(.x, .y))

  names(kappa_all) <- c("sigma_uu_i_modeled", "sigma_uu_avg", "sigma_zero")
  
  kappa_all
}

# impute covariates ----------------------------------------------------------------------------------------
transform_data <- function(original_data, variables, kappa_list) {
  assert_that(is.list(kappa_list[[1]]), msg = 'kappa_list must be a list of lists')
  assert_that(!is.null(names(kappa_list)), msg = 'kappa_list must be named')
  
  original_data <- arrange(ungroup(original_data), treatment, state, cpuma)
  
  scale_data <- function(data, variables) {
    scaled_data = data %>%
      ungroup() %>%
      select(variables) %>%
      scale(center = TRUE, scale = FALSE)
    
    scaled_data
  }
  
  all_variables <- sort(variables)
  
  data_process <- function(original_data, kappa, treat_level) {
    indices <- grep(treat_level, original_data$treatment)
    kappa_list <- map(kappa, ~.x[indices])
    data <- subset(original_data, treatment == treat_level)
    
    scaled_data <- scale_data(data, all_variables)
    
    center <- as.matrix(attr(scaled_data, 'scaled:center'))
    
    adjusted_data <- scaled_data %>%
      as_tibble() %>%
      mutate(id = 1:nrow(.)) %>%
      nest(-id) %>% 
      mutate_at('data', ~map(., ~t(as.matrix(.x))))
    
    for(i in 1:length(kappa_list)) {
      adjusted_data <- adjusted_data %>%
        mutate(!!(names(kappa_list)[i]) := kappa_list[[i]])
    }
    
    other_variables <- names(data)[!names(data) %in% all_variables]
    
    other_variables <- data %>%
      select(cpuma, state, all_of(other_variables))
    
    transformed_data <- adjusted_data %>%
      gather(key, kappa, -id, -data) %>%
      mutate(transformed_X = map2(data, kappa, ~as_tibble(t(center) + t(.x) %*% .y))) %>%
      mutate(transformed_X = map2(transformed_X, kappa, ~set_names(.x, colnames(.y)))) %>%
      select(-data, -kappa, -id) %>%
      unnest(cols = c(transformed_X)) %>%
      nest(-key) 
    
    transformed_data$data <- map(transformed_data$data, 
                                 ~mutate(.x, cpuma = data$cpuma, state = data$state) %>%
                                   left_join(other_variables, by = c('cpuma', 'state')) %>%
                                   mutate_if(is.numeric, ~if_else(abs(.) < 1e-10, 0, .)) %>%
                                   arrange(state, cpuma))
    
    transformed_data
  }
  
  map(unique(original_data$treatment), ~data_process(original_data, kappa_list, .x)) %>%
    invoke(rbind, .) 
}

