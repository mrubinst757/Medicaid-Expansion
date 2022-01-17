# program: 02_model-estimation.R
# purpose: pestimate weighting models
# author: max rubinstein
# date modified: december 14, 2020

# load libraries and read data ----------------------------------------------------------------
library(tidyverse)
library(assertthat)
library(optweight)

# calculate SBW/HSBW weights
generate_sbw_weights <- function(X_1, tols, targets, sigma2.y = 1, re = 0,
                                 state_num = NULL, max_iter) {
  variable_names <- names(tols)
  if(is.null(state_num)) state_num = rep(1, nrow(X_1)) 
  
  balancing_formula <- as.formula(paste0(" ~ ", paste0(variable_names, collapse = "+")))

  calc_sbw_weights <- partial(optweight.svy, formula = balancing_formula, std.cont = FALSE, 
                              std.binary = FALSE,
                             targets = targets, sigma2.y = sigma2.y, re = re, 
                             group_n = state_num, max_iter = max_iter)
  
  results.sbw <- calc_sbw_weights(data = X_1, tols = tols)
  
  iter = 0
  
  # adjust error tolerance if the program does not converge
  while(results.sbw$info$status != "solved" & iter < 10) {
    if(tols[grep("child|married", variables)][1] < 5 & results.sbw$info$status != "solved") {
      tols[grep("child|married", variables)] = tols[grepl("child|married", variables)] + 0.5
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("educ|student|inc|urban", variables)][1] < 5 & results.sbw$info$status != "solved") {
      tols[grep("educ|student|inc|urban", variables)] = tols[grepl("educ|student|inc|urban", variables)] + 0.5
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("white", variables)][1] < 5 & results.sbw$info$status != "solved") {
      tols[grep("black", variables)] = tols[grepl("black", variables)] + 0.25
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("disability|foreign|hispanic|citiz", variables)][1] < 5 & results.sbw$info$status != "solved") {
      tols[grep("disability|foreign|hispanic|citiz", variables)] = tols[grepl("disability|foreign|hispanic|citiz", variables)] + 1
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("growth|hh", variables)][1] < 1 & results.sbw$info$status != "solved") {
      tols[grep("growth|hh", variables)] = tols[grepl("growth", variables)] + 0.1
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("unemp", variables)][1] <= 0.25 & results.sbw$info$status != "solved") {
      tols[grep("unemp", variables)] = tols[grepl("unemp", variables)] + 0.025
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("unins", variables)][1] <= 0.25 & results.sbw$info$status != "solved") {
      tols[grep("unins", variables)] = tols[grepl("unins", variables)] + 0.025
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    if(tols[grep("repub", variables)][1] <= 50 & results.sbw$info$status != "solved") {
      tols[grep("repub", variables)] = tols[grepl("repub", variables)] + 5
      results.sbw = calc_sbw_weights(data = X_1, tols = tols)
    }
    iter = iter + 1
  }
  list(weights = results.sbw$weights, tols = results.sbw$tols, status = results.sbw$info$status)
}

# calculate BC weights
debias_sbw_weights <- function(sbw_weights, X_1, tols, targets, distance, autotune = TRUE, 
                               lambda = 0, stop_criterion = NULL, re = 0,
                               state_num = NULL, max_iter) {
  # modify variables for data subset
  variable_names <- names(tols)
  variables_to_remove <- names(tols)[grepl(100, tols)]
  assert_that(identical(names(targets), names(tols)))
  targets <- targets[!grepl(100, tols)]
  X_1 <- X_1 %>%
    select(all_of(variable_names)) %>%
    select(-variables_to_remove)
  
  assert_that(ncol(X_1) == length(targets))
  
  # calculate imbalance from weighting
  weights <- sbw_weights$weights/sum(sbw_weights$weights)
  X_1 <- as.matrix(X_1)
  orig_imbalance <- t(targets - t(X_1) %*% weights)
  imb_measure <- distance(orig_imbalance)
  
  if (re == 0) {
    sigma <- diag(rep(1, nrow(X_1)))
  }
  
  if (re > 0) {
    assert_that(!is.null(state_num))
    sigma <- map(state_num, ~matrix(rep(re, .x^2), .x, .x)) %>%
      Matrix::bdiag()
    diag(sigma) <- 1
  }
  
  sigma_inv <- solve(sigma)
  
  if (autotune == TRUE) {
    lambda <- 1e7
    
    if (imb_measure < stop_criterion) {
      w_new <- weights
    }
    count <- 0

    imbalance <- t(targets - t(X_1) %*% weights)
    
    while((imb_measure > stop_criterion) & count <= 50) {
      w_new <- adjust_weights(weights, lambda, X_1, targets, sigma_inv)
      imbalance <- t(targets - t(X_1) %*% w_new)
      imb_measure <- distance(imbalance)
      lambda <- lambda/10
      count <- count + 1
      iter <- count %% 10
      if (iter == 0 & imb_measure > stop_criterion) {
        imb_measure <- imb_measure + 0.1
        lambda <- 1e7
      }
    }
  }
  if (autotune == FALSE)  {
    w_new <- adjust_weights(weights, lambda, X_1, targets, sigma_inv)
    imbalance <- t(targets - t(X_1) %*% w_new)
    imb_measure = distance(imbalance)
    lambda <- lambda/10
  }
  list(weights = w_new*length(weights), imbalance = imbalance, lambda = lambda*10)
}

collinearity_check <- function(X_1, targets) {
  constants <- unlist(map(1:ncol(X_1), ~var(X_1[, .x])))
  constant_index <- unlist(map(constants, ~near(.x, 0))) %>%
    grep(TRUE, .)
  
  if (length(constant_index) > 0) {
    X_1 <- X_1[, -constant_index]
    targets <- targets[-constant_index]
  }
  
  X_1 <- cbind(rep(1, nrow(X_1)), X_1)
  targets = c(1, targets)
  rank_X1 <- Matrix::rankMatrix(X_1)
  
  if (rank_X1 != ncol(X_1)) {
    qr.x <- qr(X_1)
    X_1 <- X_1[, qr.x$pivot[1:rank_X1]]
    targets <- targets[names(targets) %in% colnames(X_1)]
  }
  
  list(X_1 = X_1, targets = targets)
}

adjust_weights <- function(weights, lambda, X_1, targets, sigma_inv) {
  result <- collinearity_check(X_1, targets)
  X_1 <- result$X_1; targets <- result$targets
  I_d <- lambda*diag(rep(1, ncol(X_1)))
  I_d[1, 1] <- 0 #don"t want to regularize the summing to one
  XX <- t(X_1) %*% sigma_inv %*% X_1
  XX_inv <- solve(XX + I_d)
  adjustment1 <- t(targets) %*% XX_inv %*% t(X_1) %*% sigma_inv
  adjustment2 <- t(t(X_1) %*% weights) %*% XX_inv %*% t(X_1) %*% sigma_inv
  adjustment <- adjustment1 - adjustment2
  w_new <- as.vector(weights + adjustment)
  w_new
}

# calculate list of weights
generate_weight_list <- function(tols, data_list, targets, distance,
                                 stop_criterion = NULL, max_iter) {
  state_num <- as.numeric(table(data_list[[1]]$state))
  sbw <- list(); hsbw <- list(); hsbw1 <- list()
  tols_new <- rep(list(tols), length(data_list))
  for (i in 1:length(data_list)) {
    sbw[[i]]  <- generate_sbw_weights(data_list[[i]], tols_new[[i]], targets, sigma2.y = 1, re = 0, max_iter = max_iter)
    if (any(sbw[[i]]$tols != tols_new[[i]])) {
      tols_new[[i]] <- sbw[[i]]$tols
    }
    hsbw[[i]] <- generate_sbw_weights(data_list[[i]], tols_new[[i]], targets, sigma2.y = 1, re = 0.2, 
                                                 state_num = state_num, max_iter = max_iter)
  }
  
  no_weights <- rep(list(list(weights = rep(1, length(sbw[[1]]$weights)))), length(sbw))
  
  debiased_sbw <- map2(sbw, data_list, ~debias_sbw_weights(.x, .y, tols, 
                                                           targets, distance, stop_criterion = stop_criterion))
  
  debiased_hsbw <- map2(hsbw, data_list, ~debias_sbw_weights(.x, .y, tols, 
                                                             targets, distance, stop_criterion = stop_criterion, 
                                                             re = 1/6, state_num = state_num))

  list(SBW = sbw, HSBW = hsbw, `BC-SBW` = debiased_sbw, `BC-HSBW` = debiased_hsbw)
}

# iterate weights across a list
iterate_covariate_subsets <- function(data_list, tol_list, targets, distance, 
                                      stop_criterion = NULL, max_iter = 1e6) {
  map(tol_list, ~generate_weight_list(.x, data_list, targets, distance,
                                      stop_criterion, max_iter))
}

# function that inputs list of datasets and filters states to subset
data_subset <- function(data_list, state_name) {
  map(data_list, ~filter(.x, !state %in% state_name))
}

# read variable names by variable group
read_varnames <- function(group = "") {
  variable_names <- read_csv("../02_Specs/tol_specs.csv") %>%
    filter(`Reference Variable` == 0) %>%
    filter(!Group == group) %>%
    arrange(Variable) %>%
    .$Variable
  variable_names
}

# distance measure
linf_imbalance <- function(imbalance) {
  max(abs(imbalance))
}

remove_states <- function(data_list) {
  res <- list()
  for(i in 1:length(data_list)) {
    new_list <- map(unique(data_list[[i]]$state), ~filter(data_list[[i]], state != .x))
    res[[i]] <- append(list(data_list[[i]]), new_list)
    names(res[[i]]) <- c("Preferred", unique(data_list[[i]]$state))
  }
  return(transpose(res))
}

