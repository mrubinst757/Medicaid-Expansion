library(tidyverse)
library(optweight)

GenerateHierPopulation <- function(number_states, SS_cor, SNR.xw, 
                                   rho, X_cor, var.y = 2, var.x = 2, hetero = TRUE) {
  
  #Note: Variance of unit level X = 2; Var(Y | X) = 2
  
  # specify variance of error term e in model: W = X + e
  var.v <- (1 - SNR.xw) * 2 / SNR.xw
  
  # specify variance of y
  var.re = var.y * rho
  var.yr = var.y - var.v - var.re
  
  # specify state-level X variance
  var.s <- var.x * SS_cor
  
  # specify correlation between variables (half at state level, half county level)
  cov.s <- var.s * X_cor
  
  # specify implied additional within state variability
  var.q <- var.x - var.s
  cov.q <- var.q * X_cor 
  
  # specify correlation matrix of state-level X values
  cor.mat.s <- matrix(c(var.s, cov.s, cov.s, 
                        cov.s, var.s, cov.s, 
                        cov.s, cov.s, var.s), 3, 3)
  
  cor.mat.q <- matrix(c(var.q, cov.q, cov.q,
                        cov.q, var.q, cov.q,
                        cov.q, cov.q, var.q), 3, 3)
  
  cor.mat.x <- cor.mat.s + cor.mat.q
  
  # draw state-level X from MVN distribution
  if (var.s != 0) {
    S <- MASS::mvrnorm(number_states, c(0, 0, 0), cor.mat.s)
  }
  if (var.s == 0) {
    S <- matrix(rep(0, 3*number_states), number_states, 3)
  }
  
  p = floor(rexp(number_states, 0.1)) + 10
  data <- tibble(
    id = 1:number_states,
    S1 = S[,1],
    S2 = S[,2],
    S3 = S[,3],
    p = p,
    re = rnorm(number_states, 0, sqrt(var.re))
  ) 
  
  X <- data[,c("S1", "S2", "S3")]
  X <- map2(1:nrow(data), data$p, ~MASS::mvrnorm(.y, mu = as.numeric(X[.x,]), Sigma = cor.mat.q)) %>%
    invoke(rbind, .) %>%
    as_tibble() %>%
    set_names(c("X1", "X2", "X3"))
  
  data <- as.data.frame(lapply(data, rep, data$p)) 
  data <- bind_cols(data, X)
  
  nxvector <- runif(nrow(data), 300, 2300)
  
  if (hetero == TRUE) {
    m1x <- 1/(log(2300)/2000 - log(300)/2000)
    SigmaV <- var.v / (log(2300)/2000 - log(300)/2000)
    nxvector.true <- nxvector
  }
  
  if (hetero == FALSE) {
    m1x <- 1
    SigmaV <- var.v
    nxvector.true <- rep(1, nrow(data))
  }
  
  data <- data %>%
    mutate(nxtrue = nxvector.true,
           nxobs = nxvector,
           SigmaVV = SigmaV / nxtrue) %>%
    mutate(err.y = rnorm(nrow(.), 0, sqrt(var.yr))) %>%
    mutate(Y0 = X1 + X2 + X3,
           Y = Y0 + re + err.y,
           W1 = rnorm(nrow(.), X1, sqrt(SigmaVV)),
           W2 = rnorm(nrow(.), X2, sqrt(SigmaVV)),
           W3 = rnorm(nrow(.), X3, sqrt(SigmaVV)),
           V1 = SigmaVV,
           V2 = SigmaVV,
           V3 = SigmaVV,
           J  = rnorm(nrow(.), Y, sqrt(SigmaVV)))
  
  SigmaVV <- SigmaV / m1x
  cor.mat.w <- cor.mat.x + diag(rep(SigmaVV, 3))
  
  return(list(data = data, 
              cor.mat.s = cor.mat.s,
              cor.mat.q = cor.mat.q,
              cor.mat.x = cor.mat.x, 
              cor.mat.w = cor.mat.w, 
              rho = rho,
              var.y = var.y, 
              SigmaV = SigmaV))
}

CalcSigSS <- function(data) {
  blocks <- as.numeric(table(data$id))
  finish <- cumsum(blocks)
  start <- lag(finish) + 1; start[1] <- 1
  indices <- map2(start, finish, ~c(.x:.y)) 
  SigmaSS_list <- list()
  
  for (i in 1:length(indices)) {
    index <- gtools::combinations(n = length(indices[[i]]), r = 2, v = indices[[i]], repeats.allowed = FALSE)
    SigmaSS_list[[i]] <- map2(index[,1], index[,2], ~as.numeric(data[.x, c("W1", "W2", "W3")]) %*% 
                                t(as.numeric(data[.y, c("W1", "W2", "W3")])))
  }
  
  SigmaSS_list <- map(SigmaSS_list, ~Reduce(`+`, .x)/length(.x))
  
  return(SigmaSS_list)
}

CalcXhat <- function(data, estimate = "Truth", cor.mat.x = NULL, cor.mat.w = NULL,
                     SigmaSS = NULL, Sigma_vv_est) {
  
  if (estimate == "Homogeneous") {
    
    W0hat <- c(mean(data$W1), mean(data$W2), mean(data$W3))
    Wmat <- as.matrix(data[, c("W1", "W2", "W3")])
    Sigma_WW_est <- cov(Wmat)
    Sigma_vv_avg <- Reduce(`+`, Sigma_vv_est)/length(Sigma_vv_est)
    Sigma_XX_est <- Sigma_WW_est - Sigma_vv_avg
    kappa_hat <- solve(Sigma_WW_est) %*% Sigma_XX_est
    Xhat <- W0hat + (Wmat - W0hat) %*% t(kappa_hat)
    data <- mutate(data, Xhat.hom1 = Xhat[,1], Xhat.hom2 = Xhat[,2], Xhat.hom3 = Xhat[,3])
  }
  if (estimate == "Heterogeneous") {
    
    W0hat <- c(mean(data$W1), mean(data$W2), mean(data$W3))
    Wmat <- as.matrix(data[, c("W1", "W2", "W3")])
    Sigma_WW_est <- cov(Wmat)
    Sigma_vv_avg <- Reduce(`+`, Sigma_vv_est) / length(Sigma_vv_est)
    Sigma_XX_est <- Sigma_WW_est - Sigma_vv_avg
    
    nx_matrices <- map(1:nrow(data), ~diag(rep(data$nxobs[.x], 3)))
    Sigma_V_est <- map2(Sigma_vv_est, data$nxobs, ~.x*.y) 
    Sigma_V_est <- Reduce(`+`, Sigma_V_est) / length(Sigma_vv_est)
    
    Sigma_vv_i_est <- map(nx_matrices, ~diag(diag(Sigma_V_est)/diag(.x)))
    Sigma_WW_list  <- map(Sigma_vv_i_est, ~Sigma_XX_est + .x)
    kappa_hat_list <- map(Sigma_WW_list, ~solve(.x) %*% Sigma_XX_est)
    Xhat_list <- map2(kappa_hat_list, 1:nrow(data), ~t(W0hat + t(.x) %*% (Wmat[.y,] - W0hat)))
    Xhat <- invoke(rbind, Xhat_list)
    data <- mutate(data, Xhat.het1 = Xhat[,1], Xhat.het2 = Xhat[,2], Xhat.het3 = Xhat[,3])
  }
  if (estimate == "Correlated") {
    
    blocks <- as.numeric(table(data$id))
    finish <- cumsum(blocks)
    start <- lag(finish) + 1; start[1] <- 1
    indices <- map2(start, finish, ~c(.x:.y)) 
    
    W0hat <- c(mean(data$W1), mean(data$W2), mean(data$W3))
    Wmat <- as.matrix(data[, c("W1", "W2", "W3")])
    Sigma_WW <- cov(Wmat)
    Sigma_vv_avg <- Reduce(`+`, Sigma_vv_est) / length(Sigma_vv_est)
    Sigma_XX <- Sigma_WW - Sigma_vv_avg
    
    CreateMatrix <- function(block, Sigma) {
      res <- list()
      for (i in 1:length(block)) {
        res[[i]] <- assign_in(block, i, Sigma) 
        if (i > 1) {
          res[[i]] <- modify_at(res[[i]], 1:(i-1), t) 
        }
        res[[i]] <- invoke(rbind, res[[i]])
      }
      res <- invoke(cbind, res)
      return(res)
    }
    
    ImputeBlock <- function(Wmat, block, indices, SigmaXX, SigmaWW) {
      Wmat0 <- scale(Wmat, center = TRUE, scale = FALSE)
      nb <- length(block)
      Wrows <- map(indices, ~Wmat0[.x,]) %>%
        map(as.numeric) %>%
        unlist()
      SX <- CreateMatrix(block, SigmaXX)
      SW <- CreateMatrix(block, SigmaWW)
      kappa <- solve(SW) %*% SX
      wimp <- rep(colMeans(Wmat), length(block)) + t(kappa) %*% Wrows
      Xhat <- map(seq(1, ncol(Wmat0)*nb, ncol(Wmat0)), ~seq(.x, (.x + ncol(Wmat0) - 1), 1)) %>%
        map(~wimp[.x]) %>%
        invoke(rbind, .)
      return(Xhat)
    }
    
    SigmaSSmat <- map(blocks, ~rep(list(SigmaSS), .x)) 
    Xhat <- map2(SigmaSSmat, indices, ~ImputeBlock(Wmat, .x, .y, Sigma_XX, Sigma_WW)) %>%
      invoke(rbind, .) 
    data <- mutate(data, Xhat.cor1 = Xhat[,1], Xhat.cor2 = Xhat[,2], Xhat.cor3 = Xhat[,3])
  }
  if (estimate == "Truth") {
    Wmat <- as.matrix(data[, c("W1", "W2", "W3")])
    Sigma_WW <- cor.mat.w
    Sigma_XX <- cor.mat.x
    kappa <- solve(Sigma_WW) %*% Sigma_XX
    EX <- Wmat %*% t(kappa)
    data <- mutate(data, EX1 = EX[,1], EX2 = EX[,2], EX3 = EX[,3])
  }
  return(data)
}

EstimateSBW <- function(balance_variable, data, target, re_spec) {
  Zs <- paste0(sprintf("Z%s", 1:length(target)), collapse = "+")
  formula_string <- paste0("~ ", Zs)
  formula <- as.formula(gsub("Z", balance_variable, formula_string))
  groups <- as.numeric(table(data$id))
  s2_spec <- 1 - re_spec
  model <- optweight.svy(formula, tols = 0, data = data, target = target,
                         re = re_spec, sigma2.y = s2_spec, 
                         group_n = groups)
  return(model$weights)
}

SimulationHSBW <- function(population, sample_states, jackknife_varest = TRUE) {
  
  popdat <- population$data
  re_true <- population$rho
  var.y <- population$var.y
  re_list <- c(0, 0.25, 0.5)
  
  cor.mat.w <- population$cor.mat.w
  cor.mat.x <- population$cor.mat.x
  
  ett.pop <- 3
  
  data <- popdat %>%
    nest(-id) %>%
    sample_n(sample_states) %>%
    unnest() %>%
    arrange(id)
  
  blocks <- as.numeric(table(data$id))
  finish <- cumsum(blocks)
  start <- lag(finish) + 1; start[1] <- 1
  indices <- map2(start, finish, ~c(.x:.y))
  
  SigmaSS_list <- CalcSigSS(data)
  weights <- unlist(map(indices, length))
  weights <- weights/sum(weights)
  SigmaSS <- Reduce(`+`, map2(SigmaSS_list, weights, ~.x * .y))
  
  cor.mat.v.list <- map(1:nrow(data), ~data[.x, c("V1", "V2", "V3")] %>%
                          as.numeric(.[1,]) %>%
                          diag())
  
  Sigma_vv_est <- map(cor.mat.v.list, ~.x + diag(rnorm(3, 0, sqrt(0.001*450))))
  
  data <- data %>%
    CalcXhat(estimate = "Truth", cor.mat.x, cor.mat.w, Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Homogeneous", Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Heterogeneous", Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Correlated", SigmaSS = SigmaSS, Sigma_vv_est = Sigma_vv_est) %>%
    arrange(id)
  
  target <- c(1, 1, 1)
  
  weights.xhat.hom  <- map(re_list, ~EstimateSBW("Xhat.hom", data, target, .x))
  weights.xhat.het  <- map(re_list, ~EstimateSBW("Xhat.het", data, target, .x))
  weights.xhat.cor  <- map(re_list, ~EstimateSBW("Xhat.cor", data, target, .x))
  weights.ex        <- map(re_list, ~EstimateSBW("EX", data, target, .x))
  weights.x         <- map(re_list, ~EstimateSBW("X", data, target, .x))
  weights.w         <- map(re_list, ~EstimateSBW("W", data, target, .x))
  
  GenEsts <- function(data, weights.w, weights.x, weights.ex, weights.xhat.hom, weights.xhat.het,
                      weights.xhat.cor) {
    if (!is.data.frame(data)) {
      ests.xhat.hom <- map(data, ~sum(.x$J*weights.xhat.hom)/sum(weights.xhat.hom))
      ests.xhat.het <- map(data, ~sum(.x$J*weights.xhat.het)/sum(weights.xhat.het))
      ests.xhat.cor <- map(data, ~sum(.x$J*weights.xhat.cor)/sum(weights.xhat.cor))
      ests.ex   <- map(data, ~sum(.x$J*weights.ex)/sum(weights.ex))
      ests.x    <- map(data, ~sum(.x$J*weights.x)/sum(weights.x))
      ests.w    <- map(data, ~sum(.x$J*weights.w)/sum(weights.w))
      return(list(ests.x = ests.x, ests.w = ests.w, ests.ex = ests.ex,
                  ests.xhat.het = ests.xhat.het, ests.xhat.hom = ests.xhat.hom,
                  ests.xhat.cor = ests.xhat.cor))
    }
    if (is.data.frame(data)) {
      ests.xhat.hom <- sum(data$J*weights.xhat.hom)/sum(weights.xhat.hom)
      ests.xhat.het <- sum(data$J*weights.xhat.het)/sum(weights.xhat.het)
      ests.xhat.cor <- sum(data$J*weights.xhat.cor)/sum(weights.xhat.cor)
      ests.ex   <- sum(data$J*weights.ex)/sum(weights.ex)
      ests.x    <- sum(data$J*weights.x)/sum(weights.x)
      ests.w    <- sum(data$J*weights.w)/sum(weights.w)
      return(tibble(ests = c(ests.x, ests.w, ests.ex, ests.xhat.het, ests.xhat.hom, ests.xhat.cor),
                    Xset = c("X", "W", "EX", "Xhat.het", "Xhat.hom", "Xhat.cor")))
    }
  }
  
  GenEstsP <- partial(GenEsts, data = data)
  
  all_weights <- pmap(list(weights.x = weights.x, 
                           weights.ex = weights.ex, 
                           weights.w = weights.w, 
                           weights.xhat.hom = weights.xhat.hom,
                           weights.xhat.het = weights.xhat.het,
                           weights.xhat.cor = weights.xhat.cor), GenEstsP)
  
  transpose_results <- function(results, xvar, re_list) {
    res <- map(results, ~filter(.x, Xset == xvar)) %>%
      invoke(rbind, .) %>%
      mutate(re = re_list)
    return(res)
  }
  
  if (jackknife_varest == TRUE) {
    
    all_weights_t <- map(c("X", "W", "EX", "Xhat.hom", "Xhat.het", "Xhat.cor"), 
                         ~transpose_results(all_weights, .x, re_list))
    
    JackknifeP <- partial(Jackknife, data = data, target = target, cor.mat.w = cor.mat.w, cor.mat.x = cor.mat.x,
                          SigmaSS_list = SigmaSS_list, Sigma_vv_est = Sigma_vv_est)
    
    res <- pmap(list(estimate_df = all_weights_t, 
                     #jackest = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                     jackest = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
                     estimate = c(rep("Homogeneous", 4), "Heterogeneous", "Correlated")), JackknifeP) %>% 
      invoke(rbind, .) %>%
      mutate(N = nrow(data))
  }
  if (jackknife_varest == FALSE) {
    res <- all_weights %>%
      map2(re_list, ~mutate(.x, re = .y)) %>%
      invoke(rbind, .) %>%
      mutate(truth = 3, 
             bias = ests - 3,
             mse = (ests - 3)^2,
             N = nrow(data))
  }
  return(res)
}

SimulationHSBW.cor <- function(population, sample_states) {
  
  popdat <- population$data
  re_true <- population$rho
  var.y <- population$var.y
  re_list <- c(0, 0.25, 0.5)
  
  cor.mat.w <- population$cor.mat.w
  cor.mat.x <- population$cor.mat.x
  
  ett.pop <- 3
  
  data <- popdat %>%
    nest(-id) %>%
    sample_n(sample_states) %>%
    unnest()
  
  blocks <- as.numeric(table(data$id))
  finish <- cumsum(blocks)
  start <- lag(finish) + 1; start[1] <- 1
  indices <- map2(start, finish, ~c(.x:.y))
  
  SigmaSS_list <- CalcSigSS(data)
  weights <- unlist(map(indices, length))
  weights <- weights/sum(weights)
  SigmaSS <- Reduce(`+`, map2(SigmaSS_list, weights, ~.x * .y))
  
  cor.mat.v.list <- map(1:nrow(data), ~data[.x, c("V1", "V2", "V3")] %>%
                          as.numeric(.[1,]) %>%
                          diag())
  
  Sigma_vv_est <- map(cor.mat.v.list, ~.x + diag(rnorm(3, 0, sqrt(0.001*450))))
  
  data <- data %>%
    CalcXhat(estimate = "Truth", cor.mat.x, cor.mat.w, Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Homogeneous", Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Heterogeneous", Sigma_vv_est = Sigma_vv_est) %>%
    CalcXhat(estimate = "Correlated", SigmaSS = SigmaSS, Sigma_vv_est = Sigma_vv_est) %>%
    arrange(id)
  
  target <- c(1, 1, 1)
  
  weights.xhat.cor  <- map(re_list, ~EstimateSBW("Xhat.cor", data, target, .x))
  weights.xhat.hom  <- map(re_list, ~EstimateSBW("Xhat.het", data, target, .x))
  weights.xhat.het  <- map(re_list, ~EstimateSBW("Xhat.hom", data, target, .x))
  
  GenEsts <- function(data, weights.xhat.cor, weights.xhat.hom, weights.xhat.het) {
    if (!is.data.frame(data)) {
      ests.xhat.cor <- map(data, ~sum(.x$J*weights.xhat.cor)/sum(weights.xhat.cor))
      ests.xhat.hom <- map(data, ~sum(.x$J*weights.xhat.hom)/sum(weights.xhat.hom))
      ests.xhat.het <- map(data, ~sum(.x$J*weights.xhat.het)/sum(weights.xhat.het))
      return(list(ests.xhat.cor = ests.xhat.cor,
                  ests.xhat.hom = ests.xhat.hom,
                  ests.xhat.het = ests.xhat.het))
    }
    if (is.data.frame(data)) {
      ests.xhat.cor <- sum(data$J*weights.xhat.cor)/sum(weights.xhat.cor)
      ests.xhat.hom <- sum(data$J*weights.xhat.hom)/sum(weights.xhat.hom)
      ests.xhat.het <- sum(data$J*weights.xhat.het)/sum(weights.xhat.het)
      
      return(tibble(ests = c(ests.xhat.cor, ests.xhat.hom, ests.xhat.het),
                    Xset = c("Xhat.cor", "Xhat.hom", "Xhat.het")))
    }
  }
  
  GenEstsP <- partial(GenEsts, data = data)
  
  all_weights <- pmap(list(weights.xhat.cor = weights.xhat.cor,
                           weights.xhat.hom = weights.xhat.hom,
                           weights.xhat.het = weights.xhat.het), GenEstsP)
  
  transpose_results <- function(results, xvar, re_list) {
    res <- map(results, ~filter(.x, Xset == xvar)) %>%
      invoke(rbind, .) %>%
      mutate(re = re_list)
    return(res)
  }
  
  res <- all_weights %>%
    map2(re_list, ~mutate(.x, re = .y)) %>%
    invoke(rbind, .) %>%
    mutate(truth = 3, 
           bias = ests - 3,
           mse = (ests - 3)^2,
           N = nrow(data))
  return(res)
}

Jackknife <- function(data, estimate_df, target, jackest, cor.mat.w, cor.mat.x, 
                      SigmaSS_list, Sigma_vv_est, estimate = NULL) {
  
  re_list <- estimate_df$re
  estimates <- estimate_df$ests
  balance_variable <- unique(estimate_df$Xset)
  
  id_list <- unique(data$id)
  M <- length(id_list)
  groups <- as.numeric(table(data$id))
  data_list <- map(id_list, ~filter(data, id != .x))
  group_list <- map(1:length(groups), ~groups[-.x])
  
  blocks <- as.numeric(table(data$id))
  finish <- cumsum(blocks)
  start <- lag(finish) + 1; start[1] <- 1
  indices <- map2(start, finish, ~c(.x:.y)) 
  index_list <- map(1:length(indices), ~indices[-.x])
  weight_list <- map(index_list, ~unlist(map(.x, length))) %>%
    map(~.x/sum(.x))
  index_jack <- map(index_list, unlist)
  Sigma_vv_est_r <- map(index_jack, ~Sigma_vv_est[.x])
  SigmaSS_list_r <- map(1:length(SigmaSS_list), ~SigmaSS_list[-.x])
  SigmaSS_list_r <- map2(weight_list, SigmaSS_list_r, ~Reduce(`+`, map2(.y, .x, ~.x * .y)))
  
  if (jackest == TRUE) {
    CalcXhatp <- partial(CalcXhat, estimate = estimate, cor.mat.x = cor.mat.x, cor.mat.w = cor.mat.w)
    data_list <- pmap(list(data = data_list, 
                           SigmaSS = SigmaSS_list_r,
                           Sigma_vv_est = Sigma_vv_est_r), CalcXhatp)
  }
  
  IterDatList <- function(data_list, Xset, target, re_spec) {
    wres <- map(data_list, ~EstimateSBW(Xset, .x, target, re_spec))
    res <- map2(data_list, wres, ~sum(.x$J*.y)/sum(.y))
    return(res)
  }
  
  est_list <- map(re_list, ~IterDatList(data_list, balance_variable, target, .x))
  mean_est <- map(est_list, ~mean(unlist(.x)))
  var_est  <- map2(est_list, mean_est, ~ ((M - 1)/M) * sum((unlist(.x) - .y)^2))
  
  res <- tibble(
    truth = rep(3, length(est_list)),
    est = estimates,
    bias = estimates - truth,
    se = sqrt(unlist(var_est)),
    b.est = unlist(mean_est),
    lci = est - 1.96*se,
    uci = est + 1.96*se,
    covered = ifelse(truth > lci & truth < uci, 1, 0),
    re = re_list,
    Xset = balance_variable,
    jackest = jackest
  )
  return(res)
}

RunSims <- function(num_sims, pop_size, sample_states, X_cor, SNR.xw, rho, SS_cor, hetero) {
  population <- GenerateHierPopulation(number_states = pop_size, 
                                       X_cor = X_cor, 
                                       SNR.xw = SNR.xw, 
                                       rho = rho, 
                                       SS_cor = SS_cor, 
                                       hetero = hetero)
  
  sim <- map(1:num_sims, ~SimulationHSBW(population, sample_states, jackknife_varest = TRUE))
  return(sim)
}


RunSimsX <- function(num_sims, pop_size, sample_states, X_cor, SNR.xw, rho, SS_cor, hetero) {
  population <- GenerateHierPopulation(number_states = pop_size, 
                                       X_cor = X_cor, 
                                       SNR.xw = SNR.xw, 
                                       rho = rho, 
                                       SS_cor = SS_cor, 
                                       hetero = hetero)
  
  sim <- map(1:num_sims, ~SimulationHSBW(population, sample_states, jackknife_varest = TRUE,
                                         x.only = TRUE))
  return(sim)
}

RunCorSims <- function(population, num_states, nsims) {
  res <- map(1:nsims, ~SimulationHSBW.cor(population, num_states))
  return(res)
}
