# program: 07-validation.R
# purpose: conduct validation study on pre-treatment outcomes
# author: max rubinstein
# date modified: january 14, 2022

source("03_Analysis/02_model-estimation.R")
library(xtable); library(gt)

# placebo tests ----------------------------------------------------------------------------------------
estimates <- function(data, variables, tols) {
  tdat <- data %>%
    map(~filter(.x, treatment == 1))
  cdat <- subset(data[[3]], treatment == 0)
  targets <- colMeans(cdat[variables])
  w1 <- generate_weight_list(tols, tdat, targets, linf_imbalance, 
                             stop_criterion = 0.5, max_iter = 1e7)
  w1
}

comparison1 <- function(data, weight_result) {
  predict_oos <- function(tdat, weights) {
    map2(tdat, weights, ~mutate(.x, weights = .y$weights) %>%
           dplyr::summarize(p12 = sum(weights*hins_unins_pct_2012)/sum(weights),
                            p13 = sum(weights*hins_unins_pct_2013)/sum(weights),
                            p14 = sum(weights*hins_unins_pct_2014)/sum(weights)))
  }
  tdat <- data %>%
    map(~filter(.x, treatment == 1))
  
  map(weight_result, ~predict_oos(tdat, .x))
}

comparison2 <- function(data, weight_result) {
  predict_oos <- function(tdat, weights) {
    map2(tdat, weights, ~mutate(.x, weights = .y$weights) %>%
           dplyr::summarize(p13 = sum(weights*hins_unins_pct_2013)/sum(weights),
                            p14 = sum(weights*hins_unins_pct_2014)/sum(weights),
                            p09 = sum(weights*hins_unins_pct_2009)/sum(weights)))
  }
  tdat <- data %>%
    map(~filter(.x, treatment == 1))
  
  map(weight_result, ~predict_oos(tdat, .x))
}

comparison3 <- function(data, weight_result) {
  predict_oos <- function(tdat, weights) {
    map2(tdat, weights, ~mutate(.x, weights = .y$weights) %>%
           dplyr::summarize(p14 = sum(weights*hins_unins_pct_2014)/sum(weights),
                            p10 = sum(weights*hins_unins_pct_2010)/sum(weights),
                            p09 = sum(weights*hins_unins_pct_2009)/sum(weights)))
  }
  tdat <- data %>%
    map(~filter(.x, treatment == 1))
  
  map(weight_result, ~predict_oos(tdat, .x))
}

read_vars <- function(var_name) {
  read_csv("../02_Specs/tol_specs.csv") %>%
    filter(`Reference Variable` == 0) %>%
    arrange(Variable) %>%
    .[[var_name]] %>%
    sort()
}

cov_errors <- function(data, weight_results, variables) {
  tdat <- data %>%
    map(~filter(.x, treatment == 1))
  cdat <- subset(data[[3]], treatment == 0)[variables] 
  prey <- variables[grep("hins_unins", variables)]
  targets <- colMeans(cdat[variables])
  preymeans <- colMeans(cdat[prey])
  
  predict_oos <- function(data, weight_results, variables) {
    map2(data, weight_results, ~mutate(.x, weights = .y$weights) %>%
           dplyr::summarize_at(vars(all_of(variables)), ~sum(.*weights)/sum(weights))) %>%
      invoke(rbind, .)
  }
  map(weight_results, ~predict_oos(tdat, .x, variables)) %>%
    map(~t(.)) %>%
    map(~cbind(.x, targets))
}

# set names ---------------------------------------------------------------------------------------
covariate_group <- c("Unadjusted")

sigma_estimator <- c("sigma_uu_i", "sigma_uu_avg", "sigma_zero")

estimator <- c("SBW", "H-SBW", "BC-SBW", "BC-HSBW")

variables <- read_vars("Variable")
vars1 <- read_vars("var_valid")
vars2 <- read_vars("var_test")

pretxy <- variables[grep("unins", variables)]

tol_list <- read_csv("../02_Specs/tol_specs.csv") 

tol_list <- map(0, ~tol_list %>% mutate(`Base Tol` = if_else(grepl(.x, Group), 100, `Base Tol`))) %>%
  map(~filter(.x, Variable %in% variables)) %>%
  map(~arrange(.x, Variable)) %>%
  map(~.x$`Base Tol`) %>%
  map(~set_names(.x, variables))

tols1 <- tol_list[[1]]
names(tols1) <- vars1
tols2 <- tols1
names(tols2) <- vars2

imputed_dat_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.rds") %>%
  unnest() %>%
  nest(-key, -set)

imputed_dat_c2 <- readRDS("../01_ProcessedData/calibrated-data-c2.rds") %>%
  unnest() %>%
  nest(-key, -set)

# estimate HSBW weights --------------------------------------------------------
weights1 <- estimates(imputed_dat_c1$data[7:9], vars1, tols1)
weights2 <- estimates(imputed_dat_c1$data[4:6], vars2, tols2)
weights3 <- estimates(imputed_dat_c1$data[1:3], variables, tol_list[[1]])

weights1.c2 <- estimates(imputed_dat_c2$data[7:9], vars1, tols1)
weights2.c2 <- estimates(imputed_dat_c2$data[4:6], vars2, tols2)
weights3.c2 <- estimates(imputed_dat_c2$data[1:3], variables, tol_list[[1]])

dat1 <- imputed_dat_c1$data[7:9]
dat2 <- imputed_dat_c1$data[4:6]
dat3 <- imputed_dat_c1$data[1:3]

dat1.c2 <- imputed_dat_c2$data[7:9]
dat2.c2 <- imputed_dat_c2$data[4:6]
dat3.c2 <- imputed_dat_c2$data[1:3]

cdat <- subset(dat1[[3]], treatment == 0)
tdat <- subset(dat1[[3]], treatment == 1)
cdat.c2 <- subset(dat1.c2[[3]], treatment == 0)
tdat.c2 <- subset(dat1.c2[[3]], treatment == 1)

py <- c("hins_unins_pct_2009", "hins_unins_pct_2010",
        "hins_unins_pct_2011", "hins_unins_pct_2012",
        "hins_unins_pct_2013", "hins_unins_pct_2014")

truth <- colMeans(cdat[, py])
names(truth) <- c("2009", "2010", "2011", "2012", "2013", "2014")
txtruth <- colMeans(tdat[, py])
txtruth.c2 <- colMeans(tdat.c2[, py])

trends_table <- round(rbind(truth, txtruth, txtruth.c2), 2) %>%
  as_tibble() %>%
  mutate(Treatment = c("Non-expansion", "Expansion (primary dataset)", "Expansion (early excluded)")) %>%
  select(Treatment, everything())

print(xtable::xtable(trends_table,caption = "Mean non-elderly adult uninsurance rates, 2009-2014"), 
      type="latex", caption.placement = "top",
      include.rownames = FALSE)

# estimate GLS weights ---------------------------------------------------------
CalculateGLSEsts <- function(data_list) {
  sigma <- as.numeric(table(data_list$data[[1]]$state[data_list$data[[1]]$treatment == 1]))
  sigma <- Matrix::bdiag(map(sigma, ~matrix(rep(1/6, .x^2), .x, .x)))
  diag(sigma) <- 1
  sigma <- Matrix::solve(sigma)
  sigma <- as.matrix(sigma)
  
  res <- data_list %>%
    mutate(outcome = map2(data, rep(sprintf("hins_unins_pct_201%s", 4:2), each = 3), ~.x[[.y]]),
           variables = list(variables, variables, variables, vars2, vars2, vars2, vars1, vars1, vars1),
           txdat = map2(data, variables, ~as.matrix(cbind(1, filter(.x, treatment == 1) %>% select(.y)))),
           ctdat = map2(data, variables, ~as.matrix(cbind(1, filter(.x, treatment == 0) %>% select(.y)))),
           ctmeans = map2(ctdat, variables, ~c(1, as.matrix(colMeans(.x[,.y])))),
           weights.ols = map2(txdat, ctmeans, ~t(t(.y) %*% solve(t(.x) %*% .x) %*% t(.x))),
           weights.gls = map2(txdat, ctmeans, ~t(t(.y) %*% solve(t(.x) %*% sigma %*% .x) %*% t(.x) %*% sigma)),
           outcomes.1 = map2(data, outcome, ~.y[.x$treatment == 1]),
           outcomes.0 = map2(data, outcome, ~mean(.y[.x$treatment == 0])),
           estimate.ols = map2(weights.ols, outcomes.1, ~sum(.x*.y)),
           estimate.gls = map2(weights.gls, outcomes.1, ~sum(.x*.y)),
           OLS = map2(estimate.ols, outcomes.0, ~.x - .y),
           GLS = map2(estimate.gls, outcomes.0, ~.x - .y)) %>%
    select(key, set, OLS, GLS) %>% 
    unnest() 
  
  return(res)
}

gls.weights <- CalculateGLSEsts(imputed_dat_c1)
gls.weights.c2 <- CalculateGLSEsts(imputed_dat_c2)

# comparison tables ------------------------------------------------------------
GenCompTable <- function(dat_list, weight_list) {
  dat1 <- dat_list[[1]]; weights1 <- weight_list[[1]]
  dat2 <- dat_list[[2]]; weights2 <- weight_list[[2]]
  dat3 <- dat_list[[3]]; weights3 <- weight_list[[3]]
  
  c1 <- comparison1(dat1, weights1) %>%
    map(~invoke(rbind, .x)) %>%
    map(~mutate(.x, sigma_estimate = sigma_estimator)) %>%
    invoke(rbind, .) %>%
    mutate(estimator = rep(estimator, each = 3)) %>%
    mutate(truth12 = truth["2012"], truth13 = truth["2013"], truth14 = truth["2014"]) %>%
    mutate(err.12 = p12 - truth12, err.13 = p13 - truth13, err.14 = p14 - truth14)
  
  c2 <- comparison2(dat2, weights2) %>%
    map(~invoke(rbind, .x)) %>%
    map(~mutate(.x, sigma_estimate = sigma_estimator)) %>%
    invoke(rbind, .) %>%
    mutate(estimator = rep(estimator, each = 3)) %>%
    mutate(truth13 = truth["2013"], truth14 = truth["2014"], truth09 = truth["2009"]) %>%
    mutate(err.13 = p13 - truth13, err.14 = p14 - truth14, err.09 = p09 - truth09) 
  
  c3 <- comparison3(dat3, weights3) %>%
    map(~invoke(rbind, .x)) %>%
    map(~mutate(.x, sigma_estimate = sigma_estimator)) %>%
    invoke(rbind, .) %>%
    mutate(estimator = rep(estimator, each = 3)) %>%
    mutate(truth14 = truth["2014"], truth10 = truth["2010"], truth09 = truth["2009"]) %>%
    mutate(err.14 = p14 - truth14, err.10 = p10 - truth10, err.09 = p09 - truth09) 
  
  list(c1 = c1, c2 = c2, c3 = c3)
}

GenFinTab <- function(table_list, gls_results) {
  c1a <- table_list$c1 %>%
    select(sigma_estimate, estimator, err.12)
  
  c2a <- table_list$c2 %>%
    select(sigma_estimate, estimator, err.13)
  
  fintab <- left_join(c1a, c2a, by = c("sigma_estimate", "estimator")) %>%
    mutate(rmse = sqrt((err.12^2 + err.13^2)/2)) %>%
    arrange(rmse) %>%
    filter(!grepl("SBW1", estimator)) %>%
    mutate_at("sigma_estimate", ~stringr::str_replace_all(.,
                                                          c("sigma_uu_avg" = "Homogeneous",
                                                            "sigma_uu_i"   = "Heterogeneous",
                                                            "sigma_zero"   = "Unadjusted"))) %>%
    select(`Sigma estimate` = sigma_estimate, Estimator = estimator,
           `2012 error` = err.12, `2013 error` = err.13, RMSE = rmse) 
  
  gls.tab <- gls_results %>%
    mutate_at("key", ~stringr::str_replace_all(.,c("sigma_uu_avg" = "Homogeneous",
                                                   "sigma_uu_i_modeled" = "Heterogeneous",
                                                   "sigma_zero" = "Unadjusted"))) %>%
    filter(set != "true") %>%
    mutate_at("set", ~stringr::str_replace_all(., c("valid" = "2012 error",
                                                    "test"  = "2013 error"))) %>%
    gather(Estimator, value, OLS, GLS) %>%
    spread(set, value) %>%
    rename(`Sigma estimate` = key) %>%
    mutate(RMSE = sqrt(0.5 * (`2012 error`^2 + `2013 error`^2)))
  
  fintab <- fintab %>%
    bind_rows(gls.tab) %>%
    arrange(RMSE)
  
  return(fintab)
}

tables_c1 <- GenCompTable(list(dat1, dat2, dat3), list(weights1, weights2, weights3))
tables_c2 <- GenCompTable(list(dat1.c2, dat2.c2, dat3.c2), list(weights1.c2, weights2.c2, weights3.c2))

fintab <- GenFinTab(tables_c1, gls.weights)
fintab.c2 <- GenFinTab(tables_c2, gls.weights.c2)

# output tables for paper ------------------------------------------------------
tot.tab <- fintab %>%
  left_join(
    fintab.c2 %>%
      rename(`2012 error1` = `2012 error`,
             `2013 error1` = `2013 error`,
             RMSE1 = RMSE)
  )

final <- tot.tab %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  gt() %>%
  tab_spanner(
    label = "Primary data",
    columns = c("2012 error", "2013 error", "RMSE")
  ) %>%
  tab_spanner(
    label = "Early expansion \n excluded",
    columns = c("2012 error1", "2013 error1", "RMSE1")
  ) %>%
  cols_label(
    `2012 error1` = "2012 error",
    `2013 error1` = "2013 error",
    RMSE = "RMSE"
  ) %>%
  as_latex()

final.alt <- bind_cols(
  tot.tab %>%
    mutate(`Mean Error` = (`2012 error` + `2013 error`) / 2) %>%
    select(`Sigma estimate1` = `Sigma estimate`, Estimator1 = Estimator, `Mean Error`, RMSE) %>%
    filter(!grepl("GLS|OLS", Estimator1)) %>%
    arrange(RMSE),
  tot.tab %>%
    mutate(`Mean Error1` = (`2012 error1` + `2013 error1`) / 2) %>%
    select(`Sigma estimate`, Estimator, `Mean Error1`, RMSE1) %>%
    filter(!grepl("GLS|OLS", Estimator)) %>%
    arrange(RMSE1)
)

final.alt.full <- bind_cols(
  tot.tab %>%
    mutate(`Mean Error` = (`2012 error` + `2013 error`) / 2) %>%
    select(`Sigma estimate1` = `Sigma estimate`, Estimator1 = Estimator, `Mean Error`, RMSE) %>%
    arrange(RMSE),
  tot.tab %>%
    mutate(`Mean Error1` = (`2012 error1` + `2013 error1`) / 2) %>%
    select(`Sigma estimate`, Estimator, `Mean Error1`, RMSE1) %>%
    arrange(RMSE1)
)

final.alt.gt <- final.alt %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  gt() %>%
  tab_spanner(
    label = "Primary data",
    columns = c("Sigma estimate1", "Estimator1", "Mean Error", "RMSE")
  ) %>%
  tab_spanner(
    label = "Early expansion \n excluded",
    columns = c("Sigma estimate", "Estimator", "Mean Error1", "RMSE1")
  ) %>%
  cols_label(
    `Sigma estimate1` = "Sigma estimate",
    `Estimator1` = "Estimator",
    RMSE1 = "RMSE",
    `Mean Error1` = "Mean Error"
  ) %>%
  as_latex()

final.alt.full.gt <- final.alt.full %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  gt() %>%
  tab_spanner(
    label = "Primary data",
    columns = c("Sigma estimate1", "Estimator1", "Mean Error", "RMSE")
  ) %>%
  tab_spanner(
    label = "Early expansion \n excluded",
    columns = c("Sigma estimate", "Estimator", "Mean Error1", "RMSE1")
  ) %>%
  cols_label(
    `Sigma estimate1` = "Sigma estimate",
    `Estimator1` = "Estimator",
    RMSE1 = "RMSE",
    `Mean Error1` = "Mean Error"
  ) %>%
  as_latex() 

rmse.order <- fintab %>%
  filter(!Estimator %in% c("OLS", "GLS")) %>%
  mutate(order = paste(`Sigma estimate`, Estimator)) %>%
  .$order

rmse.order.c2 <- fintab.c2 %>%
  filter(!Estimator %in% c("OLS", "GLS")) %>%
  mutate(order = paste(`Sigma estimate`, Estimator)) %>%
  .$order

# ordered by RMSE on primary data
print(xtable::xtable(fintab.c2),
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

final %>%
  as.character() %>%
  cat()

# mean error / RMSE versions
final.alt.gt %>%
  as.character() %>%
  cat()

final.alt.full.gt %>%
  as.character() %>%
  cat()


# additional plots -------------------------------------------------------------
FullWeightPlot <- function(data, weight_list, adjustment_set) {
  data %>%
    filter(treatment == 1) %>%
    mutate(SBW = weight_list$SBW[[adjustment_set]]$weights,
           `H-SBW` = weight_list$HSBW[[adjustment_set]]$weights,
           `BC-SBW` = weight_list$`BC-SBW`[[adjustment_set]]$weights,
           `BC-HSBW` = weight_list$`BC-HSBW`[[adjustment_set]]$weights) %>%
    select(cpuma, state, SBW:`BC-HSBW`) %>%
    gather(weight, weight_value, SBW:`BC-HSBW`) %>%
    group_by(weight) %>%
    mutate(weight_value = 100 * (weight_value / sum(weight_value))) %>%
    mutate_at("weight_value", 
              funs(Positive = if_else(. > 0, ., 0),
                   Negative = if_else(. <= 0, -1*., 0))) %>%
    select(state, weight, Positive, Negative) %>%
    gather(`Weight sign`, value, -state, -weight) %>%
    group_by(state, weight, `Weight sign`) %>% 
    dplyr::summarize(value = sum(value)) %>%
    mutate(value = if_else(`Weight sign` == "Positive", value, -value)) %>%
    mutate_at("state", factor) %>%
    mutate_at("weight", ~factor(., levels = c("H-SBW", "BC-HSBW", "SBW", "BC-SBW"))) %>%
    ggplot(aes(x = state, y = value, fill = `Weight sign`)) +
    geom_bar(stat = "identity") +
    facet_wrap(~weight) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    scale_x_discrete(limits = levels("state")) +
    xlab("") + ylab("Distribution of positive and negative weights within state – total sum of each category")
}

PaperWeightPlot <- function(data, weight_list, adjustment_set) {
  data %>%
    filter(treatment == 1) %>%
    mutate(SBW = weight_list$SBW[[adjustment_set]]$weights,
           `H-SBW` = weight_list$HSBW[[adjustment_set]]$weights,
           `BC-HSBW` = weight_list$`BC-HSBW`[[adjustment_set]]$weights) %>%
    select(cpuma, state, SBW:`BC-HSBW`) %>%
    gather(weight, weight_value, SBW:`BC-HSBW`) %>%
    group_by(weight) %>%
    mutate(weight_value = 100 * (weight_value / sum(weight_value))) %>%
    mutate_at("weight_value", 
              funs(Positive = if_else(. > 0, ., 0),
                   Negative = if_else(. <= 0, -1*., 0))) %>%
    select(state, weight, Positive, Negative) %>%
    gather(`Weight sign`, value, -state, -weight) %>%
    group_by(state, weight, `Weight sign`) %>% 
    dplyr::summarize(value = sum(value)) %>%
    mutate(value = if_else(`Weight sign` == "Positive", value, -value)) %>%
    mutate_at("state", factor) %>%
    mutate_at("weight", ~factor(., levels = c("H-SBW", "BC-HSBW", "SBW"))) %>%
    ggplot(aes(x = state, y = value, fill = `Weight sign`)) +
    geom_bar(stat = "identity") +
    facet_wrap(~weight) +
    theme_minimal() +
    scale_fill_grey() +
    coord_flip() +
    scale_x_discrete(limits = levels("state")) +
    xlab("") + ylab("Distribution of positive and negative weights within state – total sum of each category")
}

weight.plot.c1 <- FullWeightPlot(imputed_dat_c1$data[[2]], weights3, 2)  
weight.plot.c2 <- FullWeightPlot(imputed_dat_c2$data[[2]], weights3.c2, 2)  
weight.plot.c1.paper <- PaperWeightPlot(imputed_dat_c1$data[[1]], weights3, 2)

ggsave("../../02_Paper/01_Plots/weights-by-state-sbw-hsbw-c1.png", weight.plot.c1.paper,
       width = 10, height = 6, type = "cairo")
ggsave("../../02_Paper/01_Plots/weights-by-state-c1-all.png", weight.plot.c1,
       width = 10, height = 6, type = "cairo")
ggsave("../../02_Paper/01_Plots/weights-by-state-c2-all.png", weight.plot.c2,
       width = 10, height = 6, type = "cairo")

# additional investigations ----------------------------------------------------
Y2013.1 <- imputed_dat_c1$data[[3]]$hins_unins_pct_2013[imputed_dat_c1$data[[3]]$treatment == 1]
sbw     <- map(weights3$SBW, ~sum(.x$weights*Y2013.1)/sum(.x$weights)) %>% unlist()
hsbw    <- map(weights3$HSBW, ~sum(.x$weights*Y2013.1)/sum(.x$weights)) %>% unlist()
bc.sbw   <- map(weights3$`BC-SBW`, ~sum(.x$weights*Y2013.1)/sum(.x$weights)) %>% unlist()
bc.hsbw  <- map(weights3$`BC-HSBW`, ~sum(.x$weights*Y2013.1)/sum(.x$weights)) %>% unlist()

dat <- imputed_dat_c1$data[[3]] %>%
  filter(treatment == 1) %>%
  select(cpuma, state, J_2013 = hins_unins_pct_2013, J_2014 = hins_unins_pct_2014) %>%
  mutate(
    sbw_2014_hom = weights3$SBW[[2]]$weights,
    sbw_2014_non = weights3$SBW[[3]]$weights,
    sbw_2013_hom = weights2$SBW[[2]]$weights,
    sbw_2013_non = weights2$SBW[[3]]$weights,
    bc.sbw_2014_hom = weights3$`BC-SBW`[[3]]$weights,
    bc.sbw_2014_non = weights3$`BC-SBW`[[2]]$weights,
    bc.sbw_2013_hom = weights2$`BC-SBW`[[3]]$weights,
    bc.sbw_2013_non = weights2$`BC-SBW`[[2]]$weights
  ) 

dat %>%
  group_by(state) %>%
  mutate(np = n()/nrow(.)) %>%
  ungroup() %>%
  gather(weights, weight_value, sbw_2014_hom:bc.sbw_2013_non) %>%
  mutate(extrap = ifelse(weight_value < 0, 1, 0)) %>%
  filter(extrap == 1) %>%
  group_by(weights) %>%
  mutate(weight_value = weight_value/sum(weight_value)) %>%
  group_by(weights, extrap, state, np) %>%
  summarize(weight_value = sum(weight_value)) %>%
  filter(extrap == 1, weights == "bc.sbw_2014_hom") %>%
  arrange(-weight_value)

dat.c2 <- imputed_dat_c2$data[[3]] %>%
  filter(treatment == 1) %>%
  select(cpuma, state, J_2013 = hins_unins_pct_2013, J_2014 = hins_unins_pct_2014) %>%
  mutate(
    sbw_2014_hom = weights3.c2$SBW[[2]]$weights,
    sbw_2014_non = weights3.c2$SBW[[3]]$weights,
    sbw_2013_hom = weights2.c2$SBW[[2]]$weights,
    sbw_2013_non = weights2.c2$SBW[[3]]$weights,
    bc.sbw_2014_hom = weights3.c2$`BC-SBW`[[3]]$weights,
    bc.sbw_2014_non = weights3.c2$`BC-SBW`[[2]]$weights,
    bc.sbw_2013_hom = weights2.c2$`BC-SBW`[[3]]$weights,
    bc.sbw_2013_non = weights2.c2$`BC-SBW`[[2]]$weights
  ) 

itab.1 <- map(tables_c1, ~.x %>%
            select(sigma_estimate, estimator, txfx_2014 = p14)) %>%
  map2(c("2009-2011", "2010-2012", "2011-2013"), ~mutate(.x, covariates = .y)) %>%
  invoke(rbind, .) %>%
  spread(covariates, txfx_2014) %>%
  arrange(rev(estimator)) %>%
  mutate_at("sigma_estimate", ~stringr::str_replace_all(.,
                                                        c("sigma_uu_avg" = "Homogeneous",
                                                          "sigma_uu_i"   = "Heterogeneous",
                                                          "sigma_zero"   = "Unadjusted"))) %>%
  rename(Adjustment = sigma_estimate, Estimator = estimator) %>%
  mutate(order = paste(Adjustment, Estimator)) %>%
  mutate_at("order", ~factor(., levels = rmse.order)) %>%
  arrange(order) %>%
  select(-order)

itab.1.c2 <- map(tables_c2, ~.x %>%
                select(sigma_estimate, estimator, txfx_2014 = p14)) %>%
  map2(c("2009-2011", "2010-2012", "2011-2013"), ~mutate(.x, covariates = .y)) %>%
  invoke(rbind, .) %>%
  spread(covariates, txfx_2014) %>%
  arrange(rev(estimator)) %>%
  mutate_at("sigma_estimate", ~stringr::str_replace_all(.,
                                                        c("sigma_uu_avg" = "Homogeneous",
                                                          "sigma_uu_i"   = "Heterogeneous",
                                                          "sigma_zero"   = "Unadjusted"))) %>%
  rename(Adjustment = sigma_estimate, Estimator = estimator) %>%
  mutate(order = paste(Adjustment, Estimator)) %>%
  mutate_at("order", ~factor(., levels = rmse.order.c2)) %>%
  arrange(order) %>%
  select(-order)

itab.2 <- tibble(
  Adjustment = rep(c("Heterogeneous", "Homogeneous", "Unadjusted"), 4),
  Estimator = rep(c("SBW", "H-SBW", "BC-SBW", "BC-HSBW"), each = 3),
  `J2013*gamma` = c(sbw, hsbw, bc.sbw, bc.hsbw),
  Target = truth[["2013"]]
) %>%
  mutate(order = paste(Adjustment, Estimator)) %>%
  mutate_at("order", ~factor(., levels = rmse.order)) %>%
  arrange(order) %>%
  select(-order)

print(xtable(itab.1,caption = "J_2014*gamma for all years"), 
      type="latex", caption.placement = "top",
      include.rownames = FALSE)

print(xtable(itab.2 ,caption = "J_2013*gamma_2014"), 
      type="latex", caption.placement = "top",
      include.rownames = FALSE)

# QA check: verify weights equal to primary results
c1_results <- readRDS("../04_Output/c1-results.rds") 
c2_results <- readRDS("../04_Output/c2-results.rds") 

map2_lgl(c1_results$Preferred$None$SBW, weights3$SBW, ~identical(.x$weights, .y$weights))
map2_lgl(c1_results$Preferred$None$`H-SBW`, weights3$HSBW, ~identical(.x$weights, .y$weights))
map2_lgl(c1_results$Preferred$None$`BC-HSBW`, weights3$`BC-HSBW`, ~identical(.x$weights, .y$weights))
map2_lgl(c1_results$Preferred$None$`BC-SBW`, weights3$`BC-SBW`, ~identical(.x$weights, .y$weights))

map2_lgl(c2_results$`Early expansion`$None$SBW, weights3.c2$SBW, ~identical(.x$weights, .y$weights))
map2_lgl(c2_results$`Early expansion`$None$`H-SBW`, weights3.c2$HSBW, ~identical(.x$weights, .y$weights))
map2_lgl(c2_results$`Early expansion`$None$`BC-HSBW`, weights3.c2$`BC-HSBW`, ~identical(.x$weights, .y$weights))
map2_lgl(c2_results$`Early expansion`$None$`BC-SBW`, weights3.c2$`BC-SBW`, ~identical(.x$weights, .y$weights))
