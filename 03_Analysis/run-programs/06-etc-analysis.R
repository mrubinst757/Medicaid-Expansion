# program: 06-etc-analysis.R
# purpose: run etc-analysis programs
# author: max rubinstein
# date modified: december 14, 2020

source("03_Analysis/05_etc-analysis.R")
source("03_Analysis/02_model-estimation.R")
library(gt)

# read data and estimated models ------------------------------------------------------------------
imputed_dat_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.rds")

all_dat_c1 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c1-all.rds") 

all_dat_c2 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c2-all.rds") 

c1_results <- readRDS("../04_Output/c1-results.rds") 

c2_results <- readRDS("../04_Output/c2-results.rds") 

c1_jackknife <- readRDS("../04_Output/etc-jackknife-c1.rds")

c2_jackknife <- readRDS("../04_Output/etc-jackknife-c2.rds")

# set names ---------------------------------------------------------------------------------------
covariate_group <- c("None")

sigma_estimator <- c("sigma_uu_avg", "sigma_uu_i", "sigma_zero")

estimator <- c("H-SBW", "BC-HSBW", "SBW", "BC-SBW")

adjust_replace <- c("Homogeneous", "Heterogeneous", "Unadjusted")
names(adjust_replace) <- sigma_estimator

variables <- read_csv("../02_Specs/tol_specs.csv") %>%
  filter(`Reference Variable` == 0) %>%
  arrange(Variable) %>%
  .$Variable

#variables <- variables[-grep("total", variables)]
model_formula <- paste0("hins_unins_pct_2014 ~", paste0(variables, collapse = "+")) 

cdat <- imputed_dat_c1$data[[6]]
tdat <- imputed_dat_c1$data[[3]]

targets <- colMeans(cdat[variables])
targets_ett <- colMeans(tdat[variables])

c1_jackknife$data <- map(c1_jackknife$data, ~.x$data)
c2_jackknife$data <- map(c2_jackknife$data, ~.x$data)

all_tdat_c1 <- all_dat_c1 %>%
  map(~.x[1:3])

all_tdat_c2 <- all_dat_c2 %>%
  map(~.x[1:3])

all_cdat_c1 <- all_dat_c1 %>%
  map(~.x[4:6]) %>%
  .$Preferred %>%
  remove_states() 

# calculate y0estimate -----------------------------------------------------------------------------
m0 <- lm(as.formula(model_formula), data = cdat)
y0est <- c(1, targets) %*% coef(m0)
y0varest <- c(1, targets) %*% vcovCR(m0, type = "CR2", cluster = cdat$state) %*% c(1, targets)

m1 <- lm(as.formula(model_formula), data = tdat)
y1est <- c(1, targets_ett) %*% coef(m1)
y1varest <- c(1, targets_ett) %*% vcovCR(m1, type = "CR2", cluster = tdat$state) %*% c(1, targets)

m1.c2 <- lm(as.formula(model_formula), data = filter(tdat, !state %in% c("CA", "CT", "MN", "NJ", "WA")))
y1est.c2 <- c(1, targets_ett_c2) %*% coef(m1.c2)
y1varest.c2 <- c(1, targets_ett_c2) %*% vcovCR(m1, type = "CR2", cluster = tdat$state) %*% c(1, targets_ett_c2)

dat1_c1 <- rbind(cdat, tdat) %>%
  mutate(diff = hins_unins_pct_2014 - hins_unins_pct_2013)

dat1_c2 <- filter(dat1_c1, !state %in% c("CA", "CT", "MN", "NJ", "WA"))

m1.did <- lm(diff ~ treatment, dat1_c1)
m2.did <- lm(diff ~ treatment, dat1_c2)

coef(m1.did)[2]; coef(m2.did)[2]

confint(lmtest::coeftest(m1.did, vcovCR(m1.did, type = "CR2", cluster = dat1_c1$state),
                         df = length(unique(dat1_c1$state)) - 1))[2,]

confint(lmtest::coeftest(m2.did, vcovCR(m2.did, type = "CR2", cluster = dat1_c2$state),
                         df = length(unique(dat1_c2$state)) - 1))[2,]

# output results ----------------------------------------------------------------------------------------
all_txfx <- calc_all_txfx(all_tdat_c1, c1_results, all_tdat_c2, c2_results,
                          c1_jackknife, c2_jackknife, y0est, y0varest)

ccheck <- centering_check(all_txfx) 

centered <- map(ccheck, ~.x %>% filter(variable_subset == "None"))

tables_v1 <- paper_tables(all_txfx)

produce_heatmaps(all_txfx, "")

# Leave one out state estimates
loo_c1 <- map(all_txfx$loo_proc_c1, ~filter(.x, variable_subset == "None")) %>%
  map2(names(all_txfx$loo_proc_c1), ~mutate(.x, `Left out state` = .y)) %>%
  invoke(rbind, .) %>%
  bind_rows(
    all_txfx$txfx_c1 %>%
      filter(variable_subset == "None") %>%
      mutate(`Left out state` = "All")
  ) %>%
  spread(weight_type, psihat) %>%
  select(`Left out state`, Adjustment = sigma_estimate, everything(), -variable_subset) %>%
  nest(-Adjustment, -`Left out state`) %>%
  spread(Adjustment, data) %>%
  unnest()

loo_c2 <- map(all_txfx$loo_proc_c2, ~filter(.x, variable_subset == "None")) %>%
  map2(names(all_txfx$loo_proc_c2), ~mutate(.x, `Left out state` = .y)) %>%
  invoke(rbind, .) %>%
  bind_rows(
    all_txfx$txfx_c2 %>%
      filter(variable_subset == "None") %>%
      mutate(`Left out state` = "All")
  ) %>%
  spread(weight_type, psihat) %>%
  select(`Left out state`, Adjustment = sigma_estimate, everything(), -variable_subset) %>%
  nest(-Adjustment, -`Left out state`) %>%
  spread(Adjustment, data) %>%
  unnest()

loo_c1_table <- loo_c1 %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  gt() %>%
  tab_spanner("Homogeneous", c("BC-HSBW", "BC-SBW", "H-SBW", "SBW")) %>%
  tab_spanner("Heterogeneous", c("BC-HSBW1", "BC-SBW1", "H-SBW1", "SBW1")) %>%
  tab_spanner("Unadjusted", c("BC-HSBW2", "BC-SBW2", "H-SBW2", "SBW2")) %>%
  cols_label(`BC-HSBW1` = "BC-HSBW", `BC-SBW1` = "BC-SBW", `H-SBW1` = "H-SBW",
             `SBW1` = "SBW", `BC-HSBW2` = "BC-HSBW", `H-SBW2` = "H-SBW",
             `BC-SBW2` = "BC-SBW", `SBW2` = "SBW") %>%
  as_latex()

loo_c2_table <- loo_c2 %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  gt() %>%
  tab_spanner("Homogeneous", c("BC-HSBW", "BC-SBW", "H-SBW", "SBW")) %>%
  tab_spanner("Heterogeneous", c("BC-HSBW1", "BC-SBW1", "H-SBW1", "SBW1")) %>%
  tab_spanner("Unadjusted", c("BC-HSBW2", "BC-SBW2", "H-SBW2", "SBW2")) %>%
  cols_label(`BC-HSBW1` = "BC-HSBW", `BC-SBW1` = "BC-SBW", `H-SBW1` = "H-SBW",
             `SBW1` = "SBW", `BC-HSBW2` = "BC-HSBW", `H-SBW2` = "H-SBW",
             `BC-SBW2` = "BC-SBW", `SBW2` = "SBW") %>%
  as_latex()

loo_c1_table %>% cat() %>% print()
loo_c2_table %>% cat() %>% print()

main_results_paper <- GenMainTable(tables_v1, "Heterogeneous")
main_results_appendix <- GenMainTable(tables_v1, "")

print(xtable::xtable(main_results_paper), include.rownames = FALSE,
      latex.environments = NULL, booktabs = TRUE)

print(xtable::xtable(main_results_appendix), include.rownames = FALSE,
      latex.environments = NULL, booktabs = TRUE)

tables_v1$c1_confint_filtered <- tables_v1$c1_confint_filtered %>%
  select(-`CI (states)`, `95% CI` = `CI (proc)`)

tables_v1$c2_confint_filtered <- tables_v1$c2_confint_filtered %>%
  select(-`CI (states)`, `95% CI` = `CI (proc)`)

print(xtable::xtable(tables_v1$c1_confint_filtered),
      booktabs = TRUE, include.rownames = FALSE, latex.environments = NULL)

print(xtable(tables_v1$c1_point_estimate), 
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

print(xtable::xtable(tables_v1$c2_confint_filtered),
      booktabs = TRUE, include.rownames = FALSE, latex.environments = NULL)

print(xtable(tables_v1$c2_point_estimate),
      digits = c(0, 0, 0, rep(2, 4)), 
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

# calculate pre-treatment trends -----------------------------------------------
py <- c("hins_unins_pct_2009", "hins_unins_pct_2010",
        "hins_unins_pct_2011", "hins_unins_pct_2012",
        "hins_unins_pct_2013", "hins_unins_pct_2014")

truth <- colMeans(cdat[, py])
names(truth) <- c("2009", "2010", "2011", "2012", "2013", "2014")
txtruth <- colMeans(tdat[, py])
tdat.c2 <- filter(tdat, !state %in% c("CA", "CT", "MN", "NJ", "WA"))
txtruth.c2 <- colMeans(tdat.c2[, py])

trends_table <- round(rbind(truth, txtruth, txtruth.c2), 2) %>%
  as_tibble() %>%
  mutate(Treatment = c("Non-expansion", "Expansion (primary dataset)", "Expansion (early excluded)")) %>%
  select(Treatment, everything())

print(xtable::xtable(trends_table, caption = "Mean non-elderly adult uninsurance rates, 2009-2014"), 
      type="latex", caption.placement = "top",
      include.rownames = FALSE)

# check convergence --------------------------------------------------------------------------
check_status <- function(result_list, weight_type, sigma_type) {
  map(result_list, ~map(.x, ~.x[[weight_type]][[sigma_type]]$status)) %>%
    map(unlist) %>%
    unlist() 
}

all(check_status(c1_results, "SBW", "sigma_uu_i_modeled") == "solved")
all(check_status(c1_results, "SBW", "sigma_uu_avg") == "solved")
all(check_status(c1_results, "SBW", "sigma_zero") == "solved")
all(check_status(c1_results, "H-SBW", "sigma_uu_i_modeled") == "solved")
all(check_status(c1_results, "H-SBW", "sigma_uu_avg") == "solved")
all(check_status(c1_results, "H-SBW", "sigma_zero") == "solved")

centered$c1_proc; centered$c2_proc
