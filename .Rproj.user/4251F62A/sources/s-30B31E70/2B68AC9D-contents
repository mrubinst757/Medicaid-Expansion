# program: 06-etc-analysis.R
# purpose: run etc-analysis programs
# author: max rubinstein
# date modified: december 14, 2020

source("03_Analysis/05_etc-analysis.R")

# read data and estimated models ------------------------------------------------------------------
imputed_dat_c1 <- readRDS("../01_ProcessedData/calibrated-data-all.rds")

all_tdat_c1 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c1-all.rds") 

all_tdat_c2 <- readRDS("../01_ProcessedData/calibrated-data-all-subsets-c2-all-v2.rds") 

c1_results <- readRDS("../04_Output/c1-results.rds") 

c2_results <- readRDS("../04_Output/c2-results.rds") 

c1_jackknife <- readRDS("../04_Output/etc-jackknife-c1.rds")

c2_jackknife <- readRDS("../04_Output/etc-jackknife-c2.rds")

# set names ---------------------------------------------------------------------------------------
covariate_group <- c('None', 'Republican', "Unins & Unemp", "Urb-Age-Educ-Cit-Mar-Stu-Dis-F",
                     "Race-Eth-For-Inc-Pov", "Child-PGrowth-HHRatio")

sigma_estimator <- c("sigma_uu_avg", "sigma_uu_i", "sigma_zero")

estimator <- c("H-SBW", "BC-HSBW", "SBW", "BC-SBW")

adjust_replace <- c("Homogeneous", "Heterogeneous", "None")
names(adjust_replace) <- sigma_estimator

variables <- read_csv("../02_Specs/tol_specs.csv") %>%
  filter(`Reference Variable` == 0) %>%
  arrange(Variable) %>%
  .$Variable

#variables <- variables[-grep("total", variables)]
model_formula <- paste0("hins_unins_pct_2014 ~", paste0(variables, collapse = "+")) 

cdat <- imputed_dat_c1$data[[6]]
targets <- colMeans(cdat[variables])

c1_jackknife$data <- map(c1_jackknife$data, ~.x$data)
c2_jackknife$data <- map(c2_jackknife$data, ~.x$data)

all_tdat_c1 <- all_tdat_c1 %>%
  map(~.x[1:3])

all_tdat_c2 <- all_tdat_c2 %>%
  map(~.x[1:3])

# calculate y0estimate -----------------------------------------------------------------------------
m0 <- lm(as.formula(model_formula), data = cdat)
y0est <- c(1, targets) %*% coef(m0)
y0varest <- c(1, targets) %*% vcovCR(m0, type = "CR2", cluster = cdat$state) %*% c(1, targets)

# output results ----------------------------------------------------------------------------------------
all_txfx <- calc_all_txfx(all_tdat_c1, c1_results, all_tdat_c2, c2_results,
                          c1_jackknife, c2_jackknife, y0est, y0varest)

ccheck <- centering_check(all_txfx) 

tables_v1 <- paper_tables(all_txfx) 

output_all_plots(all_tdat_c1, all_tdat_c2, c1_results, c2_results, all_txfx)

tables_v1$c2_confint_filtered %>%
  separate(`CI (states)`, c("lci", "uci"), "\\,") %>%
  mutate_at(vars(contains("ci")), ~gsub("\\(|\\)", "", .)) %>%
  mutate_at(vars(contains("ci"), Psihat), ~as.numeric(.)) %>%
  filter(Adjustment != "Heterogeneous") %>%
  ggplot(aes(x = `Weight type`, y = Psihat, fill = Adjustment)) +
  geom_bar(stat = "identity", position = position_dodge2()) +
  geom_errorbar(position = position_dodge2(padding = 0.6), aes(ymin = lci, ymax = uci)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  xlab("Estimator") +
  ylab("Effect estimate") +
  ggsave("C:/Users/mdrub/Box/Medicaid-Project/02_Paper/01_Plots/point-estimates-c2.png")

print(xtable::xtable(tables_v1$c1_confint_filtered),
      booktabs = TRUE, include.rownames = FALSE, latex.environments = NULL)

print(xtable(tables_v1$c1_point_estimate), 
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

print(xtable::xtable(tables_v1$c2_confint_filtered),
      booktabs = TRUE, include.rownames = FALSE, latex.environments = NULL)

print(xtable(tables_v1$c2_point_estimate),
             digits = c(0, 0, 0, rep(2, 4)), 
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

c1_loo_table <- loo_table(all_txfx$loo_proc_c1, all_txfx$loo_states_c1, names(all_tdat_c1)[-1])
c2_loo_table <- loo_table(all_txfx$loo_proc_c2, all_txfx$loo_states_c2, names(all_tdat_c2)[-1])

loo_tab_c1 <- loo_table(all_txfx$c1_sensitivity_proc, all_txfx$c1_sensitivity_state, 
                        names(all_txfx$loo_states_c1)) %>%
  loo_summary("H-SBW", "sigma_uu_i_modeled")

loo_tab_c2 <- loo_table(all_txfx$c2_sensitivity_proc, all_txfx$c2_sensitivity_state, 
                        names(all_txfx$loo_states_c2)) %>%
  loo_summary("H-SBW", "sigma_uu_i_modeled")

print(xtable(loo_tab_c1, digits = rep(2, ncol(loo_tab_c1) + 1)),
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

print(xtable(loo_tab_c2, digits = rep(2, ncol(loo_tab_c2) + 1)),
      latex.environments = NULL, booktabs = TRUE, include.rownames = FALSE)

# republican governance differential analysis ------------------------------------------------
rdiff_c1_state_table <- calc_rdiff_se(all_txfx$txfx_c1, all_txfx$c1_sensitivity_state, 0.05)
rdiff_c1_proc_table <- calc_rdiff_se(all_txfx$txfx_c1, all_txfx$c1_sensitivity_proc, 0.05)
rdiff_c2_state_table <- calc_rdiff_se(all_txfx$txfx_c2, all_txfx$c2_sensitivity_state, 0.05)
rdiff_c2_proc_table <- calc_rdiff_se(all_txfx$txfx_c2, all_txfx$c2_sensitivity_proc, 0.05)

rdiff_c1_final <- rdiff_table(rdiff_c1_state_table, rdiff_c1_proc_table)
rdiff_c2_final <- rdiff_table(rdiff_c2_state_table, rdiff_c2_proc_table)

c1_rdiff <- rdiff_distr(c1_loo_table, rdiff_c1_final)
c2_rdiff <- rdiff_distr(c2_loo_table, rdiff_c2_final)

print(xtable::xtable(c1_rdiff), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(c2_rdiff), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

min(c1_rdiff$`0%`, c2_rdiff$`0%`)
max(c1_rdiff$`100%`, c2_rdiff$`100%`)

# produce heatmaps --------------------------------------------------------------------------
produce_heatmaps(all_txfx, "")

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

