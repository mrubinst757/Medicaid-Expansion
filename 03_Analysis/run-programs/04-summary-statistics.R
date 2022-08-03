# program: 04-summary-statistics.R
# purpose: run basic summary statistics on unadjusted/adjusted datasets
# author: max rubinstein
# date modified: january 29, 2021

# load libraries and read data --------------------------------------------------------------------------
library(tidyverse)
library(Hmisc)
library(corrplot)
library(RColorBrewer)

varnames <- read_csv("../02_Specs/codebook.csv") %>%
  filter(Pct == 1)
var_names <- varnames$`Variable Name`
names(var_names) <- varnames$Variable

c1_data_all <- readRDS("../01_ProcessedData/calibrated-data-all.rds") %>%
  unnest() %>%
  nest(-key, -set) %>%
  .$data %>%
  map(~arrange(.x, state, cpuma))

c1_jackknife <- readRDS("../04_Output/etc-jackknife-c1.rds")

c2_jackknife <- readRDS("../04_Output/etc-jackknife-c2.rds")

read_vars <- function(col_name) {
  read_csv("../02_Specs/tol_specs.csv") %>%
    mutate(`Reference Variable` = case_when(grepl("child", Variable) ~ 0, TRUE ~ `Reference Variable`)) %>%
    filter(`Reference Variable` == 0) %>%
    .[[col_name]] %>%
    sort()
}

variables  <- read_vars("Variable")
vars_test  <- read_vars("var_test")
vars_valid <- read_vars("var_valid")

# helper funs --------------------------------------------------------------------------------------------

# calculate range of a variable
my_range <- function(x) {
  range(x)[2] - range(x)[1]
}

# produce summary statistics table
produce_tables <- function(data_list, variables) {
  summary_table <- map(data_list, ~.x %>%
                         select(treatment, all_of(variables)) %>%
                         group_by(treatment) %>%
                         summarize_all(funs(mean, IQR, my_range)) %>%
                         gather(key, value, -treatment) %>%
                         mutate(stat = case_when(
                           grepl("mean", key) ~ "mean",
                           grepl("IQR", key) ~ "IQR",
                           grepl("range", key) ~ "range")) %>%
                         mutate_at("key", ~gsub("_mean|_IQR|_my_range", "", .)) %>%
                         spread(stat, value)) %>%
    map2(c("sigma_uu_i", "sigma_uu_avg", "sigma_zero"),
         ~mutate(.x, sigma_estimate = .y)) %>%
    invoke(rbind, .) %>%
    mutate_at("key", ~stringr::str_replace_all(., var_names))
  
  final_table <- summary_table %>%
    mutate(`(mean, sd)` = paste0("(", round(mean, 1), ", ", 
                                 #round(sd, 2), ", ",
                                 round(IQR, 1), ", ",
                                 round(range, 1), ")")) %>%
    select(-mean, -IQR, -range) %>%
    spread(sigma_estimate, `(mean, sd)`) %>%
    arrange(key) %>%
    select(Treatment = treatment, Variable = key,
           `No adjustment` = sigma_zero, 
           `Heterogeneous` = sigma_uu_i, 
           `Homogeneous` = sigma_uu_avg) 
  
  paper_summary_table <- map(data_list, ~.x %>%
                               select(treatment, contains("unins_pct"), contains("unemployed_pct")) %>%
                               select(-contains("2014"), treatment) %>%
                               group_by(treatment) %>%
                               summarize_all(funs(sd)) %>%
                               gather(variable, sd, -treatment)) %>% 
    map2(c("sigma_uu_i", "sigma_uu_avg", "sigma_zero"),
         ~mutate(.x, sigma_estimate = .y)) %>%
    invoke(rbind, .) 
  
  ptab1 <- paper_summary_table %>%
    spread(sigma_estimate, sd) %>%
    filter(treatment == 1) %>%
    select(Variable = variable, `No adjustment` = sigma_zero, `Preferred` = sigma_uu_i, `Secondary` = sigma_uu_avg) %>%
    mutate_at("Variable", ~stringr::str_replace_all(., var_names))
  
  # calculate extreme value table ---------------------------------------------------------------------------
  lb_all <- data_list[[3]] %>%
    select(treatment, all_of(variables)) %>%
    group_by(treatment) %>%
    summarize_all(min) %>%
    gather(key, lb, -treatment)
  
  ub_all <- data_list[[3]] %>%
    select(treatment, all_of(variables)) %>%
    group_by(treatment) %>%
    summarize_all(max) %>%
    gather(key, ub, -treatment)
  
  support_table <- map(data_list, ~.x %>%
                         select(treatment, all_of(variables)) %>%
                         gather(key, value, -treatment) %>%
                         left_join(lb_all, by = c("key", "treatment")) %>%
                         left_join(ub_all, by = c("key", "treatment")))
  
  support_table <- map(data_list, ~.x %>%
                         select(treatment, all_of(variables)) %>%
                         gather(key, value, -treatment) %>%
                         left_join(lb_all, by = c("key", "treatment")) %>%
                         left_join(ub_all, by = c("key", "treatment")) %>%
                         group_by(key, treatment) %>%
                         dplyr::summarize(mean_support  = sum(value < lb | value > ub)) %>%
                         filter(!grepl("repub|urban", key)) %>%
                         spread(treatment, mean_support)) %>%
    map2(c("sigma_uu_i", "sigma_uu_avg", "sigma_zero"), ~mutate(.x, sigma_estimator = .y)) %>%
    invoke(rbind, .) %>%
    mutate(`Counts (control, treatment)` = paste0("(", `0`, ", ", `1`, ")")) %>%
    select(-`0`, -`1`) %>%
    spread(sigma_estimator, `Counts (control, treatment)`) %>%
    select(Variables = key, `Heterogeneous adjustment` = sigma_uu_i, 
           `Homogeneous adjustment` = sigma_uu_avg) %>%
    mutate_at("Variables", ~stringr::str_replace_all(., var_names))
  
  ptab2 <- paper_summary_table %>%
    spread(sigma_estimate, sd) %>%
    filter(treatment == 0) %>%
    select(Variable = variable, `No adjustment` = sigma_zero, `Preferred` = sigma_uu_i, `Secondary` = sigma_uu_avg) %>%
    mutate_at("Variable", ~stringr::str_replace_all(., var_names))
  
  list(final_table = final_table, support_table = support_table, 
       ptab1 = ptab1, ptab2 = ptab2)
}

# create tables ---------------------------------------------------------------------------
tables <- produce_tables(c1_data_all[1:3], variables)
tables_valid <- produce_tables(c1_data_all[7:9], variables)
tables_test <- produce_tables(c1_data_all[4:6], variables)

jtab1 <- tables$ptab1 %>%
  filter(!grepl("2009|2010", Variable))  %>%
  mutate_at("Variable", ~stringr::str_replace_all(., var_names))

jtab2 <- tables$ptab2 %>%
  filter(!grepl("2009|2010", Variable)) %>%
  mutate_at("Variable", ~stringr::str_replace_all(., var_names))

print(xtable::xtable(jtab1, digits = 2), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(jtab2, digits = 2), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(tables$final_table, digits = 0), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(tables$support_table, digits = 0), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(tables$ptab1, digits = 2), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

print(xtable::xtable(tables$ptab2, digits = 2), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

# create correlation matrix ------------------------------------------------------------------------
c1_data_all <- map(c1_data_all, ~set_names(., stringr::str_replace_all(names(.), var_names)))

cor1 <- cor(c1_data_all[[3]][var_names])

corrplot(cor1, order = "original", type = "lower", tl.pos = "ld",
         tl.cex = 0.75, tl.col = "black", tl.srt = 0.001)

plot1 <- corrplot(cor1, order = "original", type = "lower", tl.pos = "ld",
                  tl.cex = 0.75, tl.col = "black", tl.srt = 0.001)

file_path <- "../../02_Paper/01_Plots/correlation-plot-c1-sigma-zero.png"
png(height = 1600, width = 1600, file = file_path, type = "cairo")
corrplot(cor3, order = "original", type = "lower", tl.pos = "ld",
         tl.cex = 2, tl.col = "black", tl.srt = 0.001, cl.cex = 2)
dev.off()

# variable-delta table for appendix -----------------------------------------------------
tol_list <- read_csv("../02_Specs/tol_specs.csv") 

tol_table <- tol_list %>%
  select(Variable, `Delta` = `Base Tol`) %>%
  mutate_at("Variable", ~stringr::str_replace_all(., var_names)) %>%
  arrange(Delta) 

print(xtable::xtable(tol_table, digits = 2), include.rownames = FALSE,
      latex.environments = NULL, 
      booktabs = TRUE)

# state-treatment assignment table -----------------------------------------------------
excluded <- state.abb[!state.abb %in% unique(c1_data_all[[1]]$state)]
primary.1  <- unique(c1_data_all[[1]]$state[c1_data_all[[1]]$`Expansion State` == 1])
primary.0  <- unique(c1_data_all[[1]]$state[c1_data_all[[1]]$`Expansion State` == 0])
early.1 <- names(c2_jackknife$data)

paste0(primary.0, collapse =  ", ")

state.table <- tibble(
  `Treated states (primary)` = paste0(primary.1, collapse =  ", "),
  `Control states (primary)` = paste0(primary.0, collapse =  ", "),
  `Treated states (early expansion excluded)` = paste0(early.1, collapse =  ", "),
  `Always excluded` = paste0(c(excluded, "DC"), collapse =  ", ")
) 

print(xtable::xtable(state.table), include.rownames = FALSE)
