# program: 04_cpuma-analytic.R
# purpose: create wide dataset from long, aggregating some variables across years
# author: max rubinstein
# date modified: october 26, 2020

# load libraries ---------------------------------------------------------------------------------------
library(tidyverse)

# data processing functions ----------------------------------------------------------------------------

# create annual sample sizes and total sample sizes for averaged variables
process_sample_sizes <- function(person_file, household_file, year_upper, year_lower) {
  count_data <- person_file %>%
    readRDS() %>%
    select(state, cpuma, year, person_counts = sample) %>%
    left_join(
      household_file %>%
        readRDS() %>%
        select(state, cpuma, year, household_counts = sample), by = c("cpuma", "state", "year")
    )
  
  
  annual_sample_sizes = count_data %>%
    pivot_wider(id_cols = c("cpuma", "state"),
                names_from = "year",
                values_from = c("person_counts", "household_counts"))
  
  total_sample_sizes = count_data %>%
    select(cpuma, state, year, person_counts, household_counts) %>%
    filter(year <= year_upper & year >= year_lower) %>%
    select(-year) %>%
    group_by(cpuma, state) %>%
    summarize_all(sum)
  
  sample_size_data = annual_sample_sizes %>%
    left_join(total_sample_sizes, by = c("cpuma", "state"))
  
  sample_size_data
}

# aggregated fixed variables by averaging across input timeframe
aggregate_fixed_vars <- function(data, fixed_variables, year_upper, year_lower) {
  data %>%
    select(cpuma, state, year, treatment, all_of(fixed_variables)) %>%
    rename(total_households_avg = cpuma_tot_hh) %>%
    filter(year <= as.numeric(year_upper) & year >= as.numeric(year_lower)) %>%
    select(-year) %>%
    group_by(cpuma, state) %>%
    summarize_all(mean) 
}

#create time-varying variables from 2011-2014 (uninsurance and unemployment rates)
create_time_varying_vars <- function(data, time_varying_variables) {
  data %>%
    select(state, cpuma, year, all_of(time_varying_variables)) %>%
    rename(adult_pop = sample, labor_force_pop = labor_force) %>%
    mutate(hins_unins_pct = 100*hins_unins/adult_pop,
           unemployed_pct = 100*unemployed/labor_force_pop) %>%
    pivot_wider(id_cols = c("state", "cpuma"),
                names_from = "year",
                values_from = c("hins_unins", "unemployed", "adult_pop", "labor_force_pop",
                                "hins_unins_pct", "unemployed_pct"))
}

# combine all files and create percentages with correct denominator values
create_cpuma_analytic <- function(data, fixed_variables, time_varying_variables, year_upper, year_lower) {
  fixed_data <- aggregate_fixed_vars(data, fixed_variables, year_upper, year_lower)
  time_vdata <- create_time_varying_vars(data, time_varying_variables)
  sample_sizes <- process_sample_sizes(person_file, household_file, year_upper, year_lower)
  
  fixed_vars_nhh <- fixed_variables[-grep("hh|child", fixed_variables)]
  fixed_vars_hh  <- fixed_variables[grep("hh|child", fixed_variables)] %>%
    gsub("cpuma_tot_hh", "total_households_avg", .)
  
  years <- seq(year_lower, year_upper, 1)
  adult_pop_form <- paste0("(", paste0(sprintf("adult_pop_%s", years), collapse = "+"), ")/", length(years))
  pgrowth_years <- seq(year_lower, year_upper)[-1]
  pop_growth_form <- paste0("(", paste0(sprintf("pop_growth_%s", pgrowth_years), 
                                        collapse = "+"), ")/", length(pgrowth_years))
  
  final_data <- fixed_data %>%
    left_join(time_vdata, by = c("cpuma", "state")) %>%
    mutate(adult_pop_avg := !!rlang::parse_expr(adult_pop_form)) %>%
    mutate_at(vars(all_of(fixed_vars_nhh)), funs(pct = 100*./adult_pop_avg)) %>%
    mutate_at(vars(all_of(fixed_vars_hh)), funs(pct = 100*./total_households_avg)) %>%
    mutate(pop_growth_2012 = 100*adult_pop_2012/adult_pop_2011,
           pop_growth_2013 = 100*adult_pop_2013/adult_pop_2012,
           pop_growth_2014 = 100*adult_pop_2014/adult_pop_2013,
           pop_growth_2011 = 100*adult_pop_2011/adult_pop_2010,
           pop_growth_2010 = 100*adult_pop_2010/adult_pop_2009,
           pop_growth := !!rlang::parse_expr(pop_growth_form),
           urban = 100*urban,
           avg_adult_hh_ratio = 100*adult_pop_avg/total_households_avg) %>%
    left_join(sample_sizes, by = c("cpuma", "state")) %>%
    select(-urban_pct, urban_pct = urban)
  
  final_data
}

touchups <- function(analytic_file, politics) {
    pol_vars <- politics %>%
      select(state, repub_gov, repub_lower_control, repub_total_control = tot_ctrl_rep) %>%
      mutate_if(is_numeric, ~.*100) 
    
    analytic_file  %>%
      map(~select(.x, order(colnames(.))) %>%
          mutate(state = cdlTools::fips(state, to = "Abbreviation")) %>%
          select(-contains("repub")) %>%
          left_join(pol_vars, by = "state") %>%
          ungroup() %>%
          arrange(state, cpuma))
}

# read and process data --------------------------------------------------------------------------------
person_file <- "../01_ProcessedData/missing-acs-person-data.rds"
household_file <- "../01_ProcessedData/missing-acs-household-data.rds"
variable_specs <- read_csv("../02_Specs/rawvariables.csv")
data <- readRDS("../01_ProcessedData/all-data-2009-2014-long-r0-r80.rds")
politics13 <- read_csv("../../00_Data/Intermediate_Files/state_control_analytic-2013.csv")
politics11 <- read_csv("../../00_Data/Intermediate_Files/state_control_analytic-2011.csv")

time_varying_variables <- subset(variable_specs, `time-varying` == 1)$variables
fixed_variables <- subset(variable_specs, `time-varying` == 0)$variables

analytic_file_true <- map(data, ~create_cpuma_analytic(.x, fixed_variables = fixed_variables, 
                                                       time_varying_variables = time_varying_variables,
                                                       year_upper = 2013, 
                                                       year_lower = 2011))

analytic_file_test  <- map(data, ~create_cpuma_analytic(.x, fixed_variables = fixed_variables, 
                                                       time_varying_variables = time_varying_variables,
                                                       year_upper = 2012,
                                                       year_lower = 2010))

analytic_file_valid  <- map(data, ~create_cpuma_analytic(.x, fixed_variables = fixed_variables, 
                                                        time_varying_variables = time_varying_variables,
                                                        year_upper = 2011,
                                                        year_lower = 2009))

analytic_file_true  <- touchups(analytic_file_true, politics13)
analytic_file_test  <- touchups(analytic_file_test, politics11)
analytic_file_valid <- touchups(analytic_file_valid, politics11)

saveRDS(analytic_file_true,  "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-true.rds")
saveRDS(analytic_file_test,  "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-test.rds")
saveRDS(analytic_file_valid, "../01_ProcessedData/cpuma-analytic-file-2009-2014-r0-r80-valid.rds")
