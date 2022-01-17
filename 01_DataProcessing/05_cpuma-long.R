# program: 04_cpuma-long.R
# purpose: create cpuma-year-state level dataset
# author: max rubinstein
# date modified: october 26, 2020

# load libraries and set dirs ------------------------------------------------------------------------------------------
setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")

library(tidyverse)
library(fastDummies)
library(furrr)
future::plan(multiprocess)

years <- 2009:2014 
file_end <- paste0(years[1], '-', years[length(years)])

input_dir           <- "ACS/Processed_Microdata/"
input_dir_pumas     <- "PUMAS/"
input_dir_politics  <- "Intermediate_Files/"
output_dir          <- "03_RDirectory/01_Processed_Data/"

household_file_paths = list.files(input_dir, pattern = 'household') %>%
  paste0(input_dir, .)

person_file_paths = list.files(input_dir, pattern = 'person') %>%
  paste0(input_dir, .)

#####################################################################################################
######################################### read data #################################################
#####################################################################################################

# puma/cpuma crosswalks
pumas00  <- read_csv(paste0(input_dir_pumas, 'CPUMA0010_PUMA2000_assignments.csv')) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

pumas10  <- read_csv(paste0(input_dir_pumas, 'CPUMA0010_PUMA2010_components.csv')) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

# medicaid expansion data 
mdcaid <- read_csv("Medicaid/kaestnerclassification.csv") %>%
  mutate(expansion_type = case_when(
    expansion == 0 ~ 'Non-expansion',
    expansion == 1 ~ 'Prior limited expansion',
    expansion == 2 ~ 'Prior full expansion',
    expansion == 3 ~ 'Large expansion',
    expansion == 4 ~ 'Moderate expansion'
  ),
  treatment = case_when(
    expansion %in% c(0, 1) ~ 0,
    expansion %in% c(3, 4) ~ 1,
    expansion == 2 ~ 2
  )) %>%
  filter(statename != "District of Columbia") %>%
  rename(exp.descr = expansion, state = statename) %>%
  mutate_at("state", ~as.character(cdlTools::fips(., to = "FIPS"))) %>%
  mutate_at("state", ~if_else(stringr::str_length(.) == 1, paste0("0", .), .))

# state politics data 
state_control <- read_csv(paste0(input_dir_politics, "state_control_analytic-2013.csv")) %>%
  mutate_at("state", ~as.character(cdlTools::fips(., to = "FIPS"))) %>%
  mutate_at("state", ~if_else(stringr::str_length(.) == 1, paste0("0", .), .))

# acs urban classifications
urban <- read_csv("CensusTracts2010/urban_classifications.csv") %>%
  rename(urban = urban_prop)

##########################################################################################################
###################################### aggregate acs data ################################################ 
##########################################################################################################

# aggregate data to state-cpuma-year level
aggregate_household_data <- function(household_microdata_files, weights) {
  
  household_data = household_microdata_files %>%
    read_rds() %>%
    rename(weights = !!weights) %>%
    group_by(state, cpuma, year) %>%
    mutate(cpuma_tot_hh = 1) %>%
    summarize(across(matches('child|tot_hh'), ~sum(.x*weights, na.rm = T))) 
  
  household_data
}

aggregate_person_data <- function(person_microdata_file, weights) {
  
  var_string = 'age_cat|race|inc_pov|educ|hispan|citizens|foreig|married|stud|ility|female|hins|unemp|force|sample'
  person_data = person_microdata_file %>%
    read_rds() %>%
    select(-race, -age_cat, -inc_pov, -educ, -inc_pov2, -age_cat2) %>%
    rename(pweights = !!weights) %>%
    mutate(pweights = as.numeric(pweights)) %>%
    mutate(across(matches(var_string), ~as.numeric(.x))) %>%
    group_by(state, cpuma, year) %>%
    summarize(across(matches(var_string), ~sum(.x*pweights, na.rm = T)))
  
  person_data
}

gen_household_replicates <- function(household_data_file, weight_list) {
  map(weight_list, ~aggregate_household_data(household_data_file, .x))
}

gen_person_replicates <- function(person_data_files, weight_list) {
  map(weight_list, ~aggregate_person_data(person_data_files, .x))
}

weight_list <- sprintf("wgtp%s", c('', 1:80))

household_replicates <- map(household_file_paths, ~gen_household_replicates(.x, weight_list)) %>%
  transpose() %>%
  map(~invoke(rbind, .x))

person_replicates <- map(person_file_paths, ~gen_person_replicates(.x, paste0('p', weight_list))) %>%
  transpose() %>%
  map(~invoke(rbind, .x))

all_replicates <- map2(person_replicates, household_replicates, 
                      ~left_join(.x, .y, by = c("cpuma", "state", "year"))) %>%
  map(~left_join(.x, urban, by = c("cpuma", "state"))) %>%
  map(~left_join(.x, state_control, by = c("state"))) %>%
  map(~left_join(.x, mdcaid, by = c("state")))

setwd("C:/Users/mdrub/Box/Medicaid-Project/01_Programs")
saveRDS(all_replicates, "01_ProcessedData/all-data-2009-2014-long-r0-r80.rds")
