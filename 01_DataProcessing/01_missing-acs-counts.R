# program: 06_missing-acs-counts.R
# purpose: separate program to generate counts of missing acs microdata
# author: max rubinstein
# date modified: october 26, 2020

##################################################################################################
########################## setting directories & preparing files $$###############################
##################################################################################################


# load libraries and set dirs ------------------------------------------------------------------------------------------
library(tidyverse)
library(assertthat)

setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")
input_dir  <- 'ACS/Raw_Microdata_2007-2016/'
input_dir_pumas <- 'PUMAS/'
output_dir <- 'ACS/Processed_Microdata/'

# read data and prepare file names ---------------------------------------------------------------------------

# puma/cpuma crosswalks 
pumas00  <- read_csv(paste0(input_dir_pumas, 'CPUMA0010_PUMA2000_assignments.csv')) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

pumas10  <- read_csv(paste0(input_dir_pumas, 'CPUMA0010_PUMA2010_components.csv')) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

# arrange file names
gen_state_list <- function(files, state_num) {
  c(unlist(map(1:length(files), ~files[[.x]][state_num])))
}
years = 2009:2014
chara_years = stringr::str_trunc(years, 2, side = 'left', ellipsis = '')

# acs filenames
person_file_patterns <- map(chara_years, ~sprintf('ss%sp[a-z][a-z].csv', .x))
person_files         <- map(person_file_patterns, ~list.files(input_dir, pattern = .x))

household_file_patterns <- map(chara_years, ~sprintf('ss%sh[a-z][a-z].csv', .x))
household_files         <- map(household_file_patterns, ~list.files(input_dir, pattern = .x))

person_file_list <- map(1:50, ~gen_state_list(person_files, .x)) %>%
  map(~paste0(input_dir, .x)) %>%
  unlist()

household_file_list <- map(1:50, ~gen_state_list(household_files, .x)) %>%
  map(~paste0(input_dir, .x)) %>%
  unlist()

##################################################################################################
########################## count missing values in acs microdata #################################
##################################################################################################

count_missing_person_data <- function(person_file_name) {
  char_year = unique(stringr::str_extract(person_file_name, 'ss([0-9][0-9])'))
  num_year  = as.numeric(gsub('ss', '', paste0('20', char_year)))
  
  cpuma_data <- (if (num_year < 2012) pumas00 else pumas10)
  
  missing = person_file_name %>%
    read_csv(col_types = cols(
      ST = col_character(),
      PUMA  = col_character())) %>%
    set_names(tolower(names(.))) %>%
    select(state = st, puma, serialno, cit, agep, starts_with('hins'), hicov, mar, schl, schg, sex,  
           msp, nativity, povpip, rac1p, esr, hispanic = hisp, disability = dis) %>%
    mutate(missing_age_all = sum(is.na(agep))/n()) %>%
    mutate(sample = if_else(between(agep, 19, 64), 1, 0)) %>%
    filter(sample == 1) %>%
    mutate(missing_citizenship =  if_else(is.na(cit), 1, 0),
           missing_foreign_born = if_else(is.na(nativity), 1, 0),
           missing_employment = if_else(is.na(esr), 1, 0),
           missing_age_sample = if_else(is.na(agep), 1, 0),
           missing_disability = if_else(is.na(disability), 1 ,0),
           missing_hispanic = if_else(is.na(hispanic), 1, 0),
           missing_hins1 = if_else(is.na(hins1), 1, 0),
           missing_hins2 = if_else(is.na(hins2), 1, 0),
           missing_hins3 = if_else(is.na(hins3), 1, 0),
           missing_hins4 = if_else(is.na(hins4), 1, 0),
           missing_hins5 = if_else(is.na(hins5), 1, 0),
           missing_hins6 = if_else(is.na(hins6), 1, 0),
           missing_hins7 = if_else(is.na(hins7), 1, 0),
           missing_educ = if_else(is.na(schl), 1, 0),
           missing_student = if_else(is.na(schg), 1, 0),
           missing_race = if_else(is.na(rac1p), 1, 0),
           missing_inc_pov = if_else(is.na(povpip), 1, 0),
           missing_sex = if_else(is.na(sex), 1, 0),
           missing_married = if_else(is.na(mar), 1, 0),
           missing_student = if_else(is.na(schg), 1, 0)) %>%
    select(state, puma, sample, contains('missing')) %>%
    left_join(cpuma_data, by = c('state', 'puma')) %>%
    select(-puma) %>%
    mutate(year = num_year) %>%
    group_by(state, cpuma, year) %>%
    summarize_all(sum)
  
  missing
}

count_missing_household_data <- function(household_file_name) {
  char_year = unique(stringr::str_extract(household_file_name, 'ss([0-9][0-9])'))
  num_year  = as.numeric(gsub('ss', '', paste0('20', char_year)))
  
  cpuma_data <- (if (num_year < 2012) pumas00 else pumas10)
  
  household_sample <- household_file_name %>%
    read_csv(col_types = cols(
      ST = col_character(),
      PUMA  = col_character())) %>%
    set_names(tolower(names(.))) %>%
    select(state = st, puma, serialno, noc) %>%
    mutate(missing_noc = is.na(noc), year = num_year, sample = 1) %>%
    left_join(cpuma_data, by = c('state', 'puma')) %>%
    select(-puma) %>%
    group_by(state, year, cpuma) %>%
    summarize(missing_noc = sum(missing_noc), sample = sum(sample)) 
  
  household_sample
}

person_all_missing <- map(person_file_list, count_missing_person_data) %>%
  invoke(rbind, .)

household_all_missing <- map(household_file_list, count_missing_household_data) %>%
  invoke(rbind, .)

setwd("C:/Users/mdrub/Box/Medicaid-Project/01_Programs/01_ProcessedData/")

saveRDS(person_all_missing, "missing-acs-person-data.rds")
saveRDS(household_all_missing, "missing-acs-household-data.rds")
