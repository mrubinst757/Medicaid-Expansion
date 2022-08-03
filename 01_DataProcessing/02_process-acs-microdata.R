# program: 01_process-acs-microdata.R
# purpose: take raw data, select variables, basic transformation, stack files
# author: max rubinstein
# date modified: october 18, 2020

# load libraries and set dirs ------------------------------------------------------------------------------------------
library(tidyverse)
library(assertthat)

setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")
input_dir  <- "ACS/Raw_Microdata_2007-2016/"
input_dir_pumas  <- "PUMAS/"
output_dir <- "ACS/Processed_Microdata/"

# arrange file names
gen_state_list <- function(person_files, household_files, state_num) {
  assert_that(length(person_files) == length(household_files))
  c(unlist(map(1:length(person_files), ~person_files[[.x]][state_num])),
    unlist(map(1:length(household_files),  ~household_files[[.x]][state_num])))
}

years = 2009:2014
chara_years = stringr::str_trunc(years, 2, side = "left", ellipsis = "")

person_file_patterns    <- map(chara_years, ~sprintf("ss%sp[a-z][a-z].csv", .x))
household_file_patterns <- map(chara_years, ~sprintf("ss%sh[a-z][a-z].csv", .x))
person_files            <- map(person_file_patterns, ~list.files(input_dir, pattern = .x))
household_files         <- map(household_file_patterns, ~list.files(input_dir, pattern = .x))

state_file_list <- map(1:50, ~gen_state_list(person_files, household_files, .x)) %>%
  map(~paste0(input_dir, .x))

# puma/cpuma crosswalks --------------------------------------------------------------------------------------
pumas00  <- read_csv(paste0(input_dir_pumas, "CPUMA0010_PUMA2000_assignments.csv")) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

pumas10  <- read_csv(paste0(input_dir_pumas, "CPUMA0010_PUMA2010_components.csv")) %>%
  select(state = State_FIPS, puma = PUMA, cpuma = CPUMA0010) %>%
  distinct()

# data processing funs ---------------------------------------------------------------------------------------

# add variables to person-level dataset
process_acs_person <- function(data, year) {
  cpuma_data <- (if (year < 2012) pumas00 else pumas10)
  
  data %>%
    set_names(tolower(names(.))) %>%
    select(state = st, puma, serialno, cit, agep, starts_with("hins"), hicov, mar, schl, schg, sex,  
           msp, nativity, povpip, rac1p, esr, hispanic = hisp, disability = dis, contains("wgtp")) %>%
    dplyr::mutate_at(vars(povpip, agep, pwgtp, hispanic), funs(as.numeric)) %>%
    mutate(citizenship  = if_else(cit == 5, 0, 1),
           foreign_born = if_else(nativity == 2, 1, 0),
           married      = if_else(mar == 1, 1, 0),
           unemployed   = if_else(esr == 3, 1, 0),
           labor_force  = if_else(esr == 6, 0, 1), 
           sample       = if_else(between(agep, 19, 64), 1, 0),
           noneld       = if_else(agep <= 64, 1, 0),
           disability   = if_else(disability == 1, 1, 0),
           uninsured    = if_else(hins1 == 2 & hins2 == 2 & hins3 == 2 & hins4 == 2 & hins5 == 2 & hins6 == 2 & hins7 == 2, 1, 0),
           hispanic     = if_else(hispanic == 1, 0, 1)
           )  %>%
    filter(sample == 1) %>%
    dplyr::mutate_at(vars(contains("hins")), funs(if_else(. == 2, 0, 1))) %>%
    dplyr::mutate(race = case_when(
      rac1p == 1 ~ "white",
      rac1p == 2 ~ "black",
      rac1p %in% c(3:9) ~ "other"
    ),
    educ = case_when(
      schl %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15")  ~ "less_than_hs",
      schl %in% c("16", "17")  ~ "hs_degree",
      schl %in% c("18", "19", "20")  ~ "some_college",
      schl %in% c("21", "22", "23", "24")  ~ "college"
    ),
    student = if_else(is.na(schg), 0, 1),
    inc_pov = case_when(
      povpip < 100 ~ "100",
      between(povpip, 100, 199) ~ "200",
      between(povpip, 200, 299) ~ "300",
      between(povpip, 300, 399) ~ "400",
      between(povpip, 400, 499) ~ "500",
      povpip >= 500 ~ "500_plus",
      is.na(povpip) ~ "missing"
    ),
    age_cat = case_when(
      between(agep, 19, 24) ~ "19_24",
      between(agep, 25, 29) ~ "25_29",
      between(agep, 30, 34) ~ "30_34",
      between(agep, 35, 39) ~ "35_39",
      between(agep, 40, 44) ~ "40_44",
      between(agep, 45, 49) ~ "45_49",
      between(agep, 50, 54) ~ "50_54",
      between(agep, 55, 59) ~ "55_59",
      between(agep, 60, 64) ~ "60_64"
    ),
    coarse_race_white = if_else(rac1p == 1, 1, 0),
    coarse_age_cat    = if_else(between(agep, 35, 64), 1, 0),
    coarse_inc_pov    = if_else(povpip <= 138, 1, 0),
    coarse_educ       = if_else(schl %in% c("18", "19", "20", "21", "22", "23", "24"), 1, 0),
    age_cat2 = case_when(
      between(agep, 19, 29) ~ "19_29",
      between(agep, 30, 39) ~ "30_39",
      between(agep, 40, 49) ~ "40_49",
      between(agep, 50, 64) ~ "50_64"
    ),
    inc_pov2 = case_when(
      povpip <= 138 ~ "138",
      between(povpip, 139, 299) ~ "139_299",
      between(povpip, 300, 499) ~ "300_499",
      povpip >= 500 ~ "500_plus"
    ),
    hins_any   = if_else(hins1 == 1 | hins2 == 1 | hins3 == 1 | hins4 == 1 | hins5 == 1 | hins6 == 1 | hins7 == 1, 1, 0),
    hins_unins = if_else(hins_any == 1, 0, 1),
    hins_priv  = if_else((hins1 == 1 | hins2 == 1) & (hins3 == 0 & hins4 == 0 & hins5 == 0 & hins6 == 0 & hins7 == 0), 1, 0),
    hins_mdcd  = if_else(hins4 == 1, 1, 0),
    hins_othg  = if_else(hins_priv == 0 & hins_mdcd == 0 & uninsured == 0, 1, 0),
    female = if_else(sex == 2, 1, 0)
    ) %>%
    fastDummies::dummy_cols(select_columns = c("educ", "inc_pov", "age_cat", "race", "inc_pov2", "age_cat2")) %>%
    left_join(cpuma_data, by = c("state", "puma")) %>%
    dplyr::mutate(year = year)
}

process_acs_household <- function(data, year) {

  cpuma_data <- (if (year < 2012) pumas00 else pumas10)
  
  data %>%
    set_names(tolower(names(.))) %>%
    select(state = st, serialno, puma, noc, contains("wgtp")) %>%
    mutate_at(vars(noc, contains("wgtp")), funs(as.numeric)) %>%
    left_join(cpuma_data, by = c("state", "puma")) %>%
    mutate(year = year,
           one_child = if_else(noc == 1, 1, 0),
           two_child = if_else(noc == 2, 1, 0),
           three_plus_child = if_else(noc >= 3, 1, 0),
           any_children = if_else(noc >= 1, 1, 0),
           missing_children = if_else(is.na(noc), 1, 0))
}

# process data and save to disk ------------------------------------------------------------------------------
annual_process <- function(file_names) {
  chara_years = unique(stringr::str_extract(file_names, "ss([0-9][0-9])"))
  num_years = as.numeric(gsub("ss", "", paste0("20", chara_years)))
  
  person_files <- file_names[grep("[0-9]p[a-z]", file_names)]
  household_files <- file_names[grep("[0-9]h[a-z]", file_names)]
  
  output_file_str <- stringr::str_trunc(person_files, 6, ellipsis = "", side = "left")[1] %>%
    gsub("csv", "rds", .)
  
  person_data <- map(person_files, ~read_csv(.x, col_types = cols(
    ST = col_character(),
    PUMA  = col_character()
  ))) %>%
    map2(num_years, process_acs_person) 

  household_data <- map(household_files, ~read_csv(.x, col_types = cols(
    ST = col_character(),
    PUMA  = col_character()
  ))) %>%
    map2(num_years, process_acs_household) %>%
    invoke(rbind, .)
  
  person_data <- invoke(rbind, person_data)
  
  write_rds(person_data, paste0(output_dir, "person_acs_", output_file_str))
  write_rds(household_data, paste0(output_dir, "household_acs_", output_file_str))  
}

# iterate program through all states and save to disc --------------------------------------------
map(state_file_list, annual_process)
