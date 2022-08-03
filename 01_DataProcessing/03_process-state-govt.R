#program: 02_process-state-govt.R
#purpose: process state legislature and gov control data
#author: max rubinstein
#date modified: october 18, 2020

# load libraries and set dirs ------------------------------------------------
library(tidyverse)

setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")
input_dir  <- "ElectionData/"
output_dir <- "Intermediate_Files/"

# specify file names and parameters ------------------------------------------
file <- "NCSL Partisan Composition 1997-2018.xlsx" %>%
  paste0(input_dir, .)

skips = c(rep(7, 6), rep(0, 14), rep(1, 2))

new_names  <- c("State", "Tot_Seats", "Tot_Sen", "Sen_D",
               "Sen_R", "Sen_I", "Sen_O", "Tot_H",
               "H_D", "H_R", "H_I", "H_O", "Leg_Ctrl", "Gov")

new_names1 <- new_names[!new_names %in% c("H_I", "Sen_I")] 
new_names2 <- c(new_names1, "State_Ctrl")

full_names <- c(rep(list(new_names), 6), list(new_names1), rep(list(new_names2), 15))

# process data and save to disc --------------------------------------------------------
full_state_control <- map2(c(1:22), skips, ~readxl::read_excel(file, sheet = .x, skip = .y)) %>%
  map(~select_if(.x, ~sum(!is.na(.)) > 0)) %>% 
  map(~select(.x, -contains("X"), -matches("^...[0-9]?[0-9]$"))) %>% 
  map2(full_names, ~set_names(.x, .y)) %>%
  map(~mutate_at(.x, vars(-State, -contains("Ctrl"), -Gov), as.numeric)) %>%
  map(~mutate_at(.x, vars(-State, -contains("Ctrl"), -Gov), ~replace(.x, is.na(.), 0))) %>%
  map_at(1:6, ~mutate(.x, H_O = H_O + H_I) %>% select(-H_I)) %>%
  map_at(1:6, ~mutate(.x, Sen_O = Sen_O + Sen_I) %>% select(-Sen_I)) %>%
  map2(c(1997:2018), ~mutate(.x, Year = .y)) %>%
  map(~mutate_at(.x, "State", trimws)) %>%
  invoke(bind_rows, .) %>%
  filter(!is.na(State), !(State %in% c("Total", "Total States"))) %>%
  set_names(tolower(names(.))) %>%
  mutate(tot_ctrl = case_when(
    leg_ctrl == "Rep" & gov == "Rep" ~ "Rep",
    leg_ctrl == "Dem" & gov == "Dem" ~ "Dem",
    state == "Nebraska" ~ "Rep",
    state == "District of Columbia" ~ "Dem",
    TRUE ~ "Split"
  )) %>%
  select(-state_ctrl) %>%
  mutate(pct_repub_leg   = 100*(sen_r + h_r)/tot_seats,
         pct_repub_house = 100*(h_r/tot_h),
         pct_repub_senat = 100*(sen_r/tot_sen))  %>%
  write_csv(paste0(output_dir, "full_state_control.csv"))

create_state_analytic <- function(full_state_control, year) {
  full_state_control %>%
    filter(state %in% c(state.name, "District of Columbia")) %>%
    mutate(state = usdata::state2abbr(state)) %>%
    mutate_at("year", as.character) %>% 
    filter(year == !!year) %>%
    mutate(repub_gov = if_else(gov == "Rep", 1, 0)) %>%
    dplyr::select(matches("state"), contains("repub"), tot_ctrl) %>%
    mutate(repub_lower_control = if_else(pct_repub_house >= 50 & state != "NE" | state == "NE", 1, 0)) %>%
    fastDummies::dummy_cols(select_columns = c("tot_ctrl"))  %>%
    set_names(tolower(names(.))) %>%
    write_csv(paste0(output_dir, "state_control_analytic-", year, ".csv"))
}

t1 <- create_state_analytic(full_state_control, "2013")
t2 <- create_state_analytic(full_state_control, "2012")
t3 <- create_state_analytic(full_state_control, "2011")

