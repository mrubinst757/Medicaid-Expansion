#program: 03_process-urban.R
#purpose: calculate proportion urban for each cpuma
#author: max rubinstein
#date modifed: october 18, 2020

# load libraries ------------------------------------------------------------------------------------------------
library(tidyverse)
setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")

filename <- 'DEC_10_SF1_P2_with_ann.csv'

# process data --------------------------------------------------------------------------------------------------
data <- map(sprintf('02_DataDirectory/CensusTracts2010/state_%s/', 1:51), ~paste0(.x, filename)) %>%
  map(read_csv, skip = 1) %>%
  map(~set_names(.x, tolower(gsub(':', '', names(.x))))) %>%
  map(~select(.x, id, id2, geography, total, urban, rural))  %>%
  invoke(rbind, .) 

xwalk <- read.table('02_DataDirectory/CensusTracts2010/tractpumaxwalk2010.txt', sep = ',', header = FALSE) %>%
  set_names(c('state', 'county', 'tract', 'puma')) %>%
  slice(-1) %>%
  mutate(id2 = paste0(state, county, tract)) %>%
  mutate_if(is.factor, as.character)

cpumas <- read_csv('02_DataDirectory/PUMAS/CPUMA0010_PUMA2010_components.csv', col_types = 'ccccccc') %>%
  set_names(tolower(names(.))) %>%
  select(-state_fips) %>%
  mutate(state = cdlTools::fips(state_name, to = "FIPS"),
         state = gsub('(^[0-9]$)', '0\\1', state),
         state = if_else(is.na(state), '10', state))

# cpuma level urban rural classification (lose granularity but can use 2011-14)
final <- data %>% 
  left_join(xwalk, by = 'id2') %>%
  left_join(cpumas, by = c('puma', 'state')) %>%
  group_by(state, state_name, cpuma0010) %>%
  summarize(rural2010 = sum(rural), urban2010 = sum(urban)) %>%
  mutate(urban_prop = urban2010/(urban2010 + rural2010)) %>%
  select(-rural2010, -urban2010) %>%
  rename(cpuma = cpuma0010) %>%
  write_csv('02_DataDirectory/CensusTracts2010/urban_classifications.csv')
