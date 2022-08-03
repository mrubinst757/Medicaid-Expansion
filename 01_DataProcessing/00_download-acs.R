#program: 00_download-acs.R
#purpose: download acs microdata files
#author: max rubinstein
#date modified: october 16, 2020

#note: modify years to download data from other years

library(downloader)

setwd("C:/Users/mdrub/Box/Medicaid-Project/00_Data/")

output_dir <- "ACS/Raw_Microdata_2007-2016/"
years <- stringr::str_trunc(c("2009", "2010", "2011", "2012", "2013", "2014"), 2, side = "left", ellipsis = "") 
urls <- sprintf("https://www2.census.gov/programs-surveys/acs/data/pums/20%s/1-year/", years)

household.files      <- sprintf("csv_h%s", tolower(state.abb))
household.files.src  <- unlist(map(urls, ~paste0(.x, household.files, ".zip")))
household.files.dest <- unlist(map(c(2009:2014), ~paste0(household.files, .x, ".zip")))

person.files        <- sprintf("csv_p%s", tolower(state.abb))
person.files.src    <- unlist(map(urls, ~paste0(.x, person.files, ".zip")))
person.files.dest   <- unlist(map(c(2009:2014), ~paste0(person.files, .x, ".zip")))

map2(household.files.src, household.files.dest, ~downloader::download(.x, .y))
map2(person.files.src, person.files.dest, ~download(.x, .y))

map(household.files.dest, ~unzip(.x, exdir = output_dir))
map(person.files.dest, ~unzip(.x, exdir = output_dir))
