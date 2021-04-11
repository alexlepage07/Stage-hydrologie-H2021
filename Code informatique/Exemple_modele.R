rm(list=ls())
library(tidyverse)
library(readr)
library(lubridate)
library(Kendall)
library(copula)
library(VineCopula)

# Importation des fonctions
setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/Stage-hydrologie-H2021/Code informatique")
source("Rainfall_functions.R")


# Importation des données ----
setwd("C:/Users/alexl/Google Drive/1. Université/2. Stage hydrologie/data")
data <- readxl::read_xlsx('Clearwater river.xlsx', sheet = 'T089R09W4')

data <- data %>% transmute(
   DATE = as.Date(Date),
   PRCP = `Precip. (mm)`
) %>% filter(month(DATE) %in% c(4,5,6))
data %>% summary()

data %>% filter(year(DATE) %in% (2009:2013)) %>%
   group_by(year=year(DATE)) %>% select(-DATE) %>% mutate(l=1:91) %>%
   pivot_wider(names_from=year, values_from=PRCP) %>%
   top_n(-15, l)%>% t() %>% xtable::xtable() %>% print(include.rownames=FALSE)

data2 <- preprocess_data(  # Appuyer la méthode de NA imputing avec la littérature
   data, 
   NA.neighbour_days = 2,
   C.separeting_days = 1,
   n_print=0)

data2 <- data2 %>% filter(year(start_date) %in% (2009:2013)) %>%
   mutate(year=year(start_date))

data2 %>% filter(prcp_agg >= 1, year %in% c(2010)) %>% slice(1:2) %>%
   select(start_date) %>% t() %>% xtable::xtable() %>% print(include.rownames=FALSE)

data.extremes <- calculate_excedences(data2, 1, calculate_W=T)
data.extremes %>% filter(year(start_date) %in% c(2013)) %>% slice(1:3)

