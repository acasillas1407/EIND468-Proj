############################################################
rm(list=ls())

library(tidyverse)
library(cowplot)
library(tidyverse)
library(GetBCBData)
library(forecast)

###################################################
first_date <- '2003-12-01' 
last_date <- '2023-11-30' 

series_name <- 'Copper'
df = read_csv("copper_daily.csv")

df$Close = df$Close / 1000 # IN USD per Kilogram 

# select columns and calculated log_ret and arim_ret
df_prices <- df  %>%
  mutate(log_ret = log(Close/dplyr::lag(Close) ),
         arim_ret = Close/dplyr::lag(Close) - 1,
         series_name = series_name) %>%
  na.omit() # remove all NA values

# save data into file
df_out <- 'GARCH-Data.rds'
write_rds(df_prices, rds_out)


#################################################
# OPTIONS
n_largest <- 10 # number of largest absolute returns to plot

# END OPTIONS

# close all existing plot windows
graphics.off()

# make sure folder "fig" exists
if (!dir.exists('figs')) dir.create('figs')

# source functions 
source('garch_fcts.R')

# get price data
df_prices <- read_rds('data/RAC-GARCH-Data.rds')
series_name <- df_prices$series_name[1]

# get inflation data - using the Broad National Consumer Index 
df_inflation <- gbcbd_get_series(id = 433, first.date = min(df_prices$Date), 
                                 last.date = max(df_prices$Date)) %>%
  mutate(inf_index  = cumprod(1+value/100))



