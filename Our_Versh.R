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
write_rds(df_prices, df_out)


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
df_prices <- read_rds('data/GARCH-Data.rds')
series_name <- df_prices$series_name[1]

df_prices$Date <- as.Date(df_prices$Date, format = "%m/%d/%Y")

# get inflation data - using the Broad National Consumer Index 
df_inflation <- gbcbd_get_series(id = 433, first.date = min(df_prices$Date), 
                                 last.date = max(df_prices$Date)) %>% 
  mutate(inf_index  = cumprod(1+value/100))

total_cop_ret <- last(df_prices$Close)/first(df_prices$Close)-1
total_inflation <- last(df_inflation$inf_index)/first(df_inflation$inf_index) - 1 
n_years <- as.numeric(max(df_prices$Date) - min(df_prices$Date))/365
ret_cop_year = (1+total_cop_ret)^(1/n_years) - 1
ret_inflation_year = (1+total_inflation)^(1/n_years) - 1

real_ret_cop <- (1+total_cop_ret)/(1+total_inflation) - 1
real_ret_cop_year <- (1+ret_cop_year)/(1+ret_inflation_year) - 1

# create first plot
p1 <- ggplot(df_prices, aes(x = Date, y = Close)) + 
  geom_line() + 
  labs(title = paste0('Prices of ', series_name),
       subtitle = paste0('Total nominal arithmetic return equals to ', 
                         my_perc(total_cop_ret),
                         ' (', my_perc(ret_cop_year), ' per year)\n',
                         'Total real return, adjusted for inflation, equals to ',
                         my_perc(real_ret_cop), 
                         ' (', my_perc(real_ret_cop_year), ' per year)'),
       x = '',
       y = 'Index Value',
       caption = 'Data from Yahoo Finance') + 
  theme_bw(base_family = "TT Times New Roman") 
p1

largest_tab <- df_prices %>%
  top_n(abs(log_ret), n = n_largest)

# create second plot
p2 <- ggplot(df_prices, 
             aes(x = Date, y = log_ret)) + 
  geom_line() + 
  labs(title = paste0('Nominal Daily Log Returns of ', series_name),
       subtitle = paste0('Red circles represent the largest ', n_largest, 
                         ' absolute price variations in the sample'),
       x = '',
       y = 'Log Returns',
       caption = 'Data from Business Insider') + 
  theme_bw(base_family = "TT Times New Roman") +
  geom_point(data = largest_tab, aes(x = Date, y = log_ret), 
             size = 3, color = 'red'  ) +
  scale_y_continuous(labels = scales::percent) + 
  labs(size = 'Absolute Price Variation') # + 
scale_color_brewer(palette = 'BrBG')

p2

pb <- plot_grid(p1, p2, nrow = 2, 
               labels = 'AUTO')

pb

p <- ggAcf(x = df_prices$log_ret, lag.max = 10) +
  labs(title = paste0('Autocorrelogram for the Log Returns of ', series_name)) +
  theme_bw(base_family = "TT Times New Roman")

p
