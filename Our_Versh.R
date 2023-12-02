############################################################
rm(list=ls())

library(tidyverse)
library(cowplot)
library(tidyverse)
library(GetBCBData)
library(forecast)
library(knitr)
library(kableExtra)
library(writexl)
library(FinTS)
library(texreg)
library(rugarch)

###################################################
first_date <- '2003-12-01' 
last_date <- '2023-11-30' 

series_name <- 'Copper'
df = read_csv("copper_daily.csv")

df$Close = df$Close / 1000 # IN USD per Kilogram 

# select columns and calculated log_ret and arim_ret
df_prices <- df  %>%
  mutate(log_ret = log(Close/dplyr::lag(Close) ),   # Log Returns
         arim_ret = Close/dplyr::lag(Close) - 1,  # Arithmetic Returns
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

# get inflation data - using the United States CPI
df_inflation <- gbcbd_get_series(id = 3794, first.date = min(df_prices$Date), 
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

ggsave(filename = paste0('figs/fig02a_', series_name, '_prices.png'), 
       plot = p1)

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

ggsave(filename = paste0('figs/fig02b_', series_name, '_returns.png'), 
       plot = p2)


pb <- plot_grid(p1, p2, nrow = 2, 
               labels = 'AUTO')

pb


p <- ggAcf(x = df_prices$log_ret, lag.max = 10) +
  labs(title = paste0('Autocorrelogram for the Log Returns of ', series_name)) +
  theme_bw(base_family = "TT Times New Roman")

p

ggsave(filename = paste0('figs/fig02b_', series_name, '_ACF.png'), 
       plot = p)

##################################################################################

max_lag <- 5
my_html_file <- 'tabs/tab03-Arch_Test.html'
my_xlsx_file <- 'tabs/tab03-Arch_Test.xlsx'

if (!dir.exists(dirname(my_html_file))) dir.create(dirname(my_html_file))

tab_out <- do_arch_test(x = df_prices$log_ret, max_lag = max_lag)

tab_out

# remove attributes of table so it can be correctly parsed in html
tab_out <- as.data.frame(
  lapply(tab_out, function(x) { attributes(x) <- NULL; x })
)
str(tab_out)

rownames(tab_out) <- NULL

# save table in html
my_tbl <- knitr::kable(tab_out, format = 'html' ) %>%
  kable_styling(bootstrap_options = c("striped"), 
                full_width = FALSE ) 

my_tbl

cat(my_tbl, file = my_html_file)  

# write to excel
write_xlsx(x = tab_out, path = my_xlsx_file)
########################################################################

# OPTIONS
ar_lag <- 0 # lag used for ar term in mean equation (0 in paper)
ma_lag <- 0 # lag used for ma term in mean equation (0 in paper)
arch_lag <- 1 # lag in arch effect (1 in paper)
garch_lag <- 1 # lag in garch effect (1 in paper)
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch manual for more
distribution_to_estimate <- 'norm' # distribution used in all models
my_html_file <- 'tabs/tab04-estimation_garch.html' # where to save html file?

# close all opened windows
graphics.off()

df_grid <- expand_grid(ar_lag,
                       ma_lag,
                       arch_lag,
                       garch_lag,
                       models_to_estimate,
                       distribution_to_estimate)


estimate_garch <- function(ar_lag,
                           ma_lag,
                           arch_lag,
                           garch_lag,
                           models_to_estimate,
                           distribution_to_estimate) {
  
  message('Estimating ARMA(',ar_lag,',', ma_lag, ')', '-',
          models_to_estimate, '(', arch_lag, ',', garch_lag, ') ', 
          'dist = ', distribution_to_estimate)
  
  # estimate model
  my_spec <- ugarchspec(variance.model = list(model = models_to_estimate,
                                              garchOrder = c(arch_lag, 
                                                             garch_lag)),
                        mean.model = list(armaOrder = c(ar_lag,
                                                        ma_lag)), 
                        distribution.model = distribution_to_estimate)
  
  my_garch <- ugarchfit(spec = my_spec, data = df_prices$log_ret)
  
  return(my_garch)
}

# estimate all models
l_args <- as.list(df_grid)
l_models <- pmap(.l = l_args, .f = estimate_garch)

# make sure dir "tabs" exists
if (!dir.exists('tabs')) dir.create('tabs')

# reformat models for texreg
l_models <- map(l_models, extract.rugarch, include.rsquared = FALSE)

# write custom row
custom_row <- list('Variance Model' = df_grid$models_to_estimate,
                   'Distribution' = df_grid$distribution_to_estimate)
custom_names <- paste0('Model ', 1:length(l_models))

# save to html
htmlreg(l_models, 
        file = my_html_file, 
        custom.gof.rows = custom_row,
        custom.model.names = custom_names, 
        digits = 3)

# print to screen
screenreg(l_models,
          custom.gof.rows = custom_row,
          custom.model.names = custom_names, 
          digits = 3)
#########################################################

