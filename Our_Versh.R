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
library(purrr)
library(ggtext)
library(dplyr)
library(ggplot2)

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

colnames(df_prices)[colnames(df_prices) == "Date"] <- "ref.date"
colnames(df_prices)[colnames(df_prices) == "Close"] <- "price.adjusted"
df_prices <- df_prices[nrow(df_prices):1, ]

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


df_prices$ref.date <- as.Date(df_prices$ref.date, format = "%m/%d/%Y")


# get inflation data - using the United States CPI http://www.bcb.gov.br/?sgs
df_inflation <- gbcbd_get_series(id = 3794, first.date = min(df_prices$ref.date), 
                                 last.date = max(df_prices$ref.date)) %>% 
  mutate(inf_index  = cumprod(1+value/100))

total_cop_ret <- last(df_prices$price.adjusted)/first(df_prices$price.adjusted)-1
total_inflation <- last(df_inflation$inf_index)/first(df_inflation$inf_index) - 1 
n_years <- as.numeric(max(df_prices$ref.date) - min(df_prices$ref.date))/365
ret_cop_year = (1+total_cop_ret)^(1/n_years) - 1
ret_inflation_year = (1+total_inflation)^(1/n_years) - 1

real_ret_cop <- (1+total_cop_ret)/(1+total_inflation) - 1
real_ret_cop_year <- (1+ret_cop_year)/(1+ret_inflation_year) - 1

# create first plot
p1 <- ggplot(df_prices, aes(x = ref.date, y = price.adjusted)) + 
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
       caption = 'Data from Yahoo Finance') 
p1

ggsave(filename = paste0('figs/fig02a_', series_name, '_prices.png'), 
       plot = p1)

largest_tab <- df_prices %>%
  top_n(abs(log_ret), n = n_largest)

# create second plot
p2 <- ggplot(df_prices, 
             aes(x = ref.date, y = log_ret)) + 
  geom_line() + 
  labs(title = paste0('Nominal Daily Log Returns of ', series_name),
       subtitle = paste0('Red circles represent the largest ', n_largest, 
                         ' absolute price variations in the sample'),
       x = '',
       y = 'Log Returns',
       caption = 'Data from Business Insider') + 
  geom_point(data = largest_tab, aes(x = ref.date, y = log_ret), 
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
  labs(title = paste0('Autocorrelogram for the Log Returns of ', series_name))

p

ggsave(filename = paste0('figs/fig02c_', series_name, '_ACF.png'), 
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

max_lag_AR <- 1 # used 1 in paper
max_lag_MA <- 1 # used 1 in paper
max_lag_ARCH <- 2 # used 2 in paper
max_lag_GARCH <- 1 # used 1 in paper
dist_to_use <- c('norm', 'std') # see rugarch::ugarchspec help for more

graphics.off()

out <- find_best_arch_model(x = df_prices$log_ret, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)

tab_out <- out$tab_out

# pivot table to long format (better for plotting)
df_long <- tidyr::pivot_longer(data = tab_out %>%
                                 select(model_name,
                                        type_model,
                                        type_dist,
                                        AIC, BIC),  cols = c('AIC', 'BIC'))

models_names <- unique(df_long$model_name)
best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
                 tab_out$model_name[which.min(tab_out$BIC)])

# figure out where is the best model
df_long <- df_long %>%
  mutate(order_model = if_else(model_name %in% best_models, 'Best Model', 'Not Best Model') ) %>%
  na.omit()

# make table with best models
df_best_models <- df_long %>%
  group_by(name) %>%
  summarise(model_name = model_name[which.min(value)],
            value = value[which.min(value)],
            type_model = type_model[which.min(value)])

# plot results
p3 <- ggplot(df_long %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")

 p3 ; ggsave('figs/fig04_best_garch.png')
 
 p3


# estimate best garch model by BIC (used in next section)
best_spec_bic = ugarchspec(variance.model = list(model =  out$best_bic$type_model, 
                                             garchOrder = c(out$best_bic$lag_arch,
                                                            out$best_bic$lag_garch)),
                       mean.model = list(armaOrder = c(out$best_bic$lag_ar, 
                                                       out$best_bic$lag_ma)),
                       distribution = 'std')


# estimate best garch model by AIC (used in next section)
best_spec_aic = ugarchspec(variance.model = list(model =  out$best_aic$type_model, 
                                             garchOrder = c(out$best_aic$lag_arch,
                                                            out$best_aic$lag_garch)),
                       mean.model = list(armaOrder = c(out$best_aic$lag_ar, 
                                                       out$best_aic$lag_ma)),
                       distribution = 'std')



my_best_garch <- ugarchfit(spec = best_spec_bic, 
                           data = df_prices$log_ret)

my_best_garch2 <- ugarchfit(spec = best_spec_aic, 
                           data = df_prices$log_ret)

my_best_garch

my_best_garch2



write_rds(my_best_garch, 'data/garch_model.rds')
########################################################################

my_garch = my_best_garch

set.seed(8008135) 

graphics.off()

library(foreach)
library(doParallel)

# Register parallel backend
no_cores <- detectCores() - 1  # Using one less than the total number of cores
registerDoParallel(cores = no_cores)

# Parallel processing of simulations
n_sim <- 5000
n_days_ahead <- 5*365 # Number of days ahead to simulate

results <- foreach(i_sim = 1:n_sim, .packages = c("tidyverse", "rugarch")) %dopar% {
  do_sim(n_sim = 1, 
         n_t = n_days_ahead, 
         my_garch, 
         df_prices = df_prices)
}

# Combine the results
df_sim <- bind_rows(results)
glimpse(df_sim)

# calculate probabilities of reaching peak value
tab_prob <- df_sim %>%
  group_by(ref_date) %>%
  summarise(prob = mean(sim_price > max(df_prices$price.adjusted)))


n_years_back <- 4
df_prices_temp <- df_prices %>%
  dplyr::filter(ref.date > max(ref.date) - n_years_back*365)

my_garch_name <- toupper(as.character(my_garch@model$modeldesc$vmodel))

p4 <- ggplot() + 
  geom_line(data = df_prices_temp, 
            aes(x = ref.date, y = price.adjusted), color = 'black', size = 0.75)  + 
  geom_line(data = df_sim, 
            aes(x = ref_date, 
                y = sim_price, 
                group = i_sim),
            color = 'red', 
            size = 0.25,
            alpha = 0.015) + 
  geom_hline(yintercept = max(df_prices_temp$price.adjusted)) + 
  labs(title = paste0('Price Projections of ', series_name),
       subtitle = paste0('Total of ', n_sim, ' simulations based on a ',
                         my_garch_name, 
                         ' model selected by BIC'),
       caption = 'Data from Yahoo Finance',
       x = '',
       y = 'Value') + 
  ylim(c(0.75*min(df_prices_temp$price.adjusted), 
         1.25*max(df_prices_temp$price.adjusted))) + 
  xlim(c(max(df_prices_temp$ref.date) - n_years_back*365,
         max(df_prices_temp$ref.date) + 5*365) )

ggsave('figs/fig06.png')


graphics.off()

my_idx_date <- first(which(tab_prob$prob > 0.5))
df_date <- tibble(idx = c(first(which(tab_prob$prob > 0.05)),
                          first(which(tab_prob$prob > 0.1)),
                          first(which(tab_prob$prob > 0.14)),
                          first(which(tab_prob$prob > 0.2))),
                  ref_date = tab_prob$ref_date[idx],
                  prob = tab_prob$prob[idx],
                  my_text = paste0(format(ref_date, '%m/%d/%Y'),
                                   '\nprob = ', scales::percent(prob) ) )

df_textbox <- tibble(ref_date = df_date$ref_date[2],
                     prob = 0.25,
                     label = paste0('According to the estimated _', my_garch_name, '_ model, ', 
                                    'the chances of asset **', series_name, '** to reach ',
                                    'its historical peak value of ', 
                                    format(max(df_prices$price.adjusted), 
                                           big.mark = ',',
                                           decimal.mark = '.'),
                                    ' are higher than 50% at ', format(ref_date, '%m/%d/%Y'), '.') )


p4 <- ggplot(tab_prob, aes(x = ref_date, y = prob) ) + 
  geom_line(size = 2) + 
  labs(title = paste0('Probabilities of ', series_name, ' Reaching its Historical Peak'),
       subtitle = paste0('Calculations based on simulations of ',
                         my_garch_name, 
                         ' model'),
       x = '',
       y = 'Probability') + 
  scale_y_continuous(labels = scales::percent) + 
  geom_point(data = df_date,
             aes(x = ref_date, y = prob), size = 5, color = 'red') + 
  geom_text(data = df_date, aes(x = ref_date, y = prob, 
                                label = my_text), 
            nudge_x = nrow(tab_prob)*0.085,
            nudge_y = -0.05,
            color ='red', check_overlap = TRUE) + 
  geom_textbox(data = df_textbox, 
               mapping = aes(x = ref_date, 
                             y = prob, 
                             label = label),
               width = unit(0.5, "npc"),
               #fill = "cornsilk",
               hjust = 0) + 
  theme_bw(base_family = "TT Times New Roman")


ggsave('figs/fig06_b.png')

################################################################################
#expirimental

quantiles_df <- df_sim %>%
  group_by(ref_date) %>%
  summarise(lower = quantile(sim_price, probs = 0.05),
            upper = quantile(sim_price, probs = 0.95),
            median = median(sim_price))

# Plotting
p5 <- ggplot(quantiles_df, aes(x = ref_date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = median), color = "red") +
  labs(title = "Probability Bands for Simulated Prices", 
       x = "Date", 
       y = "Price")

# Add the actual price data if available
if("price.adjusted" %in% names(df_prices)) {
  p5 <- p5 + geom_line(data = df_prices, aes(x = ref.date, y = price.adjusted), color = "black")
}

ggsave('figs/fig06_c.png')

# Print the plot
print(p5)
