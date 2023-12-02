library(dplyr)
library(ggplot2)
library(tseries)
library(xts)
library(quantmod)
library(tidyquant)
library(vars)
library(forecast)
library(zoo)
library(lubridate)

interest10_rate <- tq_get("DGS10", get = "economic.data", from = "2004-09-01", to="2020-03-30")
inflation10_rate <- tq_get("CPIAUCSL", get= "economic.data", from = "2004-08-11", to="2020-03-30")

gold <- read.csv("C:/Users/justd/Desktop/Gold2004_2023.csv", header=T)
gold$Date <- as.Date.character(gold$Date, format= "%m/%d/%Y")


gold_data <- gold %>%
  dplyr::select(Date, Price)

gold_data$Price <- as.numeric(gsub(',', '', gold_data$Price))


interest10_rate <- interest10_rate %>% 
    dplyr::select(-symbol)
inflation10_rate <- inflation10_rate %>%
  dplyr::select(-symbol)


#Performed linear interpolation for missing values 
#Since, ill be working with VAR model, I didnt want to smooth out the NA so I can minimize distortion 
#Utilized linear interpolation since the NA are short gaps, better for VAR models

interest10_rate$price <- na.approx(interest10_rate$price)

#Creating a sequence of dates to linear interpolate the first day of month   
start_date <- as.Date("2004-12-01")
end_date <- as.Date("2020-03-01")

date_seq <- seq(from = start_date, to = end_date, by = "day")


gold_ts <- data.frame(Date = date_seq) %>%
  left_join(gold_data, by = 'Date') %>%
  mutate(Price = na.approx(Price, na.rm = F))

gold_ts <- gold_ts %>%
  mutate(date = floor_date(Date, "month")) %>%
  group_by(date) %>%
  summarize(Price = mean(Price, na.rm = T)) %>%
  ungroup()



interest10_rate <- data.frame(date = date_seq) %>%
  left_join(interest10_rate, by = 'date') %>%
  mutate(price = na.approx(price, na.rm = F))

interest10_rate <- interest10_rate %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarize(price = mean(price, na.rm = T)) %>%
  ungroup()


####### Gold Timeseries Data ######
gold_xts <- xts(gold_ts$Price, order.by=gold_ts$date)
gold_xts <- na.omit(gold_xts)
plot(gold_xts)
acf(gold_xts)

#Pre-stationary checking

#Transforming that data to make it stationary
gold_xts_diff <- diff(gold_xts)
gold_xts_diff <- na.omit(gold_xts_diff)
plot(gold_xts_diff)
acf(gold_xts_diff)

#To Determine number of lags for better model robustness
ts_length <- length(gold_xts) #gets the length of observation of ts data
max_lags = floor(12 * (ts_length/100) ^ 0.25) #Schwert Formulas to decide number of maximum lag values to avoid overfitting
aic_values = rep(NA, max_lags + 1) #The aic values is the the return of max_lags

for (k in 0:max_lags) {
  tryCatch({
    fit_gold = Arima(gold_xts_diff, order = c(k, 0, 0), include.constant = F)
    aic_values[k + 1] = fit_gold$aic},
      error = function(e) {
      cat("Error: ", e$message, "\n")
      aic_values[k + 1] = NA
 })
}

#AIC testing took a while to run, future points to consider is a package parallel computing
#parallel: Basically spreads the work across different CPU cores in your PC

#Based on the best score for AIC values it looks for the optimal lag
valid_aic <- !is.na(aic_values)
if (any(valid_aic)) {
  optimal_lags = which.min(aic_values[valid_aic]) - 1
  optimal_lags
} else {
  cat("No valid AIC values were found.")
}



#The output turned out to be zero as it means that there was no need for additional lags,...
#...Meaning the ACF shows that the series is already quite stationary 
#From results the total maximum number of lags were tested which came back as all being stationary 


adf_gold <- adf.test(gold_xts_diff, k = optimal_lags)
cat("Lags: ", k, "ADF Test Statistics: ",adf_gold$statistic, "P-value: ", adf_gold$p.value, "\n")
adf_gold



#################################################

######### Inflation Time series #################
inflation10_xts <- xts(inflation10_rate$price, order.by=inflation10_rate$date)
plot(inflation10_xts)
acf(inflation10_xts)



#Make inflation ts data stationary 
inflation_log <- log(inflation10_xts)

plot(inflation_log)
acf(inflation_log)

inflation_diff <- diff(inflation_log)
inflation_diff <- na.omit(inflation_diff)

plot(inflation_diff)
acf(inflation_diff)


#############################################################################################################

for (p in 0:max_lags_infl) {
   tryCatch({
     fit_inflation = Arima(inflation_diff, order = c(p, 0, 0), include.constant = F)
     aic_values_infl[p + 1] = fit_inflation$aic},
     
     error_infl = function(e) {
       aic_values_infl[p + 1] = NA
     })
}
 
#Based on the best score for AIC values it looks for the optimal lag
valid_aic_infl <- !is.na(aic_values_infl)
if (any(valid_aic_infl)) {
   optimal_lags_infl = which.min(aic_values_infl[valid_aic_infl]) - 1
   optimal_lags_infl
 } else {
   cat("No valid AIC values were found.")
}

residuals_infl <- residuals(fit_inflation)
plot(residuals_infl)
acf(residuals_infl)

inflation_adf <- adf.test(inflation_diff, k= optimal_lags_infl) #Now testing the stationary from transformed data using optimal lag values 
inflation_adf
 
#########################################################################################################################

######## Interest Time series ######################
interest10_xts <- xts(interest10_rate$price, order.by=interest10_rate$date)
interest10_xts <- na.omit(interest10_xts)

plot(interest10_xts)
acf(interest10_xts)

interest_diff <- diff(interest10_xts)
interest_diff <- na.omit(interest_diff)
plot(interest_diff)


#Determinig optimal number
ts_length_interest <- length(interest10_xts)
max_lags_inter = floor(12 * (ts_length_interest/100) ^ 0.25)
aic_values_inter = rep(NA, max_lags_inter + 1)

for (p in 0:max_lags_inter) {
  tryCatch({
    fit_interest = Arima(interest_diff, order = c(p, 0, 0), include.constant = F)
    aic_values_inter[p + 1] = fit_interest$aic},
    
    error_int = function(e) {
      aic_values_inter[p + 1] = NA
    })
}

#Based on the best score for AIC values it looks for the optimal lag
valid_aic_inter <- !is.na(aic_values_inter)
if (any(valid_aic_inter)) {
  optimal_lags_inter = which.min(aic_values_inter[valid_aic_inter]) - 1
  optimal_lags_inter
} else {
  cat("No valid AIC values were found.")
}



interest_adf <- adf.test(interest_diff, k = optimal_lags_inter)
interest_adf 

############## Vector AutoRegressive Model Estimating ###########################
gold_diff_sub <- gold_xts_diff['2005-01-01/2020-02-01']
inflation_diff_sub <- inflation_diff['2005-01-01/2020-02-01']
interest_diff_sub <- interest_diff['2005-01-01/2020-02-01']

merged.xts.data <- merge(gold_diff_sub, inflation_diff_sub, interest_diff_sub)
summary(merged.xts.data)

merged.ts.data <- as.ts(merged.xts.data) #changed xts to ts data type for VAR


#Choosing optimal lag order 
merged.ts.maxlag<- length(merged.ts.data) #gets the length of observations of ts data
merged.lag.max <-  floor(12 * (merged.ts.maxlag/100) ^ 0.25)

lag_selection <- VARselect(merged.ts.data, lag.max = merged.lag.max, type= 'both')
lag_selection$selection

optimal_lag.merged <- lag_selection$selection["AIC(n)"]

VAR.model <- VAR(merged.ts.data, p = optimal_lag.merged, type= 'both')

########Pre-Diagnostics of the model###############

#Portmanteau (Ljun-Box Test) - Checking for white noise
serial.test(VAR.model) #Checking for serial correlations in the residuals AKA Auto Correlations
#Based on Results there is no serial correlation meaning that the model captures the time series... 
#...dynamics adequately and the residuals don't influence each other

#data:  Residuals of VAR object VAR.model
#Chi-squared = 100.15, df = 126, p-value = 0.9566
#P-values is higher than 0.05, which we fail to reject null, meaning there so significant evidence of autocorrelation

##########################
#This impliees that the model is adequately capturing the pattern since the residuals do not have predictable pattern

VAR.model_data <- VAR.model$datamat

#Impulse Response Function - This will examine the response or reaction of one variable shocks another

irf_var.inflation <- irf(VAR.model, impulse = "inflation_diff_sub", response = "gold_diff_sub", n.ahead = 10 )
plot(irf_var.inflation)

#Horizontal Axis is the number of lags (periods) after initial lag.
#Vertical Axis is the response magnitude in terms of how much is gold expected to change from inflation shock
#Solid line estimated response from gold prices 
#Dashed are Confidence intervals, if it includes zero then its not statistically significant on that particular lag. 
#Overall since most of the CIs include zero then its not statistically signifcant 

irf_var.interest <- irf(VAR.model, impulse = "interest_diff_sub", response = "gold_diff_sub", n.ahead = 10 )
plot(irf_var.interest)

#This confirms the dynamic relationship between gold and interest rate from the plot showcasing that an initial negative response aligns with traditional views...
#Rising interest rates could dampen the appeal of non yielding assets like gold as an increase in interest rates leads to a decrease in gold prices 
#Since impact stabilizes as shown after 5 then the shock is not permanent, shown by conf int. converging towards zero.

###########################
#Forecast Error Variance Decomposition (FEVD)

fevd.VAR <- fevd(VAR.model, n_ahead = 10)
plot(fevd.VAR)

#For the most part Gold is affected by its own shocks.

#####################################
#Constant variance overtime Residual Diagnostics 

arch.model <- arch.test(VAR.model)
arch.model

#data:  Residuals of VAR object VAR.model
#Chi-squared = 259.86, df = 180, p-value = 9.047e-05

#The p value is significantly below 0.05 which means there is no evidence to suggest ARCH effect
#The residuals are considered homeostatic

#Checking for model stability by seeing if the roots of the characteristic polynomial lie inside
#the unit circle 

VAR.stability <- vars::stability(VAR.model)
VAR.stability


#Out-of-samples estimations
#This will create forecast scenarios of the gold prices trend up until when COVID technically ended or stopped affected prices

n.ahead.forecast = 22 #About 2 years 

VAR.gold_estimations <- predict(VAR.model, n.ahead = n.ahead.forecast)


#extracting the forecasting values 

gold.forecast_values <- VAR.gold_estimations$fcst$gold_diff_sub[, "fcst"]


#Comparative Analysis actual observed gold prices with forecasted pre-COVID gold prices

#Creating a sequence of dates to linear interpolate the first day of month

start_covid = as.Date("2020-02-01")
end_covid = as.Date("2021-12-01")
 
date_seq_covid <- seq(from = start_covid, to = end_covid, by = "day")

gold_covid <- data.frame(Date = date_seq_covid) %>%
   left_join(gold_data, by = 'Date') %>%
   mutate(Price = na.approx(Price, na.rm = F))
 
gold_covid <- gold_covid %>%
   mutate(date = floor_date(Date, "month")) %>%
   group_by(date) %>%
   summarize(Price = mean(Price, na.rm = T)) %>%
   ungroup()

#So the predicted values outputted from the VAR model is the predictions of the changes from one period to the next
#Not the level themselves, in other words not the series' level used to forecast
#To go back and compared at the same level, we need to un-difference the predicted values

last_obs_value <- tail(gold_covid[gold_covid$date <= start_covid, "Price"], 1)
forecasted_levels <- cumsum(c(last_obs_value, gold.forecast_values))

comparative_analysis <- data.frame(
  Time = seq(from = as.Date("2020-02-01"), by = 'month', length.out = length(forecasted_levels)),
  Forecasted = forecasted_levels,
  Actual = gold_covid$Price
)

comparative_analysis <- comparative_analysis[-1, ]

plot(comparative_analysis$Time, comparative_analysis$Forecasted, type = 'l', col = 'blue', ylim = range(comparative_analysis[, 2:3]), ylab= "Gold Prices", xlab = 'Time')
lines(comparative_analysis$Time, comparative_analysis$Actual, col = 'red')
legend("topright", legend = c("Forecasted Estimates", "Actual Observed"), fill = c("blue", "red"))


#Mean Absolute Error 
#Measures the average magnitude of the errors in a set of forecasts basically calculates the accuracy in prediction 
#from the average over the test samples of absolute differences between predicted and actual at equal weight 

MAE <- function(actual, forecast) {
  mean(abs(actual - forecast))
}

attach(comparative_analysis)
mae.val <- MAE(Actual, Forecasted)
mae.val

#Mean Absolute Percentage Error (MAPE)

MAPE <- function(actual, forecast) {
  mean(abs((actual - forecast) / actual)) * 100
}

mape.val <- MAPE(Actual, Forecasted)
mape.val

#The result of 9.612435 suggests that the actual and observed values are approximately off by 9.61%
#Generally, value less than 10% is considered very good, as in fairly accurate 

