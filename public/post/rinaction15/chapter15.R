
# Chapter 15. Time series


# This chapter covers
#   Creating a time series
#   Decomposing a time series into components
#   Developing predictive models
#   Forecasting future values


# Remove most objects from the working environment
rm(list = ls())
options(stringsAsFactors = F)

# 15.1. Creating a time-series object in R
# code listing 15.1. Creating a time-series object
sales <- c(18, 33, 41,  7, 34, 35, 24, 25, 24, 21, 25, 20,
           22, 31, 40, 29, 25, 21, 22, 54, 31, 25, 26, 35)

tsales <- ts(sales, start = c(2003, 1), frequency = 12)
tsales
#      Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
# 2003  18  33  41   7  34  35  24  25  24  21  25  20
# 2004  22  31  40  29  25  21  22  54  31  25  26  35

plot(tsales)  # figure 15.2
plot(tsales, type="o", pch=17, col="blue")  # figure 15.2-1
start(tsales)
# [1] 2003    1

end(tsales)
# [1] 2004   12

frequency(tsales)
# [1] 12

tsales.subset <- window(tsales, start=c(2003, 5), end=c(2004, 6))
tsales.subset
#      Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
# 2003                  34  35  24  25  24  21  25  20
# 2004  22  31  40  29  25  21


#===============================================================
# 15.2.1. Smoothing with simple moving averages
# code listing 15.2 Simple moving averages
# install.packages("forecast")
library(forecast)
opar <- par(no.readonly = T)
par(mfrow=c(2,2))
ylim <- c(min(Nile), max(Nile))
plot(Nile, main="Raw time series")
plot(ma(Nile, 3), main="Simple Moving Averages (k=3)", ylim=ylim)
plot(ma(Nile, 3), main="Simple Moving Averages (k=7)", ylim=ylim)
plot(ma(Nile, 7), main="Simple Moving Averages (k=15)", ylim=ylim)
par(opar)


# 15.2.2. Seasonal decomposition
# code listing 15.3. Seasonal decomposition using stl()
opar <- par(no.readonly = T)
par(mfrow=c(2,1))
plot(AirPassengers)
lAirPassengers <- log(AirPassengers)
plot(lAirPassengers, ylab="log(AirPassengers)")
par(opar)

fit <- stl(lAirPassengers, s.window = "period")
plot(fit)


head(fit$time.series)
#          seasonal trend remainder
# Jan 1949   -0.092   4.8    -0.019
# Feb 1949   -0.114   4.8     0.054
# Mar 1949    0.016   4.8     0.036
# Apr 1949   -0.014   4.8     0.040
# May 1949   -0.015   4.8    -0.025
# Jun 1949    0.110   4.8    -0.043


head(exp(fit$time.series))
#          seasonal trend remainder
# Jan 1949     0.91   125      0.98
# Feb 1949     0.89   125      1.06
# Mar 1949     1.02   125      1.04
# Apr 1949     0.99   126      1.04
# May 1949     0.99   126      0.98
# Jun 1949     1.12   126      0.96

# figure 15.7
opar <- par(no.readonly = T)
par(mfrow=c(2,1))
library(forecast)
monthplot(AirPassengers, xlab = "", ylab = "")
seasonplot(AirPassengers, year.labels = "T", main = "")
par(opar)

# 15.3.1. Simple exponential smoothing
# code listing 15.4. Simple exponential smoothing
library(forecast)
fit <- ets(nhtemp, model = "ANN")
fit
# ETS(A,N,N) 
# 
# Call:
#   ets(y = nhtemp, model = "ANN") 
# 
#    Smoothing parameters:
#      alpha = 0.1819 
# 
#    Initial states:
#      l = 50.2762 
# 
#    sigma:  1.1
# 
# AIC AICc  BIC 
# 266  266  272 
nhtemp
# Time Series:
#   Start = 1912 
# End = 1971 
# Frequency = 1 
# [1] 50 52 49 51 49 48 50 51 49 52 51 50 49 51 48 51 51 51 52 53 52 51 50 50 50 52 52 51 49
# [30] 52 51 51 52 52 52 51 51 54 51 53 53 55 52 52 51 53 50 53 52 52 50 51 52 51 52 51 52 52
# [59] 52 53

forecast(fit, 1)
#      Point Forecast Lo 80 Hi 80 Lo 95 Hi 95
# 1972             52    50    53    50    54

plot(forecast(fit, 1), xlab = "Year",
     ylab = expression(paste("Temperature (", degree*F,")",)),
     main = "New Haven Annual Mean Temperature")  # figure 15.8

accuracy(fit)
#                ME RMSE MAE  MPE MAPE MASE    ACF1
# Training set 0.15  1.1 0.9 0.24  1.7 0.75 -0.0064


#=================================================================
# 15.3.2. Holt and Holt-Winters exponential smoothing
# code listing 15.5. Exponential smoothing with level, slope, and seasonal components
library(forecast)
fit <- ets(log(AirPassengers), model = "AAA")
fit
# ETS(A,A,A) 
# 
# Call:
#   ets(y = log(AirPassengers), model = "AAA") 
# 
#    Smoothing parameters:
#      alpha = 0.6975 
#       beta  = 0.0031 
#      gamma = 1e-04 
# 
#    Initial states:
#      l = 4.7925 
#      b = 0.0111 
#      s = -0.1 -0.22 -0.079 0.056 0.2 0.21
#           0.11 -0.0081 -0.0059 0.022 -0.11 -0.084
# 
#    sigma:  0.038
# 
# AIC AICc  BIC 
# -207 -202 -157

accuracy(fit)
#                   ME  RMSE   MAE    MPE MAPE MASE  ACF1
# Training set -0.0018 0.036 0.028 -0.034 0.51 0.23 0.056

pred <- forecast(fit, 5)
pred
#          Point Forecast Lo 80 Hi 80 Lo 95 Hi 95
# Jan 1961            6.1   6.1   6.2   6.0   6.2
# Feb 1961            6.1   6.0   6.2   6.0   6.2
# Mar 1961            6.2   6.2   6.3   6.1   6.3
# Apr 1961            6.2   6.1   6.3   6.1   6.3
# May 1961            6.2   6.1   6.3   6.1   6.4

plot(pred, main = "Forecast for Air Travel",
     ylab = "Log(Airpassengers)", xlab = "Time")  # figure 15.9

pred$mean <- exp(pred$mean)
pred$lower <- exp(pred$lower)
pred$upper <- exp(pred$upper)
p <- cbind(pred$mean, pred$lower, pred$upper)
dimnames(p)[[2]] <- c("mean", "Lo 80", "Lo 95", "Hi 80", "Hi 95")
p
#          mean Lo 80 Lo 95 Hi 80 Hi 95
# Jan 1961  450   429   418   473   485
# Feb 1961  443   417   404   470   485
# Mar 1961  511   477   460   548   568
# Apr 1961  502   465   446   542   565
# May 1961  506   465   445   551   576

#===============================================================
# 15.3.3. The ets() function and automated forecasting
# code lisitng Automatic exponential forecasting with ets()
library(forecast)
fit <- ets(JohnsonJohnson)
fit
# ETS(M,A,A) 
# 
# Call:
#   ets(y = JohnsonJohnson) 
# 
#    Smoothing parameters:
#      alpha = 0.2776 
#       beta  = 0.0636 
#      gamma = 0.5867 
# 
#    Initial states:
#      l = 0.6276 
#      b = 0.0165 
#      s = -0.23 0.19 -0.0074 0.045
# 
#    sigma:  0.092
# 
# AIC AICc  BIC 
# 164  166  186 

plot(forecast(fit), main = "Johnson & Johnson Forecasts",
     ylab = "Quarterly Earnings (Dollars)", xlab = "Time", flty=2)  # figure 15.10


fit2 <- ets(JohnsonJohnson, model = "MMM")
fit2
# ETS(M,M,M) 
# 
# Call:
#   ets(y = JohnsonJohnson, model = "MMM") 
# 
#    Smoothing parameters:
#      alpha = 0.2195 
#       beta  = 0.0341 
#      gamma = 0.5675 
# 
#    Initial states:
#      l = 0.625 
#      b = 1.0265 
#      s = 0.7 1.3 0.98 1.1
# 
#    sigma:  0.091
# 
# AIC AICc  BIC 
# 164  166  186 


#========================================
# 15.4.2. ARMA and ARIMA models
# code listing 15.7. Transforming the time series and assessing stationarity
library(forecast)
library(tseries)

# figure 15.11
opar <- par(no.readonly = T)
par(mfrow=c(2,1))
plot(Nile)
ndiffs(Nile)
# [1] 1

dNile <- diff(Nile)
plot(dNile)
adf.test(dNile)
#         Augmented Dickey-Fuller Test
# 
# data:  dNile
# Dickey-Fuller = -7, Lag order = 4, p-value = 0.01
# alternative hypothesis: stationary
par(opar)

# figure 15.12
opar <- par(no.readonly = T)
par(mfrow=c(2,1))
Acf(dNile)
Pacf(dNile)


# code listing 15.8. Fitting an ARIMA model
library(forecast)
fit <- arima(Nile, order = c(0, 1, 1))
fit
# Call:
# arima(x = Nile, order = c(0, 1, 1))
# 
# Coefficients:
#         ma1
#       -0.73
# s.e.   0.11
# 
# sigma^2 estimated as 20600:  log likelihood = -633,  aic = 1269

accuracy(fit)
#               ME RMSE MAE  MPE MAPE MASE ACF1
# Training set -12  143 112 -3.6   13 0.84 0.12


# code listing 15.9. Evaluating the model fit
# figure 15.13
qqnorm(fit$residuals)
qqline(fit$residuals)
Box.test(fit$residuals, type = "Ljung-Box")
#         Box-Ljung test
# 
# data:  fit$residuals
# X-squared = 1, df = 1, p-value = 0.2


# code listing 15.10. Forecasting with an ARIMA model
forecast(fit, 3)
#      Point Forecast Lo 80 Hi 80 Lo 95 Hi 95
# 1971            798   614   982   517  1080
# 1972            798   608   989   507  1090
# 1973            798   602   995   498  1099


plot(forecast(fit, 3), xlab = "Year", ylab = "Annual Flow") # figure 15.14


# 15.4.3. Automated ARIMA forecasting
# code listing 15.11. Automated ARIMA forecasting
library(forecast)
fit <- auto.arima(sunspots)
fit
# Series: sunspots 
# ARIMA(2,1,2) 
# 
# Coefficients:
#        ar1     ar2    ma1    ma2
#       1.35  -0.396  -1.77  0.810
# s.e.  0.03   0.029   0.02  0.019
# 
# sigma^2 estimated as 244:  log likelihood=-11746
# AIC=23501   AICc=23501   BIC=23531

forecast(fit, 3)
#          Point Forecast Lo 80 Hi 80 Lo 95 Hi 95
# Jan 1984             40    20    60   9.8    71
# Feb 1984             41    18    64   6.0    77
# Mar 1984             40    15    64   2.2    77

accuracy(fit)
#                  ME RMSE MAE MPE MAPE MASE   ACF1
# Training set -0.027   16  11 NaN  Inf 0.48 -0.011