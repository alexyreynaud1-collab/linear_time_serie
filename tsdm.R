install.packages("usethis")
# Packages installation 
install.packages("haven")
install.packages("tseries")

library(haven)
library(dplyr)
library(tseries)

# Data importation
path <-"valeurs_mensuelles.csv"
data <- read.csv2(path)

data <- rename(data,index = Indice.CVS.CJO.de.la.production.industrielle..base.100.en.2021....Industrie.pharmaceutique..NAF.rév..2..niveau.division..poste.21.)
data <- rename(data, date = Libellé)
data <- filter(data,Codes == "A")

install.packages("lubridate")
library(lubridate)
data$date <- ym(data$date)
data$date <- as.Date(data$date)
data$time <- 1:nrow(data)
colnames(data)

# Part I : the data
#Plot of the serie before transformations 
plot(data$date, data$index, 
     type = "l",           
     main = "Index throughout time before transformations",
     xlab = "Time",
     ylab = "Index",
     col = "steelblue")

#Detrending of the time serie
# Stochastic trend in variance
data$index <- as.numeric(data$index)
data$logindex <- log(data$index)
colnames(data)

plot(data$date, data$logindex, 
     type = "l",           
     main = "Index after the stochastic detrending in variance",
     xlab = "Time",
     ylab = "Index",
     col = "steelblue")

#Stochastic trend in mean 
data$logdiffindex <- c(NA, diff(data$logindex))
colnames(data)
plot(data$date, data$logdiffindex, 
     type = "l",           
     main = "Index after the stochastic detrending (log differenciation)",
     xlab = "Time",
     ylab = "Index",
     col = "steelblue")

#Test of the stationarity 
adf.test(na.omit(data$logdiffindex))

# H0: The series has a unit root (non-stationary)
# We want to reject H0 → p-value < 0.05
cat("=== Phillips-Perron Test ===\n")
pp.test(na.omit(data$logdiffindex))

# H0: The series is stationary
# We want NOT to reject H0 → p-value > 0.05
cat("=== KPSS Test ===\n")
kpss.test(na.omit(data$logdiffindex), null = "Level")

## ARMA structure 
par(mfrow = c(1, 2))
acf(na.omit(data$logdiffindex), main = "ACF")
pacf(na.omit(data$logdiffindex), main = "PACF")

# Estimation of candidate models
arma44 <-arima(na.omit(data$logdiffindex), order = c(4, 0, 4), method = "ML") #CSS for the OLS
arma45 <- arima(na.omit(data$logdiffindex), order = c(4, 0, 5), method = "ML") #CSS for the OLS
arma55 <- arima(na.omit(data$logdiffindex), order = c(5, 0, 5), method = "ML") #CSS for the OLS
arma54 <- arima(na.omit(data$logdiffindex), order = c(5, 0, 4), method = "ML") #CSS for the OLS



install.packages("lmtest")   # à faire une seule fois
library(lmtest)
# Visualisation of coefficient significance. Every coefficient must be significant.
cat("=== Coefficients ARMA(5,4) ===\n")
print(coeftest(arma54))

cat("=== Coefficients ARMA(5,5) ===\n")
print(coeftest(arma55))

cat("=== Coefficients ARMA(4,5) ===\n")
print(coeftest(arma45))

cat("=== Coefficients ARMA(4,4) ===\n")
print(coeftest(arma44))

# Ljung-Box test (H0 : non correlated residuals = good model)
cat("=== ARMA(5,4) ===\n")
Box.test(residuals(arma54), lag = 20, type = "Ljung-Box")

cat("=== ARMA(5,5) ===\n")
Box.test(residuals(arma55), lag = 20, type = "Ljung-Box")

cat("=== ARMA(4,5) ===\n")
Box.test(residuals(arma45), lag = 20, type = "Ljung-Box")

cat("=== ARMA(4,4) ===\n")
Box.test(residuals(arma44), lag = 20, type = "Ljung-Box")

# Comparison AIC / BIC
comparison <- data.frame(
  Modele = c("ARMA(4,5)", "ARMA(5,5)", "ARMA(4,4)", "ARMA(5,4)"),
  AIC    = c(AIC(arma45), AIC(arma55), AIC(arma44), AIC(arma54)),
  BIC    = c(BIC(arma45), BIC(arma55), BIC(arma44), BIC(arma54))
)
cat("=== Comparaison AIC / BIC ===\n")
print(comparison[order(comparison$AIC), ])

# final choice of the orders for the next section
p <- 5
q <- 5

#Estimation of the coefficients
arma_ols <- arima(na.omit(data$logdiffindex), order = c(p, 0, q), method = "CSS") #CSS for the OLS
summary(arma_ols)

cat("=== Coefficients ARMA(5,5) ===\n")
print(coeftest(arma_ols))

#Validity of the model
Box.test(residuals(arma_ols), lag = 20, type = "Ljung-Box")

# Remark : the model is not well adjusted since the coefficients aren't significant...

# The model for the series is an ARIMA(5,1,5).

