
install.packages("usethis")
# Packages installation 
install.packages("haven")
install.packages("tseries")

library(haven)
library(dplyr)
library(tseries)

# Data importation
path <-"~/Desktop/ensae/S2/serie_010767832_25042026/valeurs_mensuelles.csv"
data <- read.csv2(path)

data <- rename(data,index = Indice.CVS.CJO.de.la.production.industrielle..base.100.en.2021....Industrie.pharmaceutique..NAF.rév..2..niveau.division..poste.21.)
data <- rename(data, date = Libellé)
data <- filter(data,Codes == "A")


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

## ARMA structure 
par(mfrow = c(1, 2))
acf(na.omit(data$logdiffindex), main = "ACF")
pacf(na.omit(data$logdiffindex), main = "PACF")

#choice of the orders
p <- 5
q <- 5

#Estimation of the coefficients
arma_ols <- arima(na.omit(data$logdiffindex), order = c(p, 0, q), method = "CSS") #CSS for the OLS
summary(arma_ols)

#Validity of the model
Box.test(residuals(arma_ols), lag = 20, type = "Ljung-Box")


