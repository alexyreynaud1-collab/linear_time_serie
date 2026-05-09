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

# Estimation of candidate models (CSS method doesn't give AIC/BIC values)
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


# 95% confidence region for (Z_{T+1}, Z_{T+2})

phi    <- as.numeric(coef(arma_ols)[1:5])   # AR coefficients
theta  <- as.numeric(coef(arma_ols)[6:10])  # MA coefficients
sigma2 <- arma_ols$sigma2

# Compute psi_1 via the MA(inf) recursion:
#    psi_1 = theta_1 + phi_1 * psi_0,  psi_0 = 1
psi1 <- theta[1] + phi[1]

# Covariance matrix of (e_{T+1}, e_{T+2})
sigma1_sq <- sigma2
sigma2_sq <- sigma2 * (1 + psi1^2)
sigma12   <- sigma2 * psi1
Sigma     <- matrix(c(sigma1_sq, sigma12,
                      sigma12,   sigma2_sq), nrow = 2)

# Point forecasts
fc   <- predict(arma_ols, n.ahead = 2)
Z_T1 <- as.numeric(fc$pred[1])
Z_T2 <- as.numeric(fc$pred[2])

# Chi-squared threshold
chi2_95 <- qchisq(0.95, df = 2)

# Parametric ellipse via eigen decomposition of Sigma
eig  <- eigen(Sigma)
a    <- sqrt(chi2_95 * eig$values[1])   # semi-axis 1
b    <- sqrt(chi2_95 * eig$values[2])   # semi-axis 2
V    <- eig$vectors                      # rotation matrix

t_seq       <- seq(0, 2 * pi, length.out = 1000)
ellipse_std <- rbind(a * cos(t_seq), b * sin(t_seq))
ellipse_rot <- V %*% ellipse_std

ellipse_x <- Z_T1 + ellipse_rot[1, ]
ellipse_y <- Z_T2 + ellipse_rot[2, ]

# Individual 95% marginal intervals (for reference)
ci_z1 <- Z_T1 + c(-1, 1) * qnorm(0.975) * sqrt(sigma1_sq)
ci_z2 <- Z_T2 + c(-1, 1) * qnorm(0.975) * sqrt(sigma2_sq)

par(mfrow = c(1, 1))
plot(ellipse_x, ellipse_y,
     type = "l", lwd = 2, col = "steelblue",
     xlab = expression(Z[T+1]),
     ylab = expression(Z[T+2]),
     main = "95% Joint Confidence Region for (Z[T+1], Z[T+2])",
     asp = 1)

# Marginal CI lines
abline(v = ci_z1, lty = 2, col = "tomato")
abline(h = ci_z2, lty = 2, col = "green")

# Point forecast
points(Z_T1, Z_T2, pch = 19, col = "red", cex = 1.4)

legend("topright",
       legend = c("95% joint ellipse",
                  "Point forecast",
                  expression(paste("Marginal CI  ", Z[T+1])),
                  expression(paste("Marginal CI  ", Z[T+2]))),
       col  = c("steelblue", "red", "tomato", "green"),
       lwd  = c(2, NA, 1, 1),
       lty  = c(1, NA, 2, 2),
       pch  = c(NA, 19, NA, NA),
       bty  = "n")



# Forecast plot to compare the 40 last observations with the
#T+1 and T+2 predictions
library(lmtest)

path  <- "valeurs_mensuelles.csv"
data  <- read.csv2(path)
data  <- data[data$Codes == "A", ]
colnames(data) <- c("date", "index", "Codes")
data$date  <- as.Date(ym(data$date))
data$index <- as.numeric(data$index)

# Sort time (the CSV is newest → oldest)
data <- data[order(data$date), ]
data$time <- 1:nrow(data)

# Transformations
data$logindex     <- log(data$index)
data$logdiffindex <- c(NA, diff(data$logindex))   # Z_t = Δlog(X_t)

z_series  <- na.omit(data$logdiffindex)         
all_dates <- data$date[!is.na(data$logdiffindex)] 

T_len <- length(z_series)   # = 433, last date = 2026-02

# ARMA(5,5) model
arma_final <- arima(z_series, order = c(5, 0, 5), method = "CSS")

phi    <- as.numeric(coef(arma_final)[1:5])
theta  <- as.numeric(coef(arma_final)[6:10])
sigma2 <- arma_final$sigma2

psi1 <- theta[1] + phi[1]   # MA(∞) recursion: ψ_1 = θ_1 + φ_1·ψ_0

sigma1_sq <- sigma2
sigma2_sq <- sigma2 * (1 + psi1^2)
sigma12   <- sigma2 * psi1

# Forecast T+1 (March 2026) and T+2 (April 2026)
fc   <- predict(arma_final, n.ahead = 2)
Z_T1 <- as.numeric(fc$pred[1])
Z_T2 <- as.numeric(fc$pred[2])

date_last <- all_dates[T_len]          # 2026-02
date_T1   <- date_last %m+% months(1) # 2026-03
date_T2   <- date_last %m+% months(2) # 2026-04

# IC at 95%
z95   <- qnorm(0.975)
ci_T1 <- Z_T1 + c(-1, 1) * z95 * sqrt(sigma1_sq)
ci_T2 <- Z_T2 + c(-1, 1) * z95 * sqrt(sigma2_sq)


# Visualisation of the last 40 observations with the predicted values 

n_tail     <- 40
idx_tail   <- (T_len - n_tail + 1):T_len
z_tail     <- z_series[idx_tail]
dates_tail <- all_dates[idx_tail]

all_vals   <- c(z_tail, ci_T1, ci_T2)
ylim_range <- c(min(all_vals) * 1.3, max(all_vals) * 1.3)

xlim_range <- c(min(dates_tail), date_T2 + 150)

plot(dates_tail, z_tail,
     type  = "l", lwd  = 2, col  = "steelblue",
     xlim  = xlim_range,
     ylim  = ylim_range,
     xlab  = "Date",
     ylab  = expression(Z[t] == Delta * log(X[t])),
     main  = "Last 40 observations and 95% forecast intervals",
     xaxt  = "n")

all_x    <- c(dates_tail, date_T1, date_T2)
ax_ticks <- seq(min(dates_tail), date_T2, by = "2 months")
axis(1,
     at     = ax_ticks,
     labels = format(ax_ticks, "%b\n%Y"),
     cex.axis = 0.75, las = 2)

# Zero line and vertical separator
abline(h = 0,              lty = 1, col = "grey70")

segments(date_last, z_tail[n_tail], date_T1, Z_T1,
         lty = 2, col = "steelblue", lwd = 1.5)
segments(date_T1,   Z_T1,           date_T2, Z_T2,
         lty = 2, col = "tomato",    lwd = 1.5)

# CI on the predicted values
arrows(date_T1, ci_T1[1], date_T1, ci_T1[2],
       code = 3, angle = 90, length = 0.07,
       col = "steelblue", lwd = 2)
arrows(date_T2, ci_T2[1], date_T2, ci_T2[2],
       code = 3, angle = 90, length = 0.07,
       col = "tomato", lwd = 2)

# Predicted points
points(date_T1, Z_T1, pch = 19, col = "steelblue", cex = 1.5)
points(date_T2, Z_T2, pch = 19, col = "tomato",    cex = 1.5)

text(date_T1, Z_T1, labels = round(Z_T1, 4), pos = 4, cex = 0.8, col = "steelblue")
text(date_T2, Z_T2, labels = round(Z_T2, 4), pos = 4, cex = 0.8, col = "tomato")