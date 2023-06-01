# sensitivity for signal of standard to sample concentration of standard
k_A <- function(signalStandard, concentrationStandard) {
  return(signalStandard/concentrationStandard)
}

# C_A Analyte concentration for single point internal standard
C_A_sis <- function(C_IS, K, S_IS, S_A) {
  return((C_IS/K)*(S_A/S_IS))
}

# K for single point internal standard
K_sis <- function(C_IS, C_A, S_IS, S_A) {
  return((C_IS/C_A) * (S_A/S_IS))
}

# C_A Analyte for calibration curve, a=slope, b=y intercept
C_A_cal <- function(S_A, a, b) {
  return((S_A-b)/a)
}

C_A_cal_unknown <- function(S_A, a, b) {
  return((a*S_A)+b)
}

# stdev for calibration curve regression
sdev_cal <- function(y, yhat) {
  stopifnot(length(y)==length(yhat))
  stopifnot(length(y)>2)
  return(sqrt(sum( (y - yhat )^2 )/(length(y)-2)))
}

# stdev for b of calibration curve
sdev_cal_b <- function(sdev_cal, x) {
  mu_s <- mean(x)
  return(sqrt((sdev_cal^2)/sum( (x - mu_s)^2 )))
}

# stdev for a of calibration curve
sdev_cal_a <- function(sdev_cal, x) {
  mu_s <- mean(x)
  n <- length(x)
  return(sqrt((sdev_cal^2 * sum(x^2))/ n*sum( (x - mu_s)^2 )))
}

sdev_cal_C_A <- function(sdev_cal, a, n_sample_repl, n_cal_stds, S_avg_sample, S_avg_stds, C_stds) {
  mu_C_stds <- mean(C_stds)
  return(sdev_cal/a * (sqrt(1/n_sample_repl + 1/n_cal_stds + ( (S_avg_sample-S_avg_stds)^2 ) / (a^2)*sum( (C_stds-mu_C_stds)^2 ) )))
}

# calculate percent recovery from (corrected) Analyte concentration C_A and initial concentration C_I
perc_recovery <- function(C_A, C_I) {
  return(100*(C_A/C_I))
}