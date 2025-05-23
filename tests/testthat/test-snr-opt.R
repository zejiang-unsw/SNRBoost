# This is an example of merging 3 synthetic zero mean datasets (time series)
# NOTE:
# (1) This example has a scaling factor (a) equal to 1, which is for fair
# comparison of the two merging methods (WA and SNRopt)
# (2) Since WA does not consider 'a' in its merging process, a fair
# comparison with SNRopt is not available unless the scaling factors are
# the same sa each other. On contrary, SNRopt allows arbitrary 'a' for
# weight computation. Any non-1 scaling factor can be tested in this
# example to see how the merging results are different.
# (3) Ey2 (signal power) shoud be reasonably estimated for SNRopt
# (e.g. reanalysis or climatology). Here, 0.5 of true Ey2 is tested.
#
# REFERENCE
# For more details, see:
#
# Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021).
# Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
# IEEE Trans Geosci Remote Sens
#
# If you use the methods presented in the paper and/or this example,
# please cite this paper where appropriate.
#

rm(list = ls())
cat("\014")

# Random seed
set.seed(4)

### Step 1: Data parameters
p <- 3 # number of datasets
n <- 1000 # data length
eta <- matrix(rep(1, p), ncol = 1)

######## scaling factor (a)##########
a <- matrix(runif(p), ncol=1) # non-1-vector (for test)
######################################

ecc <- 0.3 # error cross-correlation [0,1]
SNRdB <- 0.1 # SNR in dB

### Step 2: Synthetic data generation
# Signal and error: y and e
generated_data <- dataGEN(n, p, ecc, SNRdB)
y <- generated_data$y
e <- generated_data$e
y <- matrix(y, ncol=1)

# Observation: x = a*signal + error
x <- y %*% t(a) + e

# Signal power and covariance matrices of e and x
Ey2 <- var(y)
Ey2 <- Ey2[1]
EeeT <- cov(e)
ExxT <- cov(x)
N <- EeeT / Ey2 # error-to-signal ratio

# Functions of MSE and R2 by linear combination using u
MSE_func <- function(u) apply(t(u) %*% ExxT %*% u, 2, function(col) col - (2 * Ey2 * t(u) %*% a)) + Ey2
R2_func <- function(u) Ey2 * (t(u) %*% a %*% t(a) %*% u) / (t(u) %*% ExxT %*% u)

# RMSE of observations
MSE_ori <- diag(MSE_func(diag(p)))

# Pearson correlation of observations
R2_ori <- diag(R2_func(diag(p)))

# Printing metrics
cat("+ Metrics for original data\n")
cat(" * MSE for x:", round(MSE_ori, 3), "\n")
cat(" * R2 for x:", round(R2_ori, 3), "\n")

### Step 3: Merging using true parameters
# merging by three methods
uw <- WA(EeeT) # weighted average
us <- SNRopt(N, a) # SNR-opt
ur <- maxR(a, ExxT) # maxR weight
ue <- (1/p) * matrix(rep(1, p), nrow=p, ncol=1)  # equal weight

nume <- (t(eta) %*% (solve(EeeT) %*% eta))
deno <- (1 / Ey2) + t(a) %*% (solve(EeeT) %*% a)
scaling_factor <- c(as.numeric(nume / deno)*a)
s <- diag(scaling_factor)

# RMSE of data merged by true parameters
MSE_true <- diag(MSE_func(cbind(uw, us, ur, ue)))

# Pearson correlation of data merged by true parameters
R2_true <- diag(R2_func(cbind(uw, us, ur, ue)))

# Printing metrics
cat("\n+ Metrics for merged data by 'true' parameters\n")
cat(" * MSE for WA, SNRopt, maxR, EW:", round(MSE_true, 3), "\n")
cat(" * R2 for WA, SNRopt, maxR, EW:", round(R2_true, 3), "\n")

### Step 4: Merging using estimated parameters
# Estimation of merging statistics
Ey2_est <- Ey2 * 0.5
Ey2_est <- Ey2_est[1]

ECVest_result <- ECVest(ExxT)
EeeT_est <- ECVest_result$EeeT_est
theta_est <- ECVest_result$theta_est
rho2_est <- ECVest_result$rho2_est

SNRest_result <- SNRest(ExxT, Ey2_est)
N_est <- SNRest_result$N_est
a_est <- SNRest_result$a_est

# merging by three methods
uw_est <- WA(EeeT_est)
us_est <- SNRopt(N_est, a_est)
ur_est <- maxR(a_est, ExxT)

nume <- (t(eta) %*% (solve(EeeT) %*% eta))
deno <- (1 / Ey2) + t(a_est) %*% (solve(EeeT) %*% a_est)
scaling_factor <- c(as.numeric(nume / deno)*a_est)
s_est <- diag(scaling_factor)

# merged prod
ye <- x%*%ue
yw <- x%*%uw_est
yr <- x%*%ur_est
ys <- x%*%us_est

# RMSE of data merged by estimated parameters
y_m <- cbind(yw,yr,ys,ye)

# RMSE of data merged by estimated parameters
MSE_est <- diag(MSE_func(cbind(uw_est, us_est, ur_est, ue)))

# Pearson correlation of data merged by estimated parameters
R2_est <- diag(R2_func(cbind(uw_est, us_est, ur_est, ue)))
R_est <- cor(y, y_m)^2
R_est %>% print()

# Printing metrics
cat("\n+ Metrics for merged data by 'estimated' parameters\n")
cat(" * MSE for WA, SNRopt, maxR, EW:", round(MSE_est, 3), "\n")
cat(" * R2 for WA, SNRopt, maxR, EW:", round(R2_est, 3), "\n")

### Step 5: Plotting merging results
MSE_results <- rbind(matrix(mean(MSE_ori), nrow=1, ncol=4), t(MSE_true), t(MSE_est))
R2_results <- rbind(matrix(mean(R2_ori), nrow=1, ncol=4), t(R2_true), t(R2_est))

colors <- c("red", "green", "blue")

par(mfrow = c(1, 2))

colnames(MSE_results) <- c("WA", "SNR", "maxR", "EW")
rownames(MSE_results) <- c("Ori","True prm","Est prm")
barplot(MSE_results, beside = TRUE, names.arg = rownames(data), col = colors, legend.text = colnames(data))
title(main = "MSE", xlab = "Methods", ylab = "MSE")

colnames(R2_results) <- c("WA", "SNR", "maxR", "EW")
rownames(R2_results) <- c("Ori","True prm","Est prm")
barplot(R2_results, beside = TRUE, names.arg = rownames(data), col = colors, legend.text = colnames(data))
legend("bottomleft", legend = rownames(R2_results), fill = colors, inset = c(0, 0.15), ncol = 1)
title(main = "R2", xlab = "Methods", ylab = "R2")
