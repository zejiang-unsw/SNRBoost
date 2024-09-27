rm(list = ls())
#graphics.off()


# Random seed
set.seed(101)

### Step 1: Data parameters
p <- 3 # number of datasets
n <- 1000 # data length
eta <- matrix(rep(1, p), ncol = 1)

######## scaling factor (a)##########
a <- matrix(runif(p), ncol=1) # non-1-vector (for test)

ecc <- 0.3 # error cross-correlation [0,1]
SNRdB <- 0.1 # SNR in dB

### Step 2: Synthetic data generation
# Signal and error: y and e
model <- switch(2, "rand", "sine")
generated_data <- dataGEN(n, p, ecc, SNRdB, model)
y <- generated_data$y
e <- generated_data$e
y <- matrix(y, ncol=1)

# Observation: x = a*signal + error
x <- y %*% t(a) + e

plot.ts(cbind(y,x))
#------------------------------------------------------------------------------#
# Signal power and covariance matrices of e and x
Ey2 <- var(y)
Ey2 <- Ey2[1]
EeeT <- cov(e)
ExxT <- cov(x)
N <- EeeT / Ey2 # error-to-signal ratio

Ey2_est <- Ey2 * 0.5
Ey2_est <- Ey2_est[1]

ECVest_result <- ECVest(ExxT)
EeeT_est <- ECVest_result$EeeT_est

SNRest_result <- SNRest(ExxT, Ey2_est)
N_est <- SNRest_result$N_est
a_est <- SNRest_result$a_est

# merging by three methods
ue <- (1/p) * matrix(rep(1, p), nrow=p, ncol=1)  # equal weight
uw_est <- WA(EeeT_est)
us_est <- SNRopt(N_est, a_est)

# merged prod
ye <- x%*%ue
yw <- x%*%uw_est
ys <- x%*%us_est

wf <- "haar"
boundary <- "periodic"
if(wf!="haar") v <- as.integer(parse_number(wf)/2) else v <- 1
J <- floor(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)
J %>% print()
out <- SNRopt_freq(x=x, mode="DWT", wf=wf, J=J, boundary=boundary, theta = 0.1)
us_est_WT <- out$weight
yn <- out$merged

# RMSE of data merged by estimated parameters
y_m <- cbind(y,ye,yw,ys,yn)

plot.ts(y_m)

# metrics----
MSE_est <- sapply(1:4, function(i) mean((y-y_m[,i])^2))
MSE_est %>% print()

# Pearson correlation of data merged by estimated parameters
R_est <- cor(y,y_m)
R_est %>% print()


