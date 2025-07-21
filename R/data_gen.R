#' Data Generation
#'
#' @param n sample size
#' @param p number of variable
#' @param ecc error covariance
#' @param SNRdB SNR in dB
#' @param model Sinewave or random noise model
#'
#' @return A list of y, x, and error
#' @export
#'
dataGEN <- function(n, p, ecc, SNRdB, model="sine") {
  # dataGEN function
  # This function is to generate orthogonal y and e
  #
  # INPUT
  #   n = data length (scalar)
  #   p: number of datasets (scalar)
  #   ecc = error cross-correlation (scalar, [0,1])
  #   SNRdB = signal-to-noise ratio in dB (scalar)
  #
  # OUTPUT
  #   y = signal (nx1)
  #   e = error (nxp)

  # Estimation
  SNR <- 10^(SNRdB/10)
  b <- combn(1:p, 2)

  # Error covariance matrix: EeeT
  EeeT <- diag(runif(p)) # error variances
  for (i in 1:ncol(b)) {
    EeeT[b[1,i], b[2,i]] <- sqrt(EeeT[b[1,i], b[1,i]] * EeeT[b[2,i], b[2,i]]) * ecc
  }

  # Flip the matrix to make it symmetric
  EeeT <- (EeeT + t(EeeT)) - diag(diag(EeeT))

  # Generating y and e
  Ey2 <- mean(diag(EeeT) * SNR) # Signal power based on given SNR and EeeT

  m <- matrix(0, p + 1, p + 1) # covariance matrix of [y, e]
  m[1,1] <- Ey2
  m[2:(p + 1), 2:(p + 1)] <- EeeT

  if(model=="rand"){
    ye <- matrix(rnorm(n * (p + 1)), n, p + 1) # [y, e]

  } else if(model=="sine") {
    fs = 50
    dt = 1/fs
    t = seq(0, dt*(n-1), by=dt)
    y = 1.0*sin(2*pi*t+runif(1,0,2*pi))

    e =  matrix(rnorm(n * p), n, p)
    ye = cbind(y, e)
  } else {
    message('No such model!')
  }


  ye <- scale(ye, center = TRUE)
  ye <- ye %*% chol(cov(ye))  # scale by the Cholesky factor of the covariance matrix
  ye <- ye %*% chol(m)  # scale by the Cholesky factor of m

  y <- ye[, 1]
  e <- ye[, 2:(p + 1)]

  return(list(y = y, e = e))
}
