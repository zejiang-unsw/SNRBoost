#' Discrete Cosine Transform with Penalized Least Squares
#'
#' @param x A numeric vector with potential missing values (`NA`). Represents the input time series to be smoothed or gap-filled.
#' @param s A non-negative numeric scalar specifying the smoothing parameter. Smaller values yield less smoothing. If `NULL`, the function uses generalized cross-validation (GCV) to select an optimal value. Default is `NULL`.
#' @param robust Logical indicating whether to apply robust reweighting to down-weight outliers. Default is `FALSE`.
#' @param tol Convergence tolerance for robust iteration. Currently not used but reserved for future enhancements. Default is `1e-6`.
#' @param max_iter Maximum number of iterations for robust reweighting. Default is `10`.
#'
#' @returns A numeric vector of the same length as `x`, with missing values filled and optionally smoothed.
#' @export
#'
#' @references Wang, G., Garcia, D., Liu, Y., et al. (2012). A three-dimensional gap filling method for large geophysical datasets: Application to global satellite soil moisture observations. Environmental modelling & software, 30, 139-142. doi:https://doi.org/10.1016/j.envsoft.2011.10.015
#'
#' @examples
#' # Generate synthetic time series with missing values
#' set.seed(123)
#' n <- 100
#' x <- sin(2 * pi * (1:n) / 20) + rnorm(n, sd = 0.3)
#' x[30:40] <- NA  # introduce missing values
#'
#' # Apply DCT-PLS with fixed smoothing
#' y_filled <- dct_pls(x, s = 1e-6)
#'
#' # Plot original and filled series
#' plot(x, type = "l", lty = 2, lwd = 2, col = "black", ylab = "Value", main = "DCT-PLS Example")
#' lines(y_filled, col = "blue", lwd = 1)
#' legend("topright", legend = c("Original", "DCT-PLS"), col = c("black", "blue"),
#'        lty = c(2, 1), lwd = c(1, 2))
dct_pls <- function(x, s = NULL, robust = FALSE, tol = 1e-6, max_iter = 10) {
  stopifnot(is.numeric(x))
  n <- length(x)
  is_na <- is.na(x)
  y <- x
  y[is_na] <- mean(x, na.rm = TRUE)

  # Select DCT implementation
  dct   <- function(v) fftw::DCT(v, type = 2, norm = "ortho")
  idct  <- function(v) fftw::IDCT(v, type = 2, norm = "ortho")

  # Precompute smoothing weights
  freqs <- (0:(n-1))
  lam <- (2 * cos(pi * freqs / n) - 2)
  penalty <- lam^2

  # GCV for optimal s
  if (is.null(s)) {
    gcv_fun <- function(logs) {
      s_trial <- exp(logs)
      w <- 1 / (1 + s_trial * penalty)
      y_hat <- idct(dct(y) * w)
      rss <- sum((y_hat[!is_na] - x[!is_na])^2)
      edf <- sum(w)
      gcv <- rss / (n * (1 - edf / n)^2)
      return(gcv)
    }
    opt <- stats::optimize(gcv_fun, interval = log(c(1e-6, 1e6)))
    s <- exp(opt$minimum)
  }

  # Fit with chosen s
  w <- 1 / (1 + s * penalty)
  y_hat <- idct(dct(y) * w)

  # Robust reweighting (optional)
  if (robust) {
    for (i in seq_len(max_iter)) {
      resid <- y_hat[!is_na] - x[!is_na]
      wts <- pmin(1, 4 / pmax(abs(resid), 1e-6))  # avoid division by 0

      # Update observed values only
      y[!is_na] <- wts * x[!is_na] + (1 - wts) * y_hat[!is_na]

      # Smooth again
      y_hat <- idct(dct(y) * w)
    }
  }

  return(y_hat)
}
