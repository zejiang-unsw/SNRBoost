#' SNRopt_freq
#'
#' @param x parent products
#' @param mode wavelet decomposition mode
#' @param wf wavelet filter
#' @param J max level of decomposition level
#' @param boundary boundary condition of wavelet transform
#' @param theta factor for Ey2 estimation
#'
#' @return A list of weight, merged products.
#' @export
#'
#' @import waveslim
#' @import WASP
#' @import readr
#'
SNRopt_freq <- function(y=NULL, x, mode="DWT", wf="haar", J, pad="zero", boundary="periodic", theta=0.1,
                        option="Est", scale_factor=T) {

  n <- nrow(x)
  p <- ncol(x)

  # wavelet transform
  if(wf!="haar") v <- as.integer(parse_number(wf)/2) else v <- 1
  if(!exists("J")) J <- floor(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)

  x_WTs <- array(dim=c(n, p, J+1))
  for(i_x in 1:p){

    if(mode=="DWT"){
      x1 <- padding(x[,i_x], pad = pad)
      tmp <- mra(x1, wf=wf, J=J, method="modwt", boundary=boundary)
      x_WTs[,i_x, ] <- matrix(unlist(lapply(tmp, function(z) z[1:n])), ncol = J + 1, byrow = FALSE)
    } else if(mode=="MODWT"){
      tmp <- modwt(x[,i_x], wf=wf, n.levels = J, boundary = boundary)
      #summary(tmp) %>% print()
      x_WTs[,i_x, ] <- matrix(unlist(tmp), ncol = J + 1, byrow = FALSE)
    }

  }

  if(mode=="DWT"){
    y1 <- padding(y, pad = pad)
    tmp <- mra(y1, wf=wf, J=J, method="modwt", boundary=boundary)
    y_WT <- matrix(unlist(lapply(tmp, function(z) z[1:n])), ncol = J + 1, byrow = FALSE)
  } else if(mode=="MODWT"){
    tmp <- modwt(y, wf=wf, n.levels = J, boundary = boundary)
    #summary(tmp) %>% print()
    y_WT <- matrix(unlist(tmp), ncol = J + 1, byrow = FALSE)
  }


  # weight estimation
  u_est_WT <- matrix(nrow=J+1, ncol=p)
  y_est_WT <- matrix(nrow=n, ncol=J+1)
  for(i_lev in 1:nrow(u_est_WT)){
    x_WT  <- x_WTs[,,i_lev]

    if(option=="Truth"){
      #Ey2_WT <- var(y_WT[, i_lev])
      #Ey2_WT <- sum(y_WT[, i_lev]^2)/n
      Ey2_WT <- mean(y_WT[, i_lev]^2)
      e_WT <- x_WT - y_WT[, i_lev]
      #e_WT <- x_WT - matrix(y_WT[, i_lev], nrow = n, ncol = p, byrow = TRUE)

      EeeT_WT <- cov(e_WT)
      N_WT <- EeeT_WT / Ey2_WT  # error-to-signal ratio
      #a_WT <- cor(y_WT[, i_lev], x_WT)  # truth
      a_WT <- sapply(1:p, function(i) cor(y_WT[, i_lev], x_WT[,i])*sd(y_WT[, i_lev])/sd(x_WT[,i]))

    } else if(option=="Est") {
      Ey2_WT = var(rowMeans(x_WT))*theta
      #Ey2_WT = max(apply(x_WT, 2, var))*theta
      #Ey2_WT = max(colMeans(x_WT^2))*theta

      ExxT_WT = cov(x_WT)

      SNRest_WT <- SNRest(ExxT_WT, Ey2_WT)

      N_WT <- SNRest_WT$N_est
      if(scale_factor) {
        a_WT <- SNRest_WT$a_est
      } else {
        a_WT <- rep(1, p)
      }
    } else {
      message("No such mode!")
    }

    u_WT <- SNRopt(N_WT, a_WT)

    u_est_WT[i_lev, ] <- u_WT

    y_est_WT[,i_lev] = x_WT%*%u_WT
  }

  y_est <- rowSums(y_est_WT)

  out <- list(weight=u_est_WT,
              merged_WT=y_est_WT,
              merged=y_est)
  return(out)

}

#' snr
#'
#' @param signal signal
#' @param noise noise
#'
#' @return SNR in dB
#' @export
#'
snr <- function(signal, noise) {
  # # Define a function to calculate SNR in dB

  power_signal <- sum(signal^2) / length(signal)
  power_noise <- sum(noise^2) / length(noise)
  SNR <- 10 * log10(power_signal / power_noise)
  return(SNR)
}
