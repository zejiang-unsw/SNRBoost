#' SNRopt_wav
#'
#' @param x
#' @param mode
#' @param wf
#' @param J
#' @param boundary
#' @param theta
#'
#' @return
#' @export
#' @import waveslim
#' @import WASP
#' @import readr
#'
#' @examples
SNRopt_wav <- function(x, mode="DWT", wf="haar", J, pad="zero", boundary="periodic", theta=0.1) {

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

  # weight estimation
  u_est_WT <- matrix(nrow=J+1, ncol=p)
  y_est_WT <- matrix(nrow=n, ncol=J+1)
  for(i_lev in 1:nrow(u_est_WT)){
    x_WT  <- x_WTs[,,i_lev]

    #Ey2_wT = var(rowMeans(x_WT))*theta
    Ey2_wT = max(apply(x_WT, 2, var))*theta
    ExxT_WT = cov(x_WT)

    SNRest_WT <- SNRest(ExxT_WT, Ey2_wT)

    N_WT <- SNRest_WT$N_est
    a_WT <- SNRest_WT$a_est
    #a_WT <- rep(1, p) # best estimation

    u_WT <- SNRopt(N_WT, a_WT)

    u_est_WT[i_lev, ] <- u_WT

    y_est_WT[,i_lev] = x_WT%*%u_WT
  }

  y_est <- rowMeans(y_est_WT)

  out <- list(weight=u_est_WT,
              merged=y_est_WT)
  return(out)

}
