% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snr-freq.R
\name{SNRopt_freq}
\alias{SNRopt_freq}
\title{SNRopt_freq}
\usage{
SNRopt_freq(
  y = NULL,
  x,
  mode = "DWT",
  wf = "haar",
  J,
  pad = "zero",
  boundary = "periodic",
  theta = 0.1,
  option = "Est",
  scale_factor = T
)
}
\arguments{
\item{x}{parent products}

\item{mode}{wavelet decomposition mode}

\item{wf}{wavelet filter}

\item{J}{max level of decomposition level}

\item{boundary}{boundary condition of wavelet transform}

\item{theta}{factor for Ey2 estimation}
}
\value{
A list of weight, merged products.
}
\description{
SNRopt_freq
}
