# SNRBoost

An open-source tool for multi-source data merging. It includes: 

  •	A novel frequency-based merging method resolves spectral differences in multi-source datasets.

  •	SNR optimization minimizes prediction error by recalibrating product weightage in the spectral domain.
  
  •	The proposed approach outperforms conventional methods in combining soil moisture satellite data.

## Requirements
<pre>
Dependencies:
  waveslim, geigen

Suggest:
  readr, WASP
</pre>

## Installation

You can install the package via devtools from [GitHub](https://github.com/) with:

``` r
devtools::install_github("zejiang-unsw/SNRBoost", dependencies = TRUE)
```
