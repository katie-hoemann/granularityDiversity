##KATIE HOEMANN, KU Leuven; last updated 30 April 2022

# set directories
setwd("C:/Users/Katie/Documents/R/")
base_dir <- 'C:/Users/Katie/Documents/R/'

# load required libraries (install if necessary)
library(readr)
library(meta)
library(DescTools)

# read in data: select either 'neg' or 'pos'
d <- read.csv("meta-analysis_neg.csv")

# transform t to r, maintaining directionality of effect
d$r <- sqrt(d$t^2/(d$t^2 + d$df))
d$r[d$t < 0] <- d$r[d$t < 0]*-1

# meta-analysis using Study 1 median effect size
meta_median <- metacor(data = d, cor = d$r, n = d$n, exclude = 2, sm = "ZCOR")
meta_median

# meta-analysis using Study 1 middle effect size
meta_middle <- metacor(data = d, cor = d$r, n = d$n, exclude = 1, sm = "ZCOR")
meta_middle

# sensitivity analysis (following Hedges & Pigott, 2001, pp. 207-8)
est_r <- .22 # estimated meta-analytic effect size
est_z <- FisherZ(est_r)
est_n <- mean(d$n[c(1,3,4)]) # mean sample size using Study 1 median effect size
v_i <- 1/(est_n-3) # variance of Fisher's z
v. <- v_i/meta_middle$k # variance of weighted mean of Fisher's z
lambda <- (est_z-0)/sqrt(v.)
power <- 1-pnorm(1.96-lambda)+pnorm(-1.96-lambda) # two-tailed power using the estimated effect size (should be close to .80)
power