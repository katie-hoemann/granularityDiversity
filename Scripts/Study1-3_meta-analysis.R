##KATIE HOEMANN, KU Leuven; last updated 30 April 2022

# set directories
setwd("C:/Users/Katie/Documents/R/Text_Analysis/")
base_dir <- 'C:/Users/Katie/Documents/R/Text_Analysis/'

# load required libraries (install if necessary)
library(readr)
library(meta)
library(DescTools)

## NEGATIVE EMOTIONAL GRANULARITY
# read in data
d <- read.csv("Emotional granularity_Experiential diversity_meta_neg.csv")

# transform t to r, maintaining directionality of effect
d$r <- sqrt(d$t^2/(d$t^2 + d$df))
d$r[d$t < 0] <- d$r[d$t < 0]*-1

# meta-analysis using Study 1 median effect size
meta_median_neg <- metacor(data = d, cor = d$r, n = d$n, exclude = 2, sm = "ZCOR")
meta_median_neg

# meta-analysis using Study 1 middle effect size
meta_middle_neg <- metacor(data = d, cor = d$r, n = d$n, exclude = 1, sm = "ZCOR")
meta_middle_neg

# sensitivity analysis (following Hedges & Pigott, 2001, pp. 207-8)
est_r <- .22 # estimated meta-analytic effect size
est_z <- FisherZ(est_r)
est_n <- mean(d$n[c(1,3,4)]) # mean sample size using Study 1 median effect size
v_i <- 1/(est_n-3) # variance of Fisher's z
v. <- v_i/meta_middle_neg$k # variance of weighted mean of Fisher's z
lambda <- (est_z-0)/sqrt(v.)
power <- 1-pnorm(1.96-lambda)+pnorm(-1.96-lambda) # two-tailed power using the estimated effect size (should be close to .80)
power


## POSITIVE EMOTIONAL GRANULARITY
# read in data
d <- read.csv("Emotional granularity_Experiential diversity_meta_pos.csv")

# transform t to r, maintaining directionality of effect
d$r <- sqrt(d$t^2/(d$t^2 + d$df))
d$r[d$t < 0] <- d$r[d$t < 0]*-1

# meta-analysis using Study 1 median effect size
meta_median_pos <- metacor(data = d, cor = d$r, n = d$n, exclude = 2, sm = "ZCOR")
meta_median_pos

# meta-analysis using Study 1 middle effect size
meta_middle_pos <- metacor(data = d, cor = d$r, n = d$n, exclude = 1, sm = "ZCOR")
meta_middle_pos

# sensitivity analysis (following Hedges & Pigott, 2001, pp. 207-8)
est_r <- .22 # estimated meta-analytic effect size
est_z <- FisherZ(est_r)
est_n <- mean(d$n[c(1,3,4)]) # mean sample size using Study 1 median effect size
v_i <- 1/(est_n-3) # variance of Fisher's z
v. <- v_i/meta_middle_pos$k # variance of weighted mean of Fisher's z
lambda <- (est_z-0)/sqrt(v.)
power <- 1-pnorm(1.96-lambda)+pnorm(-1.96-lambda) # two-tailed power using the estimated effect size (should be close to .80)
power

# exclude Study 3 effect size
meta_median_pos_noS3 <- metacor(data = d, cor = d$r, n = d$n, exclude = c(2,4), sm = "ZCOR") # meta-analysis using Study 1 median effect size
meta_median_pos_noS3

meta_middle_pos_noS3 <- metacor(data = d, cor = d$r, n = d$n, exclude = c(1,4), sm = "ZCOR") # meta-analysis using Study 1 middle effect size
meta_middle_pos_noS3