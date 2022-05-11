##KATIE HOEMANN, KU Leuven; last updated 30 April 2022

# set directories
setwd("C:/Users/Katie/Documents/R/")
base_dir <- 'C:/Users/Katie/Documents/R/'

# load required libraries (install if necessary)
library(pwr)

# read in data
d <- read.csv("sensitivity analyses.csv")

# compute denominator degrees of freedom
d$df <- d$n-d$total-1

# run sensitivity analyses and save minimum f2 values
Study1 <- pwr.f2.test(u = d$tested[1], v = d$df[1], f2 = NULL, sig.level = .05, power = .80)
Study2 <- pwr.f2.test(u = d$tested[2], v = d$df[2], f2 = NULL, sig.level = .05, power = .80)
Study3 <- pwr.f2.test(u = d$tested[3], v = d$df[3], f2 = NULL, sig.level = .05, power = .80)
Combined <- pwr.f2.test(u = d$tested[4], v = d$df[4], f2 = NULL, sig.level = .05, power = .80)
d$f2 <- c(Study1$f2, Study2$f2, Study3$f2, Combined$f2)

# transform f2 to standardized beta
d$B <- sqrt(d$f2/(1 + d$f2))