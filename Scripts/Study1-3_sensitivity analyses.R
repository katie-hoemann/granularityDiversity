##KATIE HOEMANN, KU Leuven; last updated 24 September 2022

# set directories
setwd("C:/Users/Katie/Documents/R/Text_Analysis")
base_dir <- 'C:/Users/Katie/Documents/R/Text_Analysis'

# load required libraries (install if necessary)
library(pwr) # documentation: https://cran.r-project.org/web/packages/pwr/pwr.pdf

# read in data
d <- read.csv("Emotional granularity_Experiential diversity_sensitivity.csv")

## TWO-TAILED TESTS
# compute denominator degrees of freedom
d$df <- d$n-d$total-1

# run sensitivity analyses and save minimum f2 values
Study1 <- pwr.f2.test(u = d$tested[1], v = d$df[1], f2 = NULL, sig.level = .05, power = .80)
Study2 <- pwr.f2.test(u = d$tested[2], v = d$df[2], f2 = NULL, sig.level = .05, power = .80)
Study3 <- pwr.f2.test(u = d$tested[3], v = d$df[3], f2 = NULL, sig.level = .05, power = .80)
Combined <- pwr.f2.test(u = d$tested[4], v = d$df[4], f2 = NULL, sig.level = .05, power = .80)
d$f2_2 <- c(Study1$f2, Study2$f2, Study3$f2, Combined$f2)

# transform f2 to standardized beta
d$B_2 <- sqrt(d$f2_2/(1 + d$f2_2))


## ONE-TAILED TESTS
# http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis#Test_on_Individual_Regression_Coefficients_.28t_Test.29
# https://stats.stackexchange.com/questions/494081/is-the-squared-cohens-f-the-same-as-cohens-f2 
# https://stats.stackexchange.com/questions/55236/prove-f-test-is-equal-to-t-test-squared 
# https://www.statisticshowto.com/cohens-f-statistic-definition-formulas/ 
Study1 <- pwr.t.test(n = d$df[1]+1, sig.level = .05, power = .80, type = "one.sample", alternative = "greater")
Study2 <- pwr.t.test(n = d$df[2]+1, sig.level = .05, power = .80, type = "one.sample", alternative = "greater")
Study3 <- pwr.t.test(n = d$df[3]+1, sig.level = .05, power = .80, type = "one.sample", alternative = "greater")
Combined <- pwr.t.test(n = d$df[4]+1, sig.level = .05, power = .80, type = "one.sample", alternative = "greater")
d$d_1 <- c(Study1$d, Study2$d, Study3$d, Combined$d)
d$f2_1 <- d$d_1^2
d$B_1 <- sqrt(d$f2_1/(1 + d$f2_1))