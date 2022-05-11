##KATIE HOEMANN, KU Leuven; last updated 22 April 2022

# set directories
setwd("C:/Users/Katie/Documents/R/")
base_dir <- 'C:/Users/Katie/Documents/R/'

# load required libraries (install if necessary)
library(readr)
library(lme4) 
library(lmerTest)

# read in data: select either 'Study 1 middle', 'Study 1 neg median', or 'Study 1 pos median'
d <- read.csv("IDA_Study 1 pos median.csv")

# set up simple effects coding, following: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/#SIMPLE 
d$StudyS <- d$Study
d$StudyS <- as.factor(d$StudyS)
c <- contr.treatment(3)
coding.matrix <- matrix(rep(1/3,6),ncol=2)
simple.coding <- c-coding.matrix
contrasts(d$StudyS) <- simple.coding

## run regression models including data from all 3 studies
# negative granularity
negGran_IDA <- lm(giniCoefThemeZ ~ StudyS + numPromptsZ + mNegativeZ + negGranZ, data=d) # can include StudyS or negEmotions (or totalEmotions)
summary(negGran_IDA)
confint(negGran_IDA)

# positive granularity
posGran_IDA <- lm(giniCoefThemeZ ~ StudyS + numPromptsZ + mPositiveZ + posGranZ, data=d) # can include StudyS or posEmotions (or totalEmotions)
summary(posGran_IDA)
confint(posGran_IDA)