##KATIE HOEMANN, KU Leuven; last updated 24 September 2022

# set directories
setwd("C:/Users/Katie/Documents/R/Text_Analysis")
base_dir <- 'C:/Users/Katie/Documents/R/Text_Analysis'

# load required libraries (install if necessary)
library(readr)
library(lme4) 
library(lmerTest)

## NEGATIVE EMOTIONAL GRANULARITY
# read in data: select either 'Study 1 neg middle' or 'Study 1 neg median'
d <- read.csv("Emotional granularity_Experiential diversity_IDA_Study 1 neg median.csv")

# set up simple effects coding, following: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/#SIMPLE 
d$StudyS <- d$Study
d$StudyS <- as.factor(d$StudyS)
c <- contr.treatment(3)
coding.matrix <- matrix(rep(1/3,6),ncol=2)
simple.coding <- c-coding.matrix
contrasts(d$StudyS) <- simple.coding

# run (hierarchical) regression models including data from all 3 studies
negGran_IDA1 <- lm(negGranZ ~ giniCoefThemeZ + StudyS + StudyS:giniCoefThemeZ, data=d)
negGran_IDA2 <- lm(negGranZ ~ giniCoefThemeZ + StudyS + StudyS:giniCoefThemeZ + mNegativeZ + numPromptsZ, data=d)
summary(negGran_IDA1)
confint(negGran_IDA1)
summary(negGran_IDA2)
confint(negGran_IDA2)


## POSITIVE EMOTIONAL GRANULARITY
# read in data: select either 'Study 1 pos middle' or 'Study 1 pos median'
d <- read.csv("Emotional granularity_Experiential diversity_IDA_Study 1 pos median.csv")

# set up simple effects coding
d$StudyS <- d$Study
d$StudyS <- as.factor(d$StudyS)
c <- contr.treatment(3)
coding.matrix <- matrix(rep(1/3,6),ncol=2)
simple.coding <- c-coding.matrix
contrasts(d$StudyS) <- simple.coding

# run (hierarchical) regression models including data from all 3 studies
posGran_IDA1 <- lm(posGranZ ~ giniCoefThemeZ + StudyS + StudyS:giniCoefThemeZ, data=d)
posGran_IDA2 <- lm(posGranZ ~ giniCoefThemeZ + StudyS + StudyS:giniCoefThemeZ + mPositiveZ + numPromptsZ, data=d)
summary(posGran_IDA1)
confint(posGran_IDA1)
summary(posGran_IDA2)
confint(posGran_IDA2)

# re-run excluding Study 3 data
d <- subset(d, Study<3)
d$StudyS_noS3 <- d$Study
d$StudyS_noS3 <- as.factor(d$StudyS_noS3)
contrasts(d$StudyS_noS3) <- c(-0.5,+0.5)
posGran_IDA2_noS3 <- lm(posGranZ ~ giniCoefThemeZ + StudyS_noS3 + StudyS_noS3:giniCoefThemeZ + mPositiveZ + numPromptsZ, data=d)
summary(posGran_IDA2_noS3)
confint(posGran_IDA2_noS3)













