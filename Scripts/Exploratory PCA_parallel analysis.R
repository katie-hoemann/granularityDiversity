#RYAN BOYD, Lancaster University, last updated on 2020-03-24; Katie Hoemann, KU Leuven, updated 2022-09-20
#----All input should be in .CSV format and MUST be in wide format for this script
#----This script requires the following packages (working as of R version 3.5.1):
# install.packages("psych")
# install.packages("MASS")
# install.packages('paran')
library(psych)
library(MASS)
library("FactoMineR")
library("factoextra")
library(data.table)
library(paran)


#SET PARAMETERS
dataSet <- "ARI" # data set corresponding to folder name
minWords <- 25 # minimum length of texts analyzed
unigrams <- 150 # number unigrams extracted
fileDate <- "2022-10-20" # date MEH file was generated
doPA <- 0 # set to 0 to disable parallel analysis


#LOAD DATA
setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/All texts_longer than ",minWords,"/",unigrams))
DF <- read.csv(paste0(fileDate,"_MEH_DTM_Binary.csv"), fileEncoding = 'UTF-8-BOM')
beginningColumn <- 5 # the first column including variables to be included in the PCA
DF[beginningColumn:length(DF)] <- apply(DF[beginningColumn:length(DF)], 2, as.character)
DF[beginningColumn:length(DF)] <- apply(DF[beginningColumn:length(DF)], 2, as.numeric)
colnames(DF)[1] = 'Filename'


#EXPLORATORY PCA & PARALLEL ANALYSIS
explorePCA <- PCA(DF[beginningColumn:length(DF)],graph=FALSE)
exploreEVs <- get_eigenvalue(explorePCA)
exploreEVs <- as.data.frame(exploreEVs)
EVgt1 <- nrow(subset(exploreEVs, exploreEVs$eigenvalue >= 1))
EVgt2 <- nrow(subset(exploreEVs, exploreEVs$eigenvalue >= 2))
fviz_eig(explorePCA, choice="eigenvalue", geom="line", linecolor = "red", ncp=EVgt1)
if (doPA == 1) paran(x=DF[beginningColumn:length(DF)],iterations=0,graph=TRUE)
  