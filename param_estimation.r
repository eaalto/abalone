# MLE param estimation for abalone
# Using continuous growth function
# Code based off Primer on bbmle2 by NT Hobbs

library(reshape2);
library(ggplot2);
library(bbmle);
# library of utility functions
source("../common/myUtil.r");
source("../common/myData.r");
source("../common/myGraphics.r");
source("../common/myStats.r");
source("../common/myIPM.r");

# Code structure:

############
# CONSTANTS
############

MODEL_NORMAL    = 1;
MODEL_LOGNORMAL = 2;
MODEL_GAMMA     = 3;
MODEL_VONBERT   = 4;

# Run the MLEs for abalone growth
# Compares several possible models
runMLEComparison <- function() {
  growthDataM = loadData();
  #plot(growthDataV, type="l");
  normRes     = runSingleMLE(y=growthDataM);
  
  # display
  print(summary(normRes));
  print(confint(normRes));
#  yearV = 1:length(growthDataV);
  #plot(yearV, growthDataV, xlab="Year", ylab="Length");
#  lines(yearV, 
}

# Runs the MLE process for a single model
runSingleMLE <- function(y, modelType=MODEL_NORMAL) {
  modelRes = mle2(minuslogl=normLLFunc, start=list(L0=4, mean=11, sd=0.9), 
                  data=list(y=y), control=list(maxit=5000));
  return(modelRes);
}

# generates a normally growing individual
generateNormData <- function(L0=5, mean=10, sd=3, numInd=1000, years=10) {
  lM    = matrix(nrow=years, ncol=numInd);
  lM[1,] = L0;
  for (t in 2:years) {
    for (i in 1:numInd) lM[t,i] = lM[t-1,i] + rnorm(1, mean=mean, sd=sd);
  }
  return(lM);
}

meanGrowth <- function(L0, mean, sd, years) {
  lV    = numeric(years);
  lV[1] = L0;
  for (t in 2:years) lV[t] = lV[t-1] + mean;
  return(lV);
}

normLLFunc <- function(y, L0, mean, sd) {
  years = nrow(y);
  meanV = meanGrowth(L0=L0, mean=mean, sd=sd, years);
  LLSum = 0;
  for (i in 1:ncol(y)) {
    LLSum = LLSum - sum(dnorm(y[,i], meanV, sd, log=TRUE));
  }
  return(LLSum);
}

# Load the abalone growth data
loadData <- function() {
  return(generateNormData());
}

# generate test data
generateTestData <- function(modelType=MODEL_NORMAL) {
  
}