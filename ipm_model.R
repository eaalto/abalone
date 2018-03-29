# IPM implementation of abalone model

# General model
# use this for analysis
library(IPMpack)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
# library of utility functions
source("../common/myUtil.r")
source("../common/myData.r")
source("../common/myGraphics.r")
source("../common/myStats.r")
source("../common/myIPM.r")

# TO DO
#  -recruitment patterns
#     -research
#  -OA effects
#  -sensitivity
#     -finer size resolution
#     -finer space resolution

# Naming conventions:
# ALL_CAPS: hard-coded constants (includes defaults, type indicators, like PREY_LARGE, PREY_SMALL, etc.)
# lowerCamelCase: variables, functions
# myXYZ: designates global utility functions from my general library (with some historical exceptions)

# Code structure:
#   CONSTANTS
#   Convenience functions
#   Model:
#     runModel():  Top-level function
#     Initialization functions
#     Presimulation function
#     Main spatial simulation function
#     Plotting

########################################
# ABALONE CONSTANTS 
########################################
###•••••••••••••••••••••  Spatial Arrangement  •••••••••••••••••••••••••••••••••••••••••••••••

# Block and coast layout
BLOCK_WIDTH  = 500                              # Width of the block in m (away from coastline)
BLOCK_LENGTH = 100                              # Length of the block in m (along coastline)
COAST_LENGTH = 15000                	        # Length of the coastline in m
NUM_BLOCKS   = round(COAST_LENGTH/BLOCK_LENGTH) # Number of blocks
BLOCK_AREA   = BLOCK_WIDTH*BLOCK_LENGTH/10000   # Area of each block in ha
TOTAL_AREA   = BLOCK_AREA*NUM_BLOCKS            # Total area in ha

# Dispersal
DISP_PERC    = 99                               # Define percentage of larvae that are retained between +- dispersal distance
DISP_QUANT   = (1-DISP_PERC/100)/2

# Set Class Limits [25 mm width adjusted to 155 mm]
# Actually, the first size class should include 1 year individuals
# Initial size of larvae in Bardos et al. 2006 is 0.3 mm
# Viana et al. measures 2000 postlarvae (unknown age): 909 +- 130 um
# Leighton et al.1981 uses postlarvae (unknown age) of 1-1.5 mm
# Leighton 1974: postlarvae of 30-40 days are 1.7-2 mm SL
# In collectors at Isla Natividad, larvae at settlment were ~500 um
# Using Bardos(2005) and starting with larvae all of 500um, >90% is between 3.5 and 27 mm
# We could simplify and say that recruitment goes all to 1st size class

###•••••••••••••••••••••  Demographic Settings  •••••••••••••••••••••••••••••••••••••••••••••••
MIN_SIZE      = 1     # in mm
MAX_SIZE      = 306#231#306;   # in mm
MIN_LAND_SIZE = 155   # set Minimum Landing Size
CLASS_WIDTH   = 5     # default to 5mm size classes
NUM_CLASSES   = (MAX_SIZE-MIN_SIZE)/CLASS_WIDTH

LARVAL_SURVIVAL  = 0.00309*1.846
# larval survival adjusted by approx. x2 to account for 
# Allee effect in Rosetto calculation (her observed densities were ~0.015/m2
# This is somewhat deprecated
LARVAL_ALLEE_ADJ = 1#1.846;
# approximation from original model
LARVAL_CONST_SIZE = 17.5

############################
# CONVENIENCE FUNCTIONS
############################

# Weight [g] at Length relationship da Shepherd 1998 [with shell]
# to be multiplied by 0.4 to obtain weight [g] without shell [From Natividad Data]
# W = a*L^b
WEIGHT_A  = 2.24 * 10^-5  
WEIGHT_B  = 3.36   
W_NOSHELL = 0.4
wForL <- function(L) { return(WEIGHT_A*(L^(WEIGHT_B))) }     

# Determine proportion dispersing to block at given position
propDispToBlock <- function(position, stdDev) { 
  return(pnorm(BLOCK_LENGTH*(position+1/2),0,stdDev) - 
         pnorm(BLOCK_LENGTH*(position-1/2),0,stdDev))
}

# Calculate dispersal SD for a given random value
getDispSD <- function(randVal) { 
  return(uniroot(function(x)qnorm(DISP_QUANT,mean=0,sd=x)+randVal,lower=0,upper=8000)$root) 
}

#####################
### Constant global vectors
#####################
# Mean class length vector
# Precompute for easy vector operations later
gMeanSizeV = rep(0, NUM_CLASSES)      # initialize below

# Vector giving distribution of new recruits across size classes
gRecruitDistribV = rep(0, NUM_CLASSES) # initialize below

# Size at sexual maturity, Rossetto et al.(2013), not accounting for 1:1 Sex Ratio
# Determines proportion mature at each size
# P(mature) = a / 1+e^(-(size-matSize)/b)
MAT_A    = 1#0.5  DON't CUT IN HALF, modeling all individuals for Allee effect
MAT_B    = 30.20
MAT_SIZE = 135.99
gPropMatureV = c() # initialize below

# Eggs per individual for each size
# Convert size (length) to weight using wForL
# Eggs only from mature individuals (use gPropMatureV)
# Eggs = a*weight*prop
EGGS_A = 3772      # from Tutschulte 1976
gEggsPerIndV = c() # initialize below

# Default catastrophe mortality
DEF_CAT_MORT = 0.75
# Default recruit K param in each block
#DEF_K = 5.43*10^7  # Tuned to unfished density of ~0.4/m2
#DEF_K = 2.59*10^7  # Tuned to unfished density of ~0.2/m2 XXX readj for fixed Allee
DEF_K = 1.29*10^7  # Tuned to unfished density of ~0.2/m2

# Dispersal constants:
DISP_GAMMA_SHAPE = 3
DISP_GAMMA_RATE  = 0.006

# Larval survival constants:
# Gamma and LN tuned to produce ~15x ratio of largest:smallest (after trimming outliers)
# Shepherd 1990
REC_GAMMA_SHAPE  = 2.5
REC_GAMMA_RATE   = 0.5
REC_GB_THRESHOLD = 0.8
REC_GB_GOOD      = 5
REC_GB_BAD       = 0.5
REC_LN_MEANLOG   = 0
REC_LN_SDLOG     = 0.7

# uses the probability of mixed-gender aggregations, based on local density (#/m2)
# Good shorthand
# From Button 2008 thesis, Fig 2.4 and Fig. 5.8
# Changed from 11.33 to 11.6 to set 80% threshold at 0.2/m^2 (well within CI on A value)
AGG_A = 11.6
AGG_B = 1#1.959  # linear fit gives 1.959, but density~0 should have agg=1

# OA EFFECTS
# OA vector indices
oaStrV  = c("Fert.", "Allee", "L.Surv.", "J.Surv", "YA.Surv", "MA.Surv.", "Growth", "MaxSize", "SizeMat", "Fec.", "MultiFJ", "MultiFS", "MultiFG", "MultiJS", "MultiJG")
OA_FERT_SUCCESS_INDEX = 1  # Overall fertilization success
OA_ALLEE_SHIFT_INDEX  = 2  # Effective density multiplier
OA_LARV_SURV_INDEX    = 3  # Larval survival modifier
OA_JUV_SURV_INDEX     = 4  # Juvenile survival index (<50mm)
OA_YA_SURV_INDEX      = 5  # Young adult survival index (>=50mm, <maturity)
OA_MA_SURV_INDEX      = 6  # Mature adult survival index (>=mature size)
OA_GROWTH_RATE_INDEX  = 7  # Growth rate multiplier
OA_MAX_SIZE_INDEX     = 8  # Maximum size multiplier for growth rate
OA_SIZE_MAT_INDEX     = 9  # Size at maturity multiplier
OA_FEC_INDEX          = 10 # Fecundity multiplier
OA_MULTI_FERT_JSURV   = 11 # Combination of fertilization and juvenile survival
OA_MULTI_FERT_SMAT    = 12 # Combination of fertilization and size at maturity
OA_MULTI_FERT_GROW    = 13 # Combination of fertilization and growth
OA_MULTI_JSURV_SMAT   = 14 # Combination of juvenile survival and size at maturity
OA_MULTI_JSURV_GROW   = 15 # Combination of juvenile survival and growth rate   

###••••••••••••••••••••• INITIALIZATION•••••••••••••••••••••••••••••••••••••••••••••••
##### INITIALIZATION CONSTANTS
NUM_REPLICATIONS    = 3     # number of replicates 
MAX_REPLICATE_GROUP = 999999   # break apart into groups
T_MAX               = 300   # simulation length
T_PRESIM            = 300   # Pre-simulation time
H_PRESIM            = 0.0  # No Pre-simulation harvest right now
INIT_ABUN_PRESIM    = 100  # Default initial density for pre-simulation in each size class (#/ha)
INIT_ABUN_NO_PRESIM = 100  # Default initial density if pre-simulation is skipped (#/ha)

#####  Allee types
ALLEE_PROB = 1
ALLEE_LIN  = 2
ALLEE_EXP  = 3

# Recreate using same randomness
RAND_SEED = 123456
# Collapse criterion: <10% of pre-catastrophe equilibrium
COLLAPSE_THRESHOLD = 0.1

# Management strategy labels
MAN_NC   = "NC"
MAN_5YC  = "5YC"
MAN_DC   = "DC"
MAN_MPA  = "MPA"

# Fraction of fMSY for non-spatial
MSY_PROP = 0.67
DISCOUNT_RATE = 0.95

# Parameters for the dynamic management strategy
DYNAMIC_DELAY  = -999
DYNAMIC_THRESH = 0.5

# Default MPA parameters
MPA_SIZE = 5
MPA_PROP = 0.2

# Parameter labels
PARAM_CAT   = "Cat"
PARAM_AGG   = "Agg"
PARAM_K     = "K"
PARAM_DISP  = "Disp"
PARAM_MPA_W = "MPASize"
PARAM_MPA_P = "MPAProp"
PARAM_MSY_P = "MSYProp"
PARAM_CLOSE = "CloseThresh"
PARAM_DISC  = "DiscRate"
PARAM_MULT  = "MultiCat"

initializeGlobalValues <- function(classWidth=CLASS_WIDTH, larvaeConst=FALSE) {
  # initialize NUM_CLASSES
  if (classWidth<=1) MIN_SIZE = 1
  else MIN_SIZE = 5;
  assign("MIN_SIZE", MIN_SIZE, envir=.GlobalEnv);
  assign("CLASS_WIDTH", classWidth, envir=.GlobalEnv);
  NUM_CLASSES   = (MAX_SIZE-MIN_SIZE)/classWidth;
  assign("NUM_CLASSES", NUM_CLASSES, envir=.GlobalEnv);
  # initialize gMeanSizeV
  gMeanSizeV = rep(0, NUM_CLASSES);
  for (i in 1:NUM_CLASSES) { gMeanSizeV[i] = MIN_SIZE + (i-1/2)*CLASS_WIDTH; } # size is average between i and i+1  
  assign("gMeanSizeV", gMeanSizeV, envir=.GlobalEnv);
  # initialize gRecruitDistribV
  gRecruitDistribV = getRecruitDistributionVector(larvaeConst);
  assign("gRecruitDistribV", gRecruitDistribV, envir=.GlobalEnv);
  # initialize gPropMatureV
  gPropMatureV = round(MAT_A/(1+exp(-(gMeanSizeV-MAT_SIZE)/MAT_B)),1)        
  assign("gPropMatureV", gPropMatureV, envir=.GlobalEnv)
  # re-initialize gEggsPerIndV
  # add 0.5 to adjust for 1:1 sex ratio
  gEggsPerIndV = 0.5*EGGS_A*wForL(gMeanSizeV)*gPropMatureV
  assign("gEggsPerIndV", gEggsPerIndV, envir=.GlobalEnv)

  # use the same seed for now
  if (RAND_SEED>0) set.seed(RAND_SEED)
}

#################################################################
####                    SETS TO RUN
#################################################################

combineReplicates <- function(paramStr=PARAM_CAT, prefix="Cont_C1", numV=seq(1,10,by=1), repsPer=100) {
  # XXX assumes equal weight right now
#  if (is.null(weightV)) weightV = rep(1, length(numV))
  
  combinedD = NULL
  for (i in 1:length(numV)) {
    numStr  = numV[i]
    fileStr = combine(combine(combine(prefix, numStr), paramStr), "metrics")
	contD   = loadMatrixFromCSVFile(fileStr, dropFirstCol=TRUE)
	
	# columns: prefix, propExtinct1,2,3; recProp1,2,3; discCatch1,2,3; equilSSB1,2,3; equilCatch1,2,3; meanClosure; sdClosure
	# average all across the files except sdClosure, which must be calculated
	# average must be weighted by non-extinct runs
	# averaging ignores number of replicates since it divides out
	if (is.null(combinedD)) combinedD = contD
	else {
	  curWeight = i-1  # curValues weighted
	  for (j in 1:nrow(combinedD)) {
	    curRowV = as.numeric(combinedD[j,-1])
		newRowV = as.numeric(contD[j,-1])
		updatedRowV = rep(0, length(curRowV))
		# three scenarios
		for (k in 1:3) {
		  curExt = curRowV[k]
		  newExt = newRowV[k]
		  
		  # averaging scenario columns with some special handling
		  # always update proportion extinct with weighting
		  updatedRowV[k]   = (curExt*curWeight + newExt) / (curWeight + 1)  # propExtinct
		  
		  if (curExt==1) {
		    # just copy new values (doesn't matter if both extinct)
		    updatedRowV[3+k] = newRowV[3+k] # recProp
			updatedRowV[6+k] = newRowV[6+k] # discCatch
			updatedRowV[9+k] = newRowV[9+k] # equilSSB
			updatedRowV[12+k] = newRowV[12+k] # equilCatch
			updatedRowV[15+k] = newRowV[15+k] # varSSB
			updatedRowV[18+k] = newRowV[18+k] # varCatch
			updatedRowV[21+k] = newRowV[21+k] # size
			updatedRowV[24+k] = newRowV[24+k] # sizeSSB
			updatedRowV[27+k] = newRowV[27+k] # meanFec
			updatedRowV[30+k] = newRowV[30+k] # totalFec
			if (k==2) {
			  updatedRowV[34] = newRowV[34] # meanClosure
			  updatedRowV[35] = newRowV[35] # sdClosure
			}
		  } else if (newExt==1) {
		    # just copy old values
		    updatedRowV[3+k] = curRowV[3+k] # recProp
			updatedRowV[6+k] = curRowV[6+k] # discCatch
			updatedRowV[9+k] = curRowV[9+k] # equilSSB
			updatedRowV[12+k] = curRowV[12+k] # equilCatch
			updatedRowV[15+k] = newRowV[15+k] # varSSB
			updatedRowV[18+k] = newRowV[18+k] # varCatch
			updatedRowV[21+k] = newRowV[21+k] # size
			updatedRowV[24+k] = newRowV[24+k] # sizeSSB
			updatedRowV[27+k] = newRowV[27+k] # meanFec
			updatedRowV[30+k] = newRowV[30+k] # totalFec
			if (k==2) {
			  updatedRowV[34] = curRowV[34] # meanClosure
			  updatedRowV[35] = curRowV[35] # sdClosure
			}
		  } else {
  		    # values to weight averaging of survival-only numbers
		    curSurvWeight = (1-curExt)*curWeight
		    newSurvWeight = (1-newExt)
		    
			# recovery proportion averaged by survivor proportion
		    updatedRowV[3+k] = (curRowV[3+k]*curSurvWeight + newRowV[3+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[6+k] = (curRowV[6+k]*curSurvWeight + newRowV[6+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[9+k] = (curRowV[9+k]*curSurvWeight + newRowV[9+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[12+k] = (curRowV[12+k]*curSurvWeight + newRowV[12+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[15+k] = (curRowV[15+k]*curSurvWeight + newRowV[15+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[18+k] = (curRowV[18+k]*curSurvWeight + newRowV[18+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[21+k] = (curRowV[21+k]*curSurvWeight + newRowV[21+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[24+k] = (curRowV[24+k]*curSurvWeight + newRowV[24+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[27+k] = (curRowV[27+k]*curSurvWeight + newRowV[27+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    updatedRowV[30+k] = (curRowV[30+k]*curSurvWeight + newRowV[30+k]*newSurvWeight) / (curSurvWeight + newSurvWeight)
		    if (k==2) {
			  # meanClosure
			  curMean = curRowV[34]
			  newMean = newRowV[34]

			  updatedRowV[34] = (curMean*curSurvWeight + newMean*newSurvWeight) / (curSurvWeight + newSurvWeight) 
			  # sdClosure
			  # Can't be averaged, more complex
			  # Multiply by repsPer to restore absolute sample size
			  n_x   = curSurvWeight*repsPer # old "sample size"
			  n_y   = newSurvWeight*repsPer # new "sample size"
			  var_x = curRowV[35]^2
			  var_y = newRowV[35]^2
			  
              # the formula: 			  
			  numer = (n_x^2)*var_x + (n_y^2)*var_y - n_y*var_x - n_y*var_y - n_x*var_x - n_x*var_y + 
			          n_y*n_x*var_x + n_y*n_x*var_y + n_y*n_x*(curMean-newMean)^2 
			  denom = (n_y + n_x - 1)*(n_x+n_y)

              newVar = numer / denom
			  updatedRowV[35] = sqrt(newVar)
			}
		  }
		
	    }
        
		# update row
  	    combinedD[j,-1] = updatedRowV
	  }
    }
  }
  saveMatrixToCSVFile(combine(combine("Combined", combine(prefix, paramStr)), "metrics"), combinedD)
}

runTest <- function(numReplic=2, maxT=300) {
#  runContinuum(numReplic=numReplic, contV=c(0.75), paramStr=PARAM_CAT, isTest=TRUE)
#  runContinuum(numReplic=numReplic, contV=c(0.9), paramStr=PARAM_AGG, isTest=TRUE)
#  runContinuum(numReplic=numReplic, contV=c(0.9), paramStr=PARAM_K, isTest=TRUE)
  kMultV = seq(0.2, 1.0, by=0.1)
  statM  = matrix(nrow=length(kMultV)+1, ncol=4)
  statM[1,] = c("", "Prop extinct", "Recovery", "Equilib.")
  statM[,1] = c("", kMultV)
  
  for (i in 1:length(kMultV)) {
    k    = kMultV[i]
    newH = getMSYValue(paramStr=PARAM_K, paramVal=k, manType=MAN_NC, propMSY=0.8)
	newK = DEF_K*k
	
	resL = runDensityTest(numReplic=numReplic, hV=c(newH), manV=c(MAN_NC), blockK=newK, showPlots=FALSE)
	statM[i+1,-1] = c(resL$propExtinct1, resL$recProp1, resL$equilSSB1)
  }
  
  print(statM)

}

runContinuum <- function(numReplic=100, maxT=300, contV=NULL, calcOnly=FALSE, paramStr=PARAM_CAT, isTest=FALSE, prefix="",
                         run5YCOnly=FALSE, defCat=DEF_CAT_MORT, setRand=RAND_SEED, manV=NULL, equalEffort=TRUE, runParallel=FALSE) {
  # Allow rerunning with diff rand
  if (setRand>0) set.seed(setRand)

  dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma")
  randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=REC_LN_SDLOG)

  if (is.null(manV)) manV = c(MAN_NC, MAN_DC, MAN_MPA)

  if (run5YCOnly) manV = c(MAN_5YC)
  if (paramStr==PARAM_CAT) {
    if (is.null(contV)) contV = seq(0.5, 0.9, by=0.05) # catMort
  } else if (paramStr==PARAM_AGG) {
    # Minimum 0.5, or population extinct
    if (is.null(contV)) contV = seq(0.5, 2.0, by=0.1) # aggA mult
  } else if (paramStr==PARAM_K) {
    if (is.null(contV)) contV = seq(0.4, 2.0, by=0.1) # blockK mult 
  } else if (paramStr==PARAM_DISP) {
    if (is.null(contV)) contV = seq(0.25, 2.0, by=0.25) # dispersal mult 
    manV = c(MAN_MPA)  # run MPA only
  } else if (paramStr==PARAM_MPA_W) {
    if (is.null(contV)) contV = c(1, 2, 5, 10, 15) # MPA width 
    manV = c(MAN_MPA)  # run MPA only
  } else if (paramStr==PARAM_MPA_P) {
    if (is.null(contV)) contV = seq(0.1, 0.8, by=0.1) # MPA proportion 
    manV = c(MAN_MPA)  # run MPA only
  } else if (paramStr==PARAM_MSY_P) {
    if (is.null(contV)) contV = c(0.5, 0.6, 0.67, 0.75, 0.8, 0.9, 1.0) # fMSY proportion 
    manV = c(MAN_NC, MAN_DC) # run non-spatial only
	if (run5YCOnly) manV = c(MAN_5YC)
  } else if (paramStr==PARAM_CLOSE) {
    if (is.null(contV)) contV = seq(0, 0.8, by=0.1) # closure threshold
    manV = c(MAN_DC)  # run DC only
  } else if (paramStr==PARAM_DISC) {
    if (is.null(contV)) contV = c(1,0.975,0.95,0.925,0.9) # discount rate
  } else if (paramStr==PARAM_MULT) {
    if (is.null(contV)) contV = c(15, 10, 5.2, 5.3)
    defCat = 0.5  # always use lightest cat
  } 

  # set the defaults
  contMort    = defCat  # Allows this one to be set independently
  contAgg     = AGG_A
  contK       = DEF_K
  dispMult    = 1.0
  mpaSize     = MPA_SIZE
  mpaProp     = MPA_PROP
  msyProp     = MSY_PROP 
  closeThresh = DYNAMIC_THRESH
  discRate    = DISCOUNT_RATE
  numCat      = 1
  catInt      = 10

  if (runParallel) {
    # let's run the replications in parallel
    theCluster = makeCluster(detectCores()/2)
    registerDoParallel(theCluster)
  
    metricD = foreach (i=1:length(contV), .combine=rbind) %dopar% {
      source("ipm_model.R")
      # loop through parameter values
      contVal = contV[i]    
	
	  # adjust changing parameter
	  if (paramStr==PARAM_CAT) contMort = contVal  
	  else if (paramStr==PARAM_AGG) contAgg = AGG_A * contVal # multiple
	  else if (paramStr==PARAM_K) contK  = DEF_K * contVal # multiple
	  else if (paramStr==PARAM_DISP) dispMult = contVal
      else if (paramStr==PARAM_MPA_W) mpaSize = contVal
	  else if (paramStr==PARAM_MPA_P) mpaProp = contVal
	  else if (paramStr==PARAM_MSY_P) msyProp = contVal
	  else if (paramStr==PARAM_CLOSE) closeThresh = contVal
	  else if (paramStr==PARAM_DISC) discRate = contVal
      else if (paramStr==PARAM_MULT) {
        if (contVal==5.2) {
          numCat = 2
          catInt = 5
        } else if (contVal==5.3) {
          numCat = 3
          catInt = 5
        } else {
          numCat = 2
          catInt = contVal
        }
      }
      
      # Fetch appropriate MSY harvest values (apply new msyProp here)
	  if ((paramStr==PARAM_CAT)||(paramStr==PARAM_MSY_P)||(paramStr==PARAM_CLOSE)||(paramStr==PARAM_DISC)||(paramStr==PARAM_MULT)) {
	    # First group: no parameter dependency for either
	    ncStr  = "Base"
	    mpaStr = "Base"
      } else if ((paramStr==PARAM_DISP)||(paramStr==PARAM_MPA_W)||(paramStr==PARAM_MPA_P)) {
        # Second group: no dependency for NC/DC
  	    ncStr = "Base"
	    mpaStr = paramStr
	  } else {
	    # Third group: dependency for both
	    ncStr  = paramStr
	    mpaStr = paramStr
      }

      runPrefix = combine(prefix, combine(paramStr, contVal))
 	
      if (calcOnly) {
	    pL = showDensityMetrics(NULL, prefix=runPrefix, manV=manV, printMetrics=FALSE)
#	    if (is.null(metricD)) metricD = as.data.frame(pL)
#  	    else metricD = rbind(metricD, as.data.frame(pL))
	  } else {
  	    # keep constant effort, so 100% at 2/3 MSY for NC = 2/3 / (1-MPA%)
	    ncH  = getMSYValue(paramStr=ncStr, paramVal=contVal, manType=MAN_NC, propMSY=msyProp)
	    # adjusts automatically if mpaProp is shifting
	    if (equalEffort) mpaH = ncH / (1 - mpaProp)
	    else mpaH = getMSYValue(paramStr=mpaStr, paramVal=contVal, manType=MAN_MPA, propMSY=1.0)
	    hV   = c(ncH, mpaH)
	
	    if (paramStr==PARAM_DISP) {
          # recalculate dispersal here as a multiple
          dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma", shape=DISP_GAMMA_SHAPE*dispMult)
        }	
	  
        # set other changing parameters here
	    pL = runDensityTest(catMort=contMort, aggA=contAgg, blockK=contK, mpaProp=mpaProp, mpaSize=mpaSize, closeThresh=closeThresh,
	                        hV=hV, prefix=runPrefix, numReplic=numReplic, maxT=maxT, manV=manV, numCat=numCat, catInt=catInt, 
	                        dispSDM=dispSDM, randRecM=randRecM, discRate=discRate, showPlots=FALSE, showMetrics=FALSE)
	  }
#	  if (is.null(metricD)) metricD = as.data.frame(pL)
#	  else metricD = rbind(metricD, as.data.frame(pL))
      contD = as.data.frame(pL)
    }

    # shut down the cluster
    stopCluster(theCluster)
  } else { # not parallel
    metricD = NULL
	for (i in 1:length(contV)) {
      # loop through parameter values
      contVal = contV[i]    
	
	  # adjust changing parameter
	  if (paramStr==PARAM_CAT) contMort = contVal  
	  else if (paramStr==PARAM_AGG) contAgg = AGG_A * contVal # multiple
	  else if (paramStr==PARAM_K) contK  = DEF_K * contVal # multiple
	  else if (paramStr==PARAM_DISP) dispMult = contVal
      else if (paramStr==PARAM_MPA_W) mpaSize = contVal
	  else if (paramStr==PARAM_MPA_P) mpaProp = contVal
	  else if (paramStr==PARAM_MSY_P) msyProp = contVal
	  else if (paramStr==PARAM_CLOSE) closeThresh = contVal
	  else if (paramStr==PARAM_DISC) discRate = contVal
      else if (paramStr==PARAM_MULT) {
        if (contVal==5.2) {
          numCat = 2
          catInt = 5
        } else if (contVal==5.3) {
          numCat = 3
          catInt = 5
        } else {
          numCat = 2
          catInt = contVal
        }
      }
      
      # Fetch appropriate MSY harvest values (apply new msyProp here)
	  if ((paramStr==PARAM_CAT)||(paramStr==PARAM_MSY_P)||(paramStr==PARAM_CLOSE)||(paramStr==PARAM_DISC)||(paramStr==PARAM_MULT)) {
	    # First group: no parameter dependency for either
	    ncStr  = "Base"
	    mpaStr = "Base"
      } else if ((paramStr==PARAM_DISP)||(paramStr==PARAM_MPA_W)||(paramStr==PARAM_MPA_P)) {
        # Second group: no dependency for NC/DC
  	    ncStr = "Base"
	    mpaStr = paramStr
	  } else {
	    # Third group: dependency for both
	    ncStr  = paramStr
	    mpaStr = paramStr
      }

      runPrefix = combine(prefix, combine(paramStr, contVal))
 	
      if (calcOnly) {
	    pL = showDensityMetrics(NULL, prefix=runPrefix, manV=manV, printMetrics=FALSE)
#	    if (is.null(metricD)) metricD = as.data.frame(pL)
#  	    else metricD = rbind(metricD, as.data.frame(pL))
	  } else {
  	    # keep constant effort, so 100% at 2/3 MSY for NC = 2/3 / (1-MPA%)
	    ncH  = getMSYValue(paramStr=ncStr, paramVal=contVal, manType=MAN_NC, propMSY=msyProp)
	    # adjusts automatically if mpaProp is shifting
	    if (equalEffort) mpaH = ncH / (1 - mpaProp)
	    else mpaH = getMSYValue(paramStr=mpaStr, paramVal=contVal, manType=MAN_MPA, propMSY=1.0)
	    hV   = c(ncH, mpaH)
	
	    if (paramStr==PARAM_DISP) {
          # recalculate dispersal here as a multiple
          dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma", shape=DISP_GAMMA_SHAPE*dispMult)
        }	
	  
        # set other changing parameters here
	    pL = runDensityTest(catMort=contMort, aggA=contAgg, blockK=contK, mpaProp=mpaProp, mpaSize=mpaSize, closeThresh=closeThresh,
	                        hV=hV, prefix=runPrefix, numReplic=numReplic, maxT=maxT, manV=manV, numCat=numCat, catInt=catInt, 
	                        dispSDM=dispSDM, randRecM=randRecM, discRate=discRate, showPlots=FALSE, showMetrics=FALSE)
	  }
	  if (is.null(metricD)) metricD = as.data.frame(pL)
	  else metricD = rbind(metricD, as.data.frame(pL))
    }  
  }
  
  if (isTest) paramStr = combine("Test", paramStr)
  else paramStr = combine(prefix, paramStr)
  saveMatrixToCSVFile(combine(combine("Cont", paramStr), "metrics"), metricD)
}

runAllSizeRuns <- function(numReps=10) {
  # To test: closest to 50% survival
#  surv1 = runMatchingOASizeRuns(OA_JUV_SURV_INDEX, 0.96, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv2 = runMatchingOASizeRuns(OA_YA_SURV_INDEX, 0.96, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv3 = runMatchingOASizeRuns(OA_MA_SURV_INDEX, 0.96, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv4 = runMatchingOASizeRuns(OA_GROWTH_RATE_INDEX, 0.96, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv5 = runMatchingOASizeRuns(OA_MAX_SIZE_INDEX, 0.96, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv6 = runMatchingOASizeRuns(OA_SIZE_MAT_INDEX, 1.04, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  surv7 = runMatchingOASizeRuns(OA_FEC_INDEX, 0.97, numReps=numReps, plotResults=FALSE, randSeed=RAND_SEED)
#  plotSizeDiffs(OA_JUV_SURV_INDEX, 0.96, surv1)
#  plotSizeDiffs(OA_YA_SURV_INDEX, 0.96, surv2)
#  plotSizeDiffs(OA_MA_SURV_INDEX, 0.96, surv3)
#  plotSizeDiffs(OA_GROWTH_RATE_INDEX, 0.96, surv4)
#  plotSizeDiffs(OA_MAX_SIZE_INDEX, 0.96, surv5)
#  plotSizeDiffs(OA_SIZE_MAT_INDEX, 1.04, surv6)
#  plotSizeDiffs(OA_FEC_INDEX, 0.97, surv7)

#  runAveragedOASizes(OA_GROWTH_RATE_INDEX, 1, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_JUV_SURV_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_YA_SURV_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_MA_SURV_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_GROWTH_RATE_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_MAX_SIZE_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_SIZE_MAT_INDEX, 1.05, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_FEC_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
#  runAveragedOASizes(OA_ALLEE_SHIFT_INDEX, 0.95, numReps=numReps, randSeed=RAND_SEED)
  
  plotMeanSizes(c(OA_GROWTH_RATE_INDEX, OA_JUV_SURV_INDEX, OA_YA_SURV_INDEX, OA_MA_SURV_INDEX, OA_GROWTH_RATE_INDEX, OA_MAX_SIZE_INDEX, OA_SIZE_MAT_INDEX, OA_FEC_INDEX, OA_ALLEE_SHIFT_INDEX), 
                c(1, 0.95, 0.95, 0.95, 0.95, 0.95, 1.05, 0.95, 0.95), survOnly=FALSE)

}

runAveragedOASizes <- function(oaIndex=OA_GROWTH_RATE_INDEX, paramVal=1, numReps=10, prefix="Sizes", maxT=300, addCat=FALSE, randSeed=RAND_SEED) {
  set.seed(randSeed)
  dispSDM  = getRandomDispSDMatrix(maxT, numReps, rType="gamma")
  randRecM = getRandomRecruitMatrix(maxT, numReps, rType="log-normal", shape=REC_LN_SDLOG)

  manV = c(MAN_NC)#, MAN_MPA)

  metricD  = NULL
  paramStr = oaStrV[oaIndex]
  myP("Running OA size models for", paramStr)

  hV = c(getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_NC, propMSY=MSY_PROP), 
	     getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_MPA, propMSY=1.0))

  # Array for size structures
  oaSizeM = matrix(nrow=(length(gMeanSizeV)+1), ncol=numReps)
  
  thisPrefix = combine(prefix, combine(paramStr, paramVal))
  # set oaEffectV
  oaEffectV  = rep(1,length(oaStrV))
  oaEffectV[oaIndex] = paramVal
  for (i in 1:numReps) {
    # loop through all available replications
    myP("Running test replicate #", i)
	
    # XXX some issues with single row matrices I'd rather not bother fixing
	# just duplicating instead
	singleDispM = cbind(dispSDM[,i], dispSDM[,i])
	singleRecM  = cbind(randRecM[,i], randRecM[,i])

	pL = runDensityTest(hV=hV, prefix=thisPrefix, numReplic=2, maxT=maxT, manV=manV, oaEffectV=oaEffectV, addCat=addCat,
	                      catMort=0.5, dispSDM=singleDispM, randRecM=singleRecM, showPlots=FALSE, showMetrics=FALSE)

	survRun = (pL$propExtinct1==0)
    # load size files and store
    fileStr = combine(thisPrefix, "NC_size")
    # load and drop first column
    sizeM = loadMatrixFromCSVFile(fileStr, dropFirstCol=TRUE)
	
	# average the central portion (avoids issues with late near-extinction)
	sizeV = apply(sizeM[,50:(maxT-50)], 1, sum)
	oaSizeM[1,i]  = ifelse(survRun, 1, 0)
    oaSizeM[-1,i] = sizeV
  }
  
  rownames(oaSizeM) <- c("Surv", gMeanSizeV)
  colnames(oaSizeM) <- 1:numReps

  saveMatrixToCSVFile(combine(thisPrefix, "repEquilSize"), oaSizeM)
}

plotMeanSizes <- function(oaIndexV=c(OA_GROWTH_RATE_INDEX), paramV=c(1), dropYOY=TRUE, survOnly=TRUE, prefix="Sizes") {
  allD = NULL
  for (i in 1:length(oaIndexV)) {
    paramD = plotMeanSize(oaIndexV[i], paramV[i], dropYOY, survOnly, prefix, dataOnly=TRUE)
	if (is.null(allD)) allD = paramD
	else allD = rbind(allD, paramD)
  }
  
  thePlot = ggplot(data=allD, aes(x=size, y=prop)) + geom_line(aes(color=param)) +
            labs(title="Comparison", x="Size", y="Proportion")
  myPPlot(thePlot)
  
  oaV   = unique(allD$param)
  baseV = allD[allD$param=="Baseline",]$prop
  allD$diff = allD$prop - baseV # baseV is repeated

  thePlot = ggplot(data=allD, aes(x=size, y=diff)) + geom_line(aes(color=param)) +
            labs(title="Difference", x="Size", y="Difference from baseline")
  myPPlot(thePlot)

}

plotMeanSize <- function(oaIndex=OA_GROWTH_RATE_INDEX, paramVal=-1, dropYOY=TRUE, survOnly=TRUE, prefix="Sizes", dataOnly=FALSE) {
  thisPrefix = combine(prefix, combine(oaStrV[oaIndex], paramVal))
  fileStr    = combine(thisPrefix, "repEquilSize")
  meanSizeM  = as.matrix(loadMatrixFromCSVFile(fileStr))
  numReps    = ncol(meanSizeM)-1

  # get survival vector, then drop
  survV     = as.numeric(meanSizeM[1,-1])
  meanSizeM = matrix(as.numeric(meanSizeM[-1,]), nrow=(nrow(meanSizeM)-1), ncol=ncol(meanSizeM))
  
  # convert to tidy data
  # drop name column
  rownames(meanSizeM) <- meanSizeM[,1]
  meanSizeM = meanSizeM[,-1]
  colnames(meanSizeM) <- 1:numReps
  
  # use only surviving runs
  if (survOnly) {
    meanSizeM = meanSizeM[,survV]
    numReps   = sum(survV)
  }
  # remove the first 2 size classes which have much larger numbers
  if (dropYOY) meanSizeM   = meanSizeM[-1:-2,]
  
  # Convert to proportions for each rep to weight equally
  for (i in 1:numReps) {
    meanSizeM[,i] = meanSizeM[,i] / sum(meanSizeM[,i])
  }
  
  # average across reps
  propV = apply(meanSizeM, 1, mean) 
  sizeV = as.numeric(rownames(meanSizeM))
  paramStr  = ifelse(paramVal==1, "Baseline", combine(oaStrV[oaIndex], paramVal))
  meanSizeD = data.frame(size=sizeV, prop=propV, param=paramStr, stringsAsFactors=FALSE)
  if (dataOnly) return(meanSizeD)
 
  thePlot = ggplot(data=meanSizeD, aes(x=size, y=prop)) + geom_line() +
            labs(title=thisPrefix, x="Size", y="Proportion")
  myPPlot(thePlot)
}

runMatchingOASizeRuns <- function(oaIndex=OA_GROWTH_RATE_INDEX, paramVal=-1, numReps=10, keepSurvOnly=TRUE, prefix="SizeDiff", maxT=300, addCat=FALSE, dropYOY=TRUE, plotResults=TRUE, randSeed=RAND_SEED) {
  set.seed(randSeed)
  dispSDM  = getRandomDispSDMatrix(maxT, numReps, rType="gamma")
  randRecM = getRandomRecruitMatrix(maxT, numReps, rType="log-normal", shape=REC_LN_SDLOG)

  manV = c(MAN_NC)#, MAN_MPA)

  metricD  = NULL
  paramStr = oaStrV[oaIndex]
  myP("Running OA size models for", paramStr)

  hV = c(getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_NC, propMSY=MSY_PROP), 
	     getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_MPA, propMSY=1.0))

  multiIndex = FALSE
		 
  if (paramVal==-1) {
    if (oaIndex==OA_SIZE_MAT_INDEX) oaValV = c(1, 1.15) 
    else if ((oaIndex==OA_GROWTH_RATE_INDEX)||(oaIndex==OA_MA_SURV_INDEX)||(oaIndex==OA_MAX_SIZE_INDEX)) oaValV = c(1, 0.95)
    else if ((oaIndex>=OA_MULTI_FERT_JSURV)) {  # XXX special case for multi effects
      multiIndex = TRUE
	  if (oaIndex==OA_MULTI_FERT_JSURV) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_JUV_SURV_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_SMAT) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 1.05)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_GROW) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_SMAT) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 1.05)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_GROW) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  }
	  oaValV  = oaValV1
	  oaIndex = index1
    } else oaValV = c(1, 0.9)
  } else oaValV = c(1, paramVal)
  
  # loop through parameter values
  # reverse these to find a surviving run
  oaValV = rev(oaValV)
  if (multiIndex) oaValV2 = rev(oaValV2)

  # Arrays for size structures
  oaSizeA   = array(dim=c(length(gMeanSizeV), maxT, numReps))
  baseSizeA = array(dim=c(length(gMeanSizeV), maxT, numReps))
  
  survRuns = 0
  for (i in 1:numReps) {
    # loop through all available replications
	useRep = FALSE
    myP("Running test replicate #", i)
	prefix1 = ""
	prefix2 = ""
    for (j in 1:length(oaValV)) {
      oaVal   = oaValV[j]
	  # XXX sloppy implementation here
	  thisPrefix = combine(prefix, combine(paramStr, oaVal))
      if (j==1) prefix1 = thisPrefix
      else if (j==2) prefix2 = thisPrefix
      # set oaEffectV
      oaEffectV = rep(1,length(oaStrV))
      oaEffectV[oaIndex] = oaVal
      if (multiIndex) {
        # Add the second effect
        oaEffectV[index2] = oaValV2[j]	  
 	  }
	
	  # XXX some issues with single row matrices I'd rather not bother fixing
	  # just duplicating instead
	  singleDispM = cbind(dispSDM[,i], dispSDM[,i])
	  singleRecM  = cbind(randRecM[,i], randRecM[,i])

	  pL = runDensityTest(hV=hV, prefix=thisPrefix, numReplic=2, maxT=maxT, manV=manV, oaEffectV=oaEffectV, addCat=addCat,
	                      catMort=0.5, dispSDM=singleDispM, randRecM=singleRecM, showPlots=FALSE, showMetrics=FALSE)

	  if (keepSurvOnly && (j==1) && (pL$propExtinct1==1)) break
      else useRep = TRUE
	}

	if (useRep) {
	  # load size files and store
      fileStr = combine(prefix1, "NC_size")
      # load and drop first column
      compM   = loadMatrixFromCSVFile(fileStr, dropFirstCol=TRUE)
      fileStr = combine(prefix2, "NC_size")
      baseM   = loadMatrixFromCSVFile(fileStr, dropFirstCol=TRUE)
	  
		  # store in appropriate array
      survRuns = survRuns + 1
	  oaSizeA[,,survRuns]   = as.matrix(compM)
	  baseSizeA[,,survRuns] = as.matrix(baseM)
	}
  }
  
  oaSizeA   = oaSizeA[,,1:survRuns]
  baseSizeA = baseSizeA[,,1:survRuns]

  sizeV = gMeanSizeV
  if (dropYOY) {
    # remove the first 2 size classes which have much larger numbers
    oaSizeA   = oaSizeA[-1:-2,,]
    baseSizeA = baseSizeA[-1:-2,,]
	sizeV     = sizeV[-1:-2]
  }
  
  # Convert to proportions
  for (i in 1:survRuns) {
    for (j in 1:maxT) {
      oaSizeA[,j,i]   = oaSizeA[,j,i] / sum(oaSizeA[,j,i])
      baseSizeA[,j,i] = baseSizeA[,j,i] / sum(baseSizeA[,j,i])
    }
  }
  
  # Mean size difference
  diffSizeA = oaSizeA - baseSizeA
  meanDiffM = apply(diffSizeA, c(1,2), mean)
  
  # 3 major classes: J <50mm, YA <136mm, MA >136
  class3M     = matrix(nrow=3, ncol=maxT)
  class3M[1,] = apply(meanDiffM[1:7,], 2, sum)
  class3M[2,] = apply(meanDiffM[8:25,], 2, sum)
  class3M[3,] = apply(meanDiffM[26:length(sizeV),], 2, sum)

  rownames(meanDiffM) <- sizeV
  colnames(meanDiffM) <- 1:maxT
  rownames(class3M) <- c("Juv.", "YA", "MA")
  colnames(class3M) <- 1:maxT

  paramStr = combine(combine(prefix, combine(paramStr, paramVal)), paste("n", survRuns, sep=""))
  saveMatrixToCSVFile(combine(paramStr, "meanDiff"), meanDiffM)
  saveMatrixToCSVFile(combine(paramStr, "3Diff"), class3M)

  if (plotResults) plotSizeDiffs(oaIndex, paramVal, survRuns)
  return(survRuns)
}

plotSizeDiffs <- function(oaIndex=OA_GROWTH_RATE_INDEX, paramVal=-1, survRuns=100, prefix="SizeDiff") {
  thisPrefix = combine(combine(prefix, combine(oaStrV[oaIndex], paramVal)), paste("n", survRuns, sep=""))
  fileStr   = combine(thisPrefix, "meanDiff")
  meanDiffM = as.matrix(loadMatrixFromCSVFile(fileStr))
  fileStr   = combine(thisPrefix, "3Diff")
  class3M   = loadMatrixFromCSVFile(fileStr)
  maxT = ncol(meanDiffM)-1

  # convert to tidy data
  # drop name column
  rownames(meanDiffM) <- meanDiffM[,1]
  meanDiffM = meanDiffM[,-1]
  colnames(meanDiffM) <- 1:maxT
  rownames(class3M) <- c("Juv.", "YA", "MA")
  class3M = as.matrix(class3M[,-1])
  colnames(class3M) <- 1:maxT
 
  meanDiffD = melt(meanDiffM, varnames=c("size", "time"), value.name="diff")
  class3D   = melt(class3M, varnames=c("size", "time"), value.name="diff")
  
  thePlot = ggplot(data=meanDiffD, aes(x=time, y=diff)) + geom_line(aes(color=size)) +
            labs(title=paste(thisPrefix, "All classes", sep=":"), x="Time", y="Diff. from base")
  myPPlot(thePlot)

  thePlot = ggplot(data=class3D, aes(x=time, y=diff)) + geom_line(aes(color=size)) +
            labs(title=paste(thisPrefix, "Juv, YA, MA", sep=":"), x="Time", y="Diff. from base")
  myPPlot(thePlot)
}

runSingleOAInstance <- function(oaIndex=OA_GROWTH_RATE_INDEX, paramVal=-1, startRep=1, runNoOA=FALSE, prefix="Single", maxT=300, repNum=1, addCat=FALSE) {
  dispSDM  = getRandomDispSDMatrix(maxT, 100, rType="gamma")
  randRecM = getRandomRecruitMatrix(maxT, 100, rType="log-normal", shape=REC_LN_SDLOG)

  manV = c(MAN_NC)#, MAN_MPA)

  metricD  = NULL
  paramStr = oaStrV[oaIndex]
  myP("Running OA models for", paramStr)

  hV = c(getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_NC, propMSY=MSY_PROP), 
	     getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_MPA, propMSY=1.0))

  multiIndex = FALSE
		 
  if (paramVal==-1) {
    if (oaIndex==OA_SIZE_MAT_INDEX) oaValV = c(1, 1.15) 
    else if ((oaIndex==OA_GROWTH_RATE_INDEX)||(oaIndex==OA_MA_SURV_INDEX)||(oaIndex==OA_MAX_SIZE_INDEX)) oaValV = c(1, 0.95)
    else if ((oaIndex>=OA_MULTI_FERT_JSURV)) {  # special case for multi effects
      multiIndex = TRUE
	  if (oaIndex==OA_MULTI_FERT_JSURV) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_JUV_SURV_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_SMAT) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 1.05)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_GROW) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_SMAT) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 1.05)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_GROW) {
	    oaValV1 = c(1, 0.95)
	    oaValV2 = c(1, 0.95)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  }
	  oaValV  = oaValV1
	  oaIndex = index1
    } else oaValV = c(1, 0.9)
  } else oaValV = c(1, paramVal)
  
  # loop through parameter values
  # reverse these to find a surviving run
  oaValV = rev(oaValV)
  if (multiIndex) oaValV2 = rev(oaValV2)
  
  repNum = startRep
  for (i in 1:length(oaValV)) {
    # only do first if !runNoOA
    if (!runNoOA && (i>1)) break
    oaVal   = oaValV[i]
    prefix1 = combine(prefix, combine(paramStr, oaVal))
    # set oaEffectV
    oaEffectV = rep(1,length(oaStrV))
    oaEffectV[oaIndex] = oaVal
    if (multiIndex) {
      # Add the second effect
      oaEffectV[index2] = oaValV2[i]	  
	}
	
	
	# loop until the param instance survives
	while (TRUE) {
	  myP("Running test replicate #", repNum)
	  # some issues with single row matrices I'd rather not bother fixing
	  # just duplicating instead
	  singleDispM = cbind(dispSDM[,repNum], dispSDM[,repNum])
	  singleRecM  = cbind(randRecM[,repNum], randRecM[,repNum])

	  pL = runDensityTest(hV=hV, prefix=prefix1, numReplic=2, maxT=maxT, manV=manV, oaEffectV=oaEffectV, addCat=addCat,
	                      catMort=0.5, dispSDM=singleDispM, randRecM=singleRecM, showPlots=FALSE, showMetrics=FALSE)
	  if ((i>1)||(repNum==100)) break
	  if (pL$propExtinct1==0) break
	  else repNum = repNum + 1
	}
	if (is.null(metricD)) metricD = as.data.frame(pL)
	else metricD = rbind(metricD, as.data.frame(pL))
  }
  
  paramStr = combine(prefix, paramStr)
  saveMatrixToCSVFile(combine(combine("SingleOA", paramStr), "metrics"), metricD)
}

showSizeDiff <- function(prefix1="MA.Surv._0.98", prefix2="NoOA", classWidth=5, dropYOY=TRUE, propOnly=TRUE) {
  fileStr = combine(combine("Single", prefix1), "NC_size")
  # load and drop first column
  compM   = loadMatrixFromCSVFile(fileStr)
  fileStr = combine(combine("Single", prefix2), "NC_size")
  baseM   = loadMatrixFromCSVFile(fileStr)
  sizeV   = gMeanSizeV  
  
  if (dropYOY) {
    baseM = baseM[-1:-2,-1]
    compM = compM[-1:-2,-1]
	sizeV = sizeV[-1:-2]
  }
  
  maxT = ncol(baseM)
  if (propOnly) {
    for (i in 1:maxT) {
      baseM[,i] = baseM[,i] / sum(baseM[,i])
      compM[,i] = compM[,i] / sum(compM[,i])
	}
  }
  
  # All size classes diff over time
  # Negative indicates lower proportion as compared to base
  allDiffM = as.matrix(compM - baseM)
  
  numClasses = length(sizeV) %/% classWidth
  if (length(sizeV) %% classWidth > 0) numClasses = numClasses + 1
  
  # every X classes (of 43)
  classXM = matrix(nrow=numClasses, ncol=maxT)
  for (i in 1:numClasses) {
    start = (i-1)*classWidth + 1
	end   = start + classWidth - 1
    if (i==numClasses) end = length(sizeV)
	
	classXM[i,] = apply(allDiffM[start:end,], 2, sum)
  }
  
  
  # 3 major classes: J <50mm, YA <136mm, MA >136
  class3M    = matrix(nrow=3, ncol=maxT)
  class3M[1,] = apply(allDiffM[1:7,], 2, sum)
  class3M[2,] = apply(allDiffM[8:25,], 2, sum)
  class3M[3,] = apply(allDiffM[26:length(sizeV),], 2, sum)
  
  # convert to tidy data
  rownames(allDiffM) <- sizeV
  colnames(allDiffM) <- 1:maxT
  rownames(classXM) <- 1:numClasses
  colnames(classXM) <- 1:maxT
  rownames(class3M) <- c("Juv.", "YA", "MA")
  colnames(class3M) <- 1:maxT
  allDiffD = melt(allDiffM, varnames=c("size", "time"), value.name="diff")
  classXD  = melt(classXM, varnames=c("size", "time"), value.name="diff")
  class3D  = melt(class3M, varnames=c("size", "time"), value.name="diff")
  
  thePlot = ggplot(data=allDiffD, aes(x=time, y=diff)) + geom_line(aes(color=size)) +
            labs(title=paste(prefix1, "All classes", sep=":"), x="Time", y="Diff. from base")
  myPPlot(thePlot)

  thePlot = ggplot(data=classXD, aes(x=time, y=diff)) + geom_line(aes(color=size)) +
            labs(title=paste(prefix1, paste("Class width=", classWidth), sep=":"), x="Time", y="Diff. from base")
  myPPlot(thePlot)

  thePlot = ggplot(data=class3D, aes(x=time, y=diff)) + geom_line(aes(color=size)) +
            labs(title=paste(prefix1, "Juv, YA, MA", sep=":"), x="Time", y="Diff. from base")
  myPPlot(thePlot)

}
# Saves a GIF showing full size structure over time, comparing two instances

saveSizeMovie <- function(prefix1="MA.Surv._0.98", prefix2="NoOA", skipV=c(), delFrames=TRUE, useLines=TRUE, inColor=FALSE, propOnly=TRUE, dropYOY=TRUE) {
  myP("Making movie...")
  fileStr = combine(combine("Single", prefix1), "NC_size")
  compD   = loadMatrixFromCSVFile(fileStr)
  fileStr = combine(combine("Single", prefix2), "NC_size")
  baseD   = loadMatrixFromCSVFile(fileStr)
  maxT    = ncol(baseD)
    
  if (dropYOY) {
    baseD = baseD[-1:-2,]
    compD = compD[-1:-2,]
  }
  
  if (propOnly) {
    for (i in 1:maxT) {
      baseD[,i] = baseD[,i] / sum(baseD[,i])
      compD[,i] = compD[,i] / sum(compD[,i])
	}
  }
  maxVal = max(baseD)
  if (max(compD)>maxVal) maxVal = max(compD)
  
  # function to draw each frame and save as PNG
  drawTimeStep <- function(ts) {
    if (ts %in% skipV) {
	  myP("Skipping frame t=", ts, "...")
	  return()
	}
    myP("Drawing frame t=", ts, "...")
	tsStr = as.character(ts)
	if (ts<10) tsStr = paste("0", tsStr, sep="")
	if (ts<100) tsStr = paste("0", tsStr, sep="")
	filename = paste(combine(prefix1, combine("timestep", tsStr)), ".png", sep="")
	png(file=filename, width=1200, height=800)
	
	tsD    = data.frame(size=gMeanSizeV[-1:-2], base=baseD[,ts], comp=compD[,ts])
	tsPlot = ggplot(data=tsD, aes(x=size)) + theme_bw()
	# choose between lines or stacked histogram
	if (useLines) {
      tsPlot = tsPlot + geom_line(aes(y=base), color="black", size=3)
      tsPlot = tsPlot + geom_line(aes(y=comp), color="gray", size=3)
	} else {
	  tsPlot = tsPlot + geom_bar(stat="identity", aes(y=base), size=2) 
	  tsPlot = tsPlot + geom_bar(stat="identity", aes(y=comp), size=2) 
	}
	tsPlot = tsPlot + scale_y_continuous(limits=c(0, maxVal))
	tsPlot = tsPlot + labs(title=combine(prefix1, paste("t=",ts,sep="")), x="Size", y="Proportion") +
	                  theme(plot.title=element_text(size=rel(5)),
		  				    legend.text=element_text(size=rel(4)), 
							legend.title=element_text(size=rel(4)),
							axis.text.x=element_text(size=rel(6), color="black"), 
	                        axis.text.y=element_text(size=rel(6), color="black"),
	   				        axis.title.x=element_text(size=rel(5)),
							axis.title.y=element_text(size=rel(5)))
    # draw cell plot
	print(tsPlot)
	dev.off()
  }

  # loop function for use in saveGIF()
  loopTimeSteps <- function() {
	lapply(1:maxT, function(t) { drawTimeStep(t) })
  }
    
  # saveGIF PATH info not working quite right, do by hand
  loopTimeSteps()
  # hardcode PATH info
  # NOTE: incorporates ALL PNG files
  #       Must delete earlier files before running again
  system(paste('"C:\\Program Files\\ImageMagick\\convert.exe" -delay 25 *.png', 
	           paste("D:\\Research\\Code\\abalone\\", combine(prefix1, "spread.gif"), sep="")))
  if (delFrames) file.remove(list.files(pattern=".png"))
#    ani.options(convert="C:/Program Files/ImageMagick-6.9.2-Q16/convert.exe");
#	saveGIF(loopTimeSteps(), interval=0.35, movie.name=combine(prefix, "spread.gif"));
	#movie.name=paste("C:\\Research\\Code\\disease\\", combine(prefix, "spread.gif"), sep=""))
}

runOAContinuum <- function(numReplic=20, maxT=300, calcOnly=FALSE, oaIndex=OA_GROWTH_RATE_INDEX, isTest=FALSE, prefix="",
                           run5YCOnly=FALSE, addCat=FALSE, randSeed=RAND_SEED, oaValV=NULL) {
  set.seed(randSeed)
  dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma")
  randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=REC_LN_SDLOG)

  manV = c(MAN_NC)#, MAN_MPA)
  if (run5YCOnly) manV = c(MAN_5YC)

  metricD  = NULL
  paramStr = oaStrV[oaIndex]
  myP("Running OA models for", paramStr)

  hV = c(getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_NC, propMSY=MSY_PROP), 
	     getMSYValue(paramStr="Base", paramVal=1.0, manType=MAN_MPA, propMSY=1.0))
	
  multiIndex = FALSE
  if (is.null(oaValV)) {
    if (oaIndex==OA_SIZE_MAT_INDEX) oaValV = c(1, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2) 
#  if (oaIndex==OA_SIZE_MAT_INDEX) oaValV = c(1, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8)  
    else if ((oaIndex==OA_GROWTH_RATE_INDEX)||(oaIndex==OA_MA_SURV_INDEX)||(oaIndex==OA_MAX_SIZE_INDEX)) oaValV = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
    else if ((oaIndex>=OA_MULTI_FERT_JSURV)) {  # XXX special case for multi effects
      multiIndex = TRUE
	  if (oaIndex==OA_MULTI_FERT_JSURV) {
	    oaValV1 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    oaValV2 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_JUV_SURV_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_SMAT) {
	    oaValV1 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    oaValV2 = c(1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.075, 1.1)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_FERT_GROW) {
	    oaValV1 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    oaValV2 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    index1  = OA_FERT_SUCCESS_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_SMAT) {
	    oaValV1 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    oaValV2 = c(1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.075, 1.1)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_SIZE_MAT_INDEX
	  } else if (oaIndex==OA_MULTI_JSURV_GROW) {
	    oaValV1 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    oaValV2 = c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9)
	    index1  = OA_JUV_SURV_INDEX
	    index2  = OA_GROWTH_RATE_INDEX
	  }
	  oaValV  = oaValV1
	  oaIndex = index1
    } else oaValV = c(1, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8)
  }

  
  # loop through parameter values
  for (i in 1:length(oaValV)) {
    oaVal   = oaValV[i]
    prefix1 = combine(prefix, combine(paramStr, oaVal))
    # set oaEffectV
    oaEffectV = rep(1,length(oaStrV))
    oaEffectV[oaIndex] = oaVal
    if (multiIndex) {
      # Add the second effect
      oaEffectV[index2] = oaValV2[i]	  
	}
	pL = runDensityTest(hV=hV, prefix=prefix1, numReplic=numReplic, maxT=maxT, manV=manV, oaEffectV=oaEffectV, addCat=addCat,
	                    catMort=0.5, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE, showMetrics=FALSE)
	if (is.null(metricD)) metricD = as.data.frame(pL)
	else metricD = rbind(metricD, as.data.frame(pL))
  }
  
  if (isTest) paramStr = combine("Test", paramStr)
  else paramStr = combine(prefix, paramStr)
  saveMatrixToCSVFile(combine(combine("OACont", paramStr), "metrics"), metricD)
  # save the replicate comparison files
  for (manType in manV) compareLikeReplicates(oaIndex, oaValV, prefix, manType)
}

compareLikeReplicates <- function(oaIndex=OA_GROWTH_RATE_INDEX, rangeV=c(1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.925, 0.9), prefix="", manType=MAN_NC) {
  paramStr = oaStrV[oaIndex]
  prefix1  = combine(combine(prefix, combine(paramStr, rangeV[1])), manType)
  fileStr  = combine(prefix1, "repData")
  baseRepD = loadMatrixFromCSVFile(fileStr)
  
  # for each oa value, calculate the following:
  # oa val, num surv, mean surv abun, abunSSB, catch, size, sizeSSB, mean fec, total fec
  # AND each of those relative to *only* the equivalent replicates in the base set
  
  repCompareM = matrix(nrow=length(rangeV), ncol=16)
  survV    = !baseRepD$collapse
  survD    = baseRepD[survV,]
  # create mean values
  baseRowV = c(rangeV[1], sum(survV), 
               mean(survD$abun), mean(survD$abunSSB), mean(survD$catch), 
			   mean(survD$size), mean(survD$sizeSSB), mean(survD$meanFec), mean(survD$totalFec))
  # add relative values (all 1s here)
  baseRowV = c(baseRowV, c(1,1,1,1,1,1,1))	
  repCompareM[1,] = baseRowV

  for (i in 2:length(rangeV)) {
    # get new oa val
    prefix2  = combine(combine(prefix, combine(paramStr, rangeV[i])), manType)
    fileStr  = combine(prefix2, "repData")
    nextRepD = loadMatrixFromCSVFile(fileStr)

	survV    = !nextRepD$collapse
	if (sum(survV)<0) {
	  # no surviving values
	  nextRowV = c(rangeV[i], rep(0,15))
	  repCompareM[i,] = nextRowV
	  next
	}
	
    survD    = nextRepD[survV,]
    # create mean values
    nextRowV = c(rangeV[i], sum(survV), 
                 mean(survD$abun), mean(survD$abunSSB), mean(survD$catch), 
		  	     mean(survD$size), mean(survD$sizeSSB), mean(survD$meanFec), mean(survD$totalFec))
	# IMPORTANT: compares only to same replicates from baseline
	baseSurvD   = baseRepD[survV,]
	baseSubsetV = c(mean(baseSurvD$abun), mean(baseSurvD$abunSSB), mean(baseSurvD$catch), 
		  	        mean(baseSurvD$size), mean(baseSurvD$sizeSSB), mean(baseSurvD$meanFec), mean(baseSurvD$totalFec))
    # add relative values 
    nextRowV  = c(nextRowV, nextRowV[3:9]/baseSubsetV)	
    repCompareM[i,] = nextRowV
  }
  
  colnames(repCompareM) <- c("paramVal", "#Surv", "abun", "abunSSB", "catch", "size", "sizeSSB", "meanFec", "totalFec",
                             "rel_abun", "rel_abunSSB", "rel_catch", "rel_size", "rel_sizeSSB", "rel_meanFec", "rel_totalFec")
  fileStr = combine(prefix, combine(combine(paramStr, manType), "repCompare"))
  saveMatrixToCSVFile(fileStr, repCompareM)
}

runAllTests <- function(numReplic=100, maxT=300) {
  dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma");
  randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=REC_LN_SDLOG);

  runBaseModel(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM)
#  runCatastropheTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);
#  runCarryingCapacityTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);
#  runMultCatDelayTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);
#  runDispersalTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);
#  runAlleeTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);
#  runSizeTests(numReplic, maxT, dispSDM=dispSDM, randRecM=randRecM);

#  runClassWidthTests(numReplic, maxT);
}

runClassWidthTests <- function(numReplic=5, maxT=300) {
  dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma");
  randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=REC_LN_SDLOG);
  randRecM = NULL;
  lowMort  = 0.5;
  medMort  = 0.75;
  highMort = 0.90;
  
  lowK     = 2.59*10^7; # 0.2/m2
  medK     = 5.43*10^7; # 0.4/m2
  highK    = 1.10*10^8; # 0.8/m2
  
  lowKhV   = c(0, 0.09, 0.1, 0.12);
  medKhV   = c(0, 0.15, 0.16, 0.22);
  highKhV  = c(0, 0.18, 0.2, 0.25);
  
  # lowK
  myP("Running cW=5 models...");
  pL1 = runDensityTest(cW=5, hV=medKhV, blockK=medK, catMort=lowMort, prefix="lowK_lowM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL2 = runDensityTest(cW=5, hV=medKhV, blockK=medK, catMort=medMort, prefix="lowK_medM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL3 = runDensityTest(cW=5, hV=medKhV, blockK=medK, catMort=highMort, prefix="lowK_highM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);

  # medK
  myP("Running cW=1 models (same h)...");
  pL4 = runDensityTest(cW=1, hV=medKhV, blockK=medK, catMort=lowMort, prefix="medK_lowM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL5 = runDensityTest(cW=1, hV=medKhV, blockK=medK, catMort=medMort, prefix="medK_medM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL6 = runDensityTest(cW=1, hV=medKhV, blockK=medK, catMort=highMort, prefix="medK_highM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);

  # highK
  myP("Running cW=1 models (1/2 h)...");
  pL7 = runDensityTest(cW=1, hV=medKhV/2, blockK=medK, catMort=lowMort, prefix="highK_lowM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL8 = runDensityTest(cW=1, hV=medKhV/2, blockK=medK, catMort=medMort, prefix="highK_medM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  pL9 = runDensityTest(cW=1, hV=medKhV/2, blockK=medK, catMort=highMort, prefix="highK_highM", numReplic=numReplic, maxT=maxT, dispSDM=dispSDM, randRecM=randRecM, showPlots=FALSE);
  
  # create 3x3 plot for density
  plotL = list(p1=pL1$densP, p2=pL2$densP, p3=pL3$densP,
			   p4=pL4$densP, p5=pL5$densP, p6=pL6$densP,
			   p7=pL7$densP, p8=pL8$densP, p9=pL9$densP);
  plot3x3Grid(plotL, "CW_dense");

  # create 3x3 plot for catch
  plotL = list(p1=pL1$cumCatchP, p2=pL2$cumCatchP, p3=pL3$cumCatchP,
			   p4=pL4$cumCatchP, p5=pL5$cumCatchP, p6=pL6$cumCatchP,
			   p7=pL7$cumCatchP, p8=pL8$cumCatchP, p9=pL9$cumCatchP);
  plot3x3Grid(plotL, "CW_cumCatch");

}

tuneUnfishedDensity <- function(targD=0.4, numReplic=3, maxT=50) {
  dispSDM = getRandomDispSDMatrix(maxT, numReplic, rType="gamma")
  # loop to find K value
  theK = DEF_K
  pL   = getParamList(hProp=0, mpaProp=0, blockK=theK, maxT=maxT, numReplic=numReplic, addCat=FALSE)
  totalRuns = 0
  aveFrom   = maxT-20
  adjBy     = 2
  wasLarger = FALSE
  while (TRUE) {
    myP("Running with K=", theK)
	pL$blockK = theK
    noHarvL   = runModelFromList(pL, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM)
	meanD     = noHarvL$means
	meanDens  = mean(meanD[meanD$t>aveFrom,]$meanSSB) / (TOTAL_AREA*10000)
	myP("Density:", meanDens)
	if (my.eq(targD, meanDens, minDiff=targD*0.01)) break
	else if (totalRuns>50) {
	  myP("Max runs reached.")
	  break
	} else if (meanDens>targD) {
	  if (!wasLarger) {
	    adjBy = adjBy * 0.5
		wasLarger = TRUE
	  }
	  theK = theK / (1+adjBy)
	} else {
	  if (wasLarger) {
	    adjBy = adjBy * 0.5
		wasLarger = FALSE
	  }
	  theK = theK * (1+adjBy)
	}
  }

  myP("Final K:", theK, ", density=", meanDens)
}

runDensityTest <- function(numReplic=5, maxT=300, hV=c(0.06,0.07), mpaSize=MPA_SIZE, blockK=DEF_K, 
                           aggA=AGG_A, randRec=FALSE, mpaProp=MPA_PROP, manV=c(MAN_NC, MAN_DC, MAN_MPA), oaEffectV=rep(1, length(oaStrV)),
						   addCat=TRUE, catMort=DEF_CAT_MORT, numCat=1, catInt=10, harvDelay=DYNAMIC_DELAY, closeThresh=DYNAMIC_THRESH,
						   gShape=DISP_GAMMA_SHAPE, rShape=REC_LN_SDLOG, prefix="", dispSDM=NULL, randRecM=NULL, discRate=DISCOUNT_RATE,
						   showPlots=TRUE, sepExtinct=TRUE, plotBoth=FALSE, showMetrics=TRUE,
						   showLegend=FALSE, zeroTime=TRUE, numbersOnly=TRUE, saveMatrix=TRUE) {
  if (is.null(dispSDM)) dispSDM = getRandomDispSDMatrix(maxT, numReplic, rType="gamma", shape=gShape)
  # if matrix is provided, use random recruitment
  if (!is.null(randRecM)) randRec = TRUE
  else randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=rShape)

  # set management parameters
  runNC  = MAN_NC %in% manV 
  runDC  = MAN_DC %in% manV 
  run5YC = MAN_5YC %in% manV
  runMPA = MAN_MPA %in% manV 

  # XXX assume here that a negative aggA value means no Allee effect
  isAllee = aggA>0
  
  if (runNC) pL1 = getParamList(hProp=hV[1], hPresim=hV[1], mpaProp=0, noHarvPostCat=0, blockK=blockK, aggA=aggA, isAllee=isAllee, oaEffectV=oaEffectV, randRec=randRec, maxT=maxT, numReplic=numReplic, addCat=addCat, catMort=catMort, numCat=numCat, catInt=catInt)
  if (runDC) pL2 = getParamList(hProp=hV[1], hPresim=hV[1], mpaProp=0, noHarvPostCat=harvDelay, closeThresh=closeThresh, blockK=blockK, aggA=aggA, isAllee=isAllee, oaEffectV=oaEffectV, randRec=randRec, maxT=maxT, numReplic=numReplic, addCat=addCat, catMort=catMort, numCat=numCat, catInt=catInt)
  else if (run5YC) pL2 = getParamList(hProp=hV[1], hPresim=hV[1], mpaProp=0, noHarvPostCat=5, closeThresh=closeThresh, blockK=blockK, aggA=aggA, isAllee=isAllee, oaEffectV=oaEffectV, randRec=randRec, maxT=maxT, numReplic=numReplic, addCat=addCat, catMort=catMort, numCat=numCat, catInt=catInt)
  # Use NC harvest for presim 
  if (runMPA) pL3 = getParamList(hProp=hV[2], hPresim=hV[1], mpaProp=mpaProp, mpaSize=mpaSize, blockK=blockK, aggA=aggA, isAllee=isAllee, randRec=randRec, maxT=maxT, numReplic=numReplic, addCat=addCat, catMort=catMort, numCat=numCat, catInt=catInt)
  
  # run models
  if (runNC) harvL1 = runModelFromList(pL1, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine(prefix, MAN_NC))
  if (runDC) harvL2 = runModelFromList(pL2, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine(prefix, MAN_DC))
  else if (run5YC) harvL2 = runModelFromList(pL2, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine(prefix, MAN_5YC))
  if (runMPA) harvL3 = runModelFromList(pL3, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine(prefix, MAN_MPA))

  # get results
  if (runNC) harvD1 = harvL1$means
  if (runDC) harvD2 = harvL2$means
  else if (run5YC) harvD2 = harvL2$means
  if (runMPA) harvD3 = harvL3$means

  # convert to mean densities  
  # adds labels, separates by reserve and fished, creates overall, survived and extinct time series
  if (runNC) densResL1 = calcSpatialDensity(harvD1, harvL1$res, pL1)
  if (runDC) densResL2 = calcSpatialDensity(harvD2, harvL2$res, pL2)
  else if (run5YC) densResL2 = calcSpatialDensity(harvD2, harvL2$res, pL2)
  if (runMPA) densResL3 = calcSpatialDensity(harvD3, harvL3$res, pL3)
  
  # Dynamic closure data
  if (runDC) closureV  = harvL2$res$closureLengthV
  
  # combine into one data frame
  allMeanD = data.frame()
  if (runNC) allMeanD = rbind(allMeanD, densResL1$densD)
  if (runDC) allMeanD = rbind(allMeanD, densResL2$densD)
  else if (run5YC) allMeanD = rbind(allMeanD, densResL2$densD)
  if (runMPA) allMeanD = rbind(allMeanD, densResL3$densD)
  allMeanD$ssbDensity = allMeanD$meanSSB / (TOTAL_AREA*10000)
  allMeanD$catchHa    = allMeanD$meanCatch / (TOTAL_AREA)
 
  if (zeroTime) {
    # set t=0 to the first catastrophe
    allMeanD$t = allMeanD$t - (maxT/2)
	maxT = maxT/2
    startTime = -5
  } else startTime = maxT/2 - 5
  
  metricL = calcDensityMetrics(allMeanD, prefix, manV=manV, discRate=discRate)
  if (runDC) {
    metricL$meanClosure = mean(closureV)
    metricL$sdClosure   = sd(closureV)
  } else if (run5YC) {
    metricL$meanClosure = 5
    metricL$sdClosure   = 0
  }
  
  # store the time series
  if (saveMatrix) saveMatrixToCSVFile(combine(prefix, "ts"), allMeanD)
  if (showPlots) plotL = plotDensityPlots(allMeanD, startTime, manV, prefix,
                                     sepExtinct=sepExtinct, plotBoth=plotBoth, numbersOnly=numbersOnly, showPlots=TRUE)
  if (showMetrics) showDensityMetrics(metricL, prefix)
  return(metricL)
}

plotDensityPlots <- function(allMeanD, startTime=-5, manV=c(MAN_NC, MAN_DC, MAN_MPA), prefix="",
                        sepExtinct=TRUE, plotBoth=FALSE, numbersOnly=FALSE, showPlots=TRUE) {
  DEF_W = 10;
  DEF_H = 7.25;
  
  pName1 = manV[1]
  pName2 = manV[2]
  pName3 = manV[3]
  
  metricL = calcDensityMetrics(allMeanD, prefix)

  maxT = max(allMeanD$t)
  # draw SSB density plot
  maxY = max(allMeanD$resDensSurv, na.rm=TRUE)*1.1;
  if (maxY<0) maxY = max(allMeanD$resDensExt, na.rm=TRUE)*1.1;
  posX = maxT*0.9;
  posY = maxY*0.15;
  if (sepExtinct || plotBoth) {
    densPlot = ggplot(allMeanD[allMeanD$t>=startTime,], aes(x=t)) + theme_bw() +
               geom_line(aes(y=resDensSurv, group=name, color=name, linetype="Reserve"), size=2) +
	   	       geom_line(aes(y=fishDensSurv, group=name, color=name, linetype="Fishery"), size=2) +
               geom_line(aes(y=resDensExt, group=name, color=name, linetype="Extinct"), size=2) +
			   geom_line(aes(y=fishDensExt, group=name, color=name, linetype="Extinct"), size=2) +
			   scale_y_continuous(limits=c(0, maxY)) +
			   scale_x_continuous() +
			   scale_color_manual(values=c("red", "green", "blue"), breaks=c(pName1, pName2, pName3)) +
			   scale_linetype_manual(values=c("Fishery"="solid", "Reserve"="dashed", "Extinct"="dotted"), 
			                         breaks=c("Fishery", "Reserve", "Extinct")) +
  		       theme(legend.text=element_text(size=rel(1.25)), 
	                legend.title=element_text(size=0),
        	        legend.position=c(.90, .2)) +
               theme(axis.text.x=element_text(size=rel(1.25)), 
                     axis.text.y=element_text(size=rel(1.25)),
 	                 axis.title.x=element_text(size=rel(1.25)),
 			         axis.title.y=element_text(size=rel(1.25))) +
			   annotate("text", label=paste(pName1, format(metricL$propExtinct1, digits=2)), x=posX*0.8, y=posY, size=rel(4), color="black") +
			   annotate("text", label=paste(pName2, format(metricL$propExtinct2, digits=2)), x=posX*0.8, y=posY*0.35, size=rel(4), color="black") +
			   annotate("text", label=paste(pName3, format(metricL$propExtinct3, digits=2)), x=posX*0.8, y=posY*0.7, size=rel(4), color="black");
	if (!numbersOnly) densPlot = densPlot + annotate("text", label="Prop. extinct:", x=posX*0.8, y=posY*1.25, size=rel(4), color="black");
	
    densPlot = densPlot + labs(x="Year", y="Density (#/m^2)", title=paste(prefix, "Spawner density"));
	
    if (showPlots) myPPlot(densPlot, w=DEF_W, h=DEF_H);
    fileName = paste("output/", combine(prefix, "ssbDensExt.png"), sep="");
    myPTiff(densPlot, fileName, h=4, w=6);    
  }
  
  if (!sepExtinct || plotBoth) {
    densPlot = ggplot(allMeanD[allMeanD$t>=startTime,], aes(x=t)) + theme_bw() +
               geom_line(aes(y=resDens, group=name, color=name, linetype="Reserve"), size=2) +
		       geom_line(aes(y=fishDens, group=name, color=name, linetype="Fishery"), size=2) +
			   scale_linetype_manual(values=c("Fishery"="solid", "Reserve"="dashed"), breaks=c("Fishery", "Reserve")) +
			   scale_y_continuous(limits=c(0, maxY)) +
  		       theme(legend.text=element_text(size=rel(1.25)), 
	                legend.title=element_text(size=0),
        	        legend.position=c(.90, .2)) +
               theme(axis.text.x=element_text(size=rel(1.25)), 
                     axis.text.y=element_text(size=rel(1.25)),
 	                 axis.title.x=element_text(size=rel(1.25)),
 			         axis.title.y=element_text(size=rel(1.25)));
    densPlot = densPlot + labs(x="Year", y="Density (#/m^2)", title=paste(prefix, "Spawner density"));
 #   if (!showLegend) densPlot = stripPlot(densPlot, stripLegend=TRUE);
    if (showPlots) myPPlot(densPlot, w=DEF_W, h=DEF_H);
    fileName = paste("output/", combine(prefix, "ssbDens.png"), sep="");
    myPTiff(densPlot, fileName, h=4, w=6);    
  }

  # mean catch overall only
  # draw catch kg plot
  if (sepExtinct || plotBoth) {
    maxY = max(allMeanD$catchSurv, na.rm=TRUE)*1.1;
    if (maxY<0) maxY = max(allMeanD$catchExt, na.rm=TRUE)*1.1;
    posY = maxY*0.15;
    catchPlot = ggplot(allMeanD, aes(x=t)) + theme_bw() +
               geom_line(aes(y=catchSurv, group=name, color=name, linetype="Fishery"), size=2) +
			   geom_line(aes(y=catchExt, group=name, color=name, linetype="Extinct"), size=2) +
			   scale_y_continuous(limits=c(0, maxY)) +
			   scale_color_manual(values=c("red", "green", "blue"), breaks=c(pName1, pName2, pName3)) +
			   scale_linetype_manual(values=c("Fishery"="solid", "Extinct"="dotted"), breaks=c("Fishery", "Extinct")) +
  		       theme(legend.text=element_text(size=rel(1.25)), 
	                legend.title=element_text(size=0),
        	        legend.position=c(.90, .2)) +
               theme(axis.text.x=element_text(size=rel(1.25)), 
                     axis.text.y=element_text(size=rel(1.25)),
 	                 axis.title.x=element_text(size=rel(1.25)),
 			         axis.title.y=element_text(size=rel(1.25))) +
			   annotate("text", label=paste(pName1, format(metricL$propExtinct1, digits=2)), x=posX*0.8, y=posY, size=rel(4), color="black") +
			   annotate("text", label=paste(pName2, format(metricL$propExtinct2, digits=2)), x=posX*0.8, y=posY*0.35, size=rel(4), color="black") +
			   annotate("text", label=paste(pName3, format(metricL$propExtinct3, digits=2)), x=posX*0.8, y=posY*0.7, size=rel(4), color="black") +
               labs(title=paste(prefix, "Catch per hectare"), x="Year", y="Catch in kg.");
	if (!numbersOnly) catchPlot = catchPlot + annotate("text", label="Prop. extinct:", x=posX*0.8, y=posY*1.25, size=rel(4), color="black");
#    if (!showLegend) catchPlot = stripPlot(catchPlot, stripLegend=TRUE);
    fileName = paste("output/", combine(prefix, "catchHaExt.png"), sep="");
    myPTiff(catchPlot, fileName, h=4, w=6);    
  }

  if (!sepExtinct || plotBoth) {
    catchPlot = ggplot(allMeanD, aes(x=t)) + theme_bw() +
               geom_line(aes(y=catchHa, group=name, color=name, linetype=name), size=2) +
  		       theme(legend.text=element_text(size=rel(1.25)), 
	                legend.title=element_text(size=0),
        	        legend.position=c(.90, .2)) +
               theme(axis.text.x=element_text(size=rel(1.25)), 
                     axis.text.y=element_text(size=rel(1.25)),
 	                 axis.title.x=element_text(size=rel(1.25)),
 			         axis.title.y=element_text(size=rel(1.25)));
    catchPlot = catchPlot + labs(title=paste(prefix, "Catch per hectare"), x="Year", y="Catch in kg.");
#    if (showPlots) myPPlot(catchPlot, w=DEF_W, h=DEF_H);
    fileName = paste("output/", combine(prefix, "catchHa.png"), sep="");
    myPTiff(catchPlot, fileName, h=4, w=6);    
  }

  # cumulative catch overall only
  # draw cumulative catch kg from catastrophe
  # Just look at 60 years after catastrophe
  postCatD = allMeanD[allMeanD$t>(startTime-5),];
  postCatD = postCatD[postCatD$t<(startTime+65),];
  if (sepExtinct || plotBoth) {
    posX = startTime + 50;
    maxY = max(postCatD$cumCatchSurv, na.rm=TRUE)*1.1;
    if (maxY<0) maxY = max(postCatD$cumCatchExt, na.rm=TRUE)*1.1;
    posY = maxY*0.15;
    cumCatchPlot = ggplot(postCatD, aes(x=t)) + theme_bw() +
               geom_line(aes(y=cumCatchSurv, group=name, color=name, linetype="Fishery"), size=2) +
			   geom_line(aes(y=cumCatchExt, group=name, color=name, linetype="Extinct"), size=2) +
			   scale_y_continuous(limits=c(0, maxY)) +
			   scale_color_manual(values=c("red", "green", "blue"), breaks=c(pName1, pName2, pName3)) +
			   scale_linetype_manual(values=c("Fishery"="solid", "Extinct"="dotted"), breaks=c("Fishery", "Extinct")) +
  		       theme(legend.text=element_text(size=rel(1.25)), 
	                legend.title=element_text(size=0),
        	        legend.position=c(.90, .2)) +
               theme(axis.text.x=element_text(size=rel(1.25)), 
                     axis.text.y=element_text(size=rel(1.25)),
 	                 axis.title.x=element_text(size=rel(1.25)),
 			         axis.title.y=element_text(size=rel(1.25))) +
			   annotate("text", label=paste(pName1, format(metricL$propExtinct1, digits=2)), x=posX*0.8, y=posY, size=rel(4), color="black") +
			   annotate("text", label=paste(pName2, format(metricL$propExtinct2, digits=2)), x=posX*0.8, y=posY*0.35, size=rel(4), color="black") +
			   annotate("text", label=paste(pName3, format(metricL$propExtinct3, digits=2)), x=posX*0.8, y=posY*0.7, size=rel(4), color="black") +
               labs(title=paste(prefix, "Cumulative catch"), x="Year", y="Catch in kg.");
	if (!numbersOnly) cumCatchPlot = cumCatchPlot + annotate("text", label="Prop. extinct:", x=posX*0.8, y=posY*1.25, size=rel(4), color="black");
#    if (!showLegend) cumCatchPlot = stripPlot(cumCatchPlot, stripLegend=TRUE);
    if (showPlots) myPPlot(cumCatchPlot, w=DEF_W, h=DEF_H);
    fileName = paste("output/", combine(prefix, "cumCatchExt.png"), sep="");
    myPTiff(cumCatchPlot, fileName, h=4, w=6);    
  }

  if (!sepExtinct || plotBoth) {
    cumCatchPlot = ggplot(postCatD, aes(x=t)) + theme_bw() +
                   geom_line(aes(y=cumCatch, group=name, color=name), size=2) +
     		       theme(legend.text=element_text(size=rel(1.25)), 
	                     legend.title=element_text(size=0),
                         legend.position=c(.90, .2)) +
                   theme(axis.text.x=element_text(size=rel(1.25)), 
                         axis.text.y=element_text(size=rel(1.25)),
 	                     axis.title.x=element_text(size=rel(1.25)),
 			             axis.title.y=element_text(size=rel(1.25)));
    cumCatchPlot = cumCatchPlot + labs(title=paste(prefix, "Cumulative catch"), x="Year", y="Catch in kg.");
#    if (!showLegend) cumCatchPlot = stripPlot(cumCatchPlot, stripLegend=TRUE);
    if (showPlots) myPPlot(cumCatchPlot, w=DEF_W, h=DEF_H);
    fileName = paste("output/", combine(prefix, "cumCatchExt.png"), sep="");
    myPTiff(cumCatchPlot, fileName, h=4, w=6);    
  }

  return(c(list(densP=densPlot, catchP=catchPlot, cumCatchP=cumCatchPlot), metricL))
}

calcDensityMetrics <- function(allMeanD, prefix, manV=c(MAN_NC, MAN_DC, MAN_MPA), aveOver=10, recPoint=20,
                               sepExtinct=TRUE, discRate=DISCOUNT_RATE) {
  runNC  = MAN_NC %in% manV 
  runDC  = MAN_DC %in% manV 
  run5YC = MAN_5YC %in% manV
  runMPA = MAN_MPA %in% manV 

  if (runNC) pName1 = combine(prefix, MAN_NC)
  if (runDC) pName2 = combine(prefix, MAN_DC)
  else if (run5YC) pName2 = combine(prefix, MAN_5YC)
  if (runMPA) pName3 = combine(prefix, MAN_MPA)

  maxT = max(allMeanD$t)
 
  # prop extinct value is time-independent
  if (runNC) propExtinct1 = allMeanD[allMeanD$name==pName1,]$propExtinct[1]
  if (runDC) propExtinct2 = allMeanD[allMeanD$name==pName2,]$propExtinct[1]
  else if (run5YC) propExtinct2 = allMeanD[allMeanD$name==pName2,]$propExtinct[1]
  if (runMPA) propExtinct3 = allMeanD[allMeanD$name==pName3,]$propExtinct[1]
  
  # calculate differences in equilibrium SSB, catch
  if (runNC) p1D = allMeanD[allMeanD$name==pName1,]
  if (runDC) p2D = allMeanD[allMeanD$name==pName2,]
  else if (run5YC) p2D = allMeanD[allMeanD$name==pName2,]
  if (runMPA) p3D = allMeanD[allMeanD$name==pName3,]

  # only use ssb density and catch in surviving runs
  if (runNC) meanSSB1   = mean(p1D[p1D$t>(maxT-aveOver),]$fishDensSurv)
  if (runDC) meanSSB2   = mean(p2D[p2D$t>(maxT-aveOver),]$fishDensSurv)
  else if (run5YC) meanSSB2   = mean(p2D[p2D$t>(maxT-aveOver),]$fishDensSurv)
  if (runMPA) meanSSB3   = mean(p3D[p3D$t>(maxT-aveOver),]$fishDensSurv)
  if (runNC) meanCatch1 = mean(p1D[p1D$t>(maxT-aveOver),]$catchSurv)
  if (runDC) meanCatch2 = mean(p2D[p2D$t>(maxT-aveOver),]$catchSurv)
  else if (run5YC) meanCatch2 = mean(p2D[p2D$t>(maxT-aveOver),]$catchSurv)
  if (runMPA) meanCatch3 = mean(p3D[p3D$t>(maxT-aveOver),]$catchSurv)

  # the variance
  if (runNC) varSSB1   = var(p1D[p1D$t>(maxT-aveOver),]$fishDensSurv)
  if (runDC) varSSB2   = var(p2D[p2D$t>(maxT-aveOver),]$fishDensSurv)
  else if (run5YC) varSSB2   = var(p2D[p2D$t>(maxT-aveOver),]$fishDensSurv)
  if (runMPA) varSSB3   = var(p3D[p3D$t>(maxT-aveOver),]$fishDensSurv)
  if (runNC) varCatch1 = var(p1D[p1D$t>(maxT-aveOver),]$catchSurv)
  if (runDC) varCatch2 = var(p2D[p2D$t>(maxT-aveOver),]$catchSurv)
  else if (run5YC) varCatch2 = var(p2D[p2D$t>(maxT-aveOver),]$catchSurv)
  if (runMPA) varCatch3 = var(p3D[p3D$t>(maxT-aveOver),]$catchSurv)

  # size
  if (runNC) size1   = mean(p1D[p1D$t>(maxT-aveOver),]$meanSizeSurv)
  if (runDC) size2   = mean(p2D[p2D$t>(maxT-aveOver),]$meanSizeSurv)
  else if (run5YC) size2   = mean(p2D[p2D$t>(maxT-aveOver),]$meanSizeSurv)
  if (runMPA) size3   = mean(p3D[p3D$t>(maxT-aveOver),]$meanSizeSurv)
  if (runNC) sizeSSB1 = mean(p1D[p1D$t>(maxT-aveOver),]$meanSizeSSBSurv)
  if (runDC) sizeSSB2 = mean(p2D[p2D$t>(maxT-aveOver),]$meanSizeSSBSurv)
  else if (run5YC) sizeSSB2 = mean(p2D[p2D$t>(maxT-aveOver),]$meanSizeSSBSurv)
  if (runMPA) sizeSSB3 = mean(p3D[p3D$t>(maxT-aveOver),]$meanSizeSSBSurv)

  # fecundity
  if (runNC) meanFec1   = mean(p1D[p1D$t>(maxT-aveOver),]$meanFecSurv)
  if (runDC) meanFec2   = mean(p2D[p2D$t>(maxT-aveOver),]$meanFecSurv)
  else if (run5YC) meanFec2   = mean(p2D[p2D$t>(maxT-aveOver),]$meanFecSurv)
  if (runMPA) meanFec3   = mean(p3D[p3D$t>(maxT-aveOver),]$meanFecSurv)
  if (runNC) totalFec1 = mean(p1D[p1D$t>(maxT-aveOver),]$totalFecSurv)
  if (runDC) totalFec2 = mean(p2D[p2D$t>(maxT-aveOver),]$totalFecSurv)
  else if (run5YC) totalFec2 = mean(p2D[p2D$t>(maxT-aveOver),]$totalFecSurv)
  if (runMPA) totalFec3 = mean(p3D[p3D$t>(maxT-aveOver),]$totalFecSurv)

  
  
  #  if (meanSSB1 > 0) {
#    ssbProp2   = meanSSB2 / meanSSB1
#    catchProp2 = meanCatch2 / meanCatch1
#    ssbProp3   = meanSSB3 / meanSSB1
#    catchProp3 = meanCatch3 / meanCatch1
#  } else {
#    myP("No surviving stock.")
#	ssbProp2 = catchProp2 = ssbProp3 = catchProp3 = 0
#  }
  
  # get recovery percentage
  if (runNC) recProp1 = mean(p1D[p1D$t==recPoint,]$fishDensSurv) / meanSSB1
  if (runDC) recProp2 = mean(p2D[p2D$t==recPoint,]$fishDensSurv) / meanSSB2
  else if (run5YC) recProp2 = mean(p2D[p2D$t==recPoint,]$fishDensSurv) / meanSSB2
  if (runMPA) recProp3 = mean(p3D[p3D$t==recPoint,]$fishDensSurv) / meanSSB3
  
  # get discounted catch 
  calcDiscCatch <- function(theD, discRate=discRate, endPoint=20) {
    totalCatchValue = 0
    for (t in 1:endPoint) {
	  tCatch = mean(theD[theD$t==t,]$catchSurv)
	  totalCatchValue = totalCatchValue + tCatch*discRate^(t-1)
	}
	return(totalCatchValue)
  }
  
  if (runNC) discCatch1 = calcDiscCatch(p1D, discRate)
  if (runDC) discCatch2 = calcDiscCatch(p2D, discRate)
  else if (run5YC) discCatch2 = calcDiscCatch(p2D, discRate)
  if (runMPA) discCatch3 = calcDiscCatch(p3D, discRate)
  # use as comparison for other two
#  if (discCatch1 > 0) {
#    discCatch2 = calcDiscCatch(p2D) / discCatch1
#    discCatch3 = calcDiscCatch(p3D) / discCatch1
#  } else {
#    myP("No catch.")
#	discCatch2 = discCatch3 = 0
#  }
  
  if (!runNC) {
    propExtinct1 = 0
	recProp1     = 0
	discCatch1   = 0
	equilSSB1    = 0
    equilCatch1  = 0	
	meanSSB1     = 0
	meanCatch1   = 0
	varSSB1      = 0
	varCatch1    = 0
	size1        = 0
	sizeSSB1     = 0
	meanFec1     = 0
	totalFec1    = 0
  }
  if (!runDC && !run5YC) {
    propExtinct2 = 0
	recProp2     = 0
	discCatch2   = 0
	equilSSB2    = 0
    equilCatch2  = 0	
	meanSSB2     = 0
	meanCatch2   = 0
	varSSB2      = 0
	varCatch2    = 0
	size2        = 0
	sizeSSB2     = 0
	meanFec2     = 0
	totalFec2    = 0
  }
  if (!runMPA) {
    propExtinct3 = 0
	recProp3     = 0
	discCatch3   = 0
	equilSSB3    = 0
    equilCatch3  = 0	
	meanSSB3     = 0
	meanCatch3   = 0
	varSSB3      = 0
	varCatch3    = 0
	size3        = 0
	sizeSSB3     = 0
	meanFec3     = 0
	totalFec3    = 0
  }
  
  return(list(prefix=prefix, 
              propExtinct1=propExtinct1, propExtinct2=propExtinct2, propExtinct3=propExtinct3,
			  recProp1=recProp1, recProp2=recProp2, recProp3=recProp3, 
			  discCatch1=discCatch1, discCatch2=discCatch2, discCatch3=discCatch3,  
			  equilSSB1=meanSSB1, equilSSB2=meanSSB2, equilSSB3=meanSSB3,  
			  equilCatch1=meanCatch1, equilCatch2=meanCatch2, equilCatch3=meanCatch3,
			  varSSB1=varSSB1, varSSB2=varSSB2, varSSB3=varSSB3,  
			  varCatch1=varCatch1, varCatch2=varCatch2, varCatch3=varCatch3,
			  size1=size1, size2=size2, size3=size3,  
			  sizeSSB1=sizeSSB1, sizeSSB2=sizeSSB2, sizeSSB3=sizeSSB3,  
			  meanFec1=meanFec1, meanFec2=meanFec2, meanFec3=meanFec3,  
			  totalFec1=totalFec1, totalFec2=totalFec2, totalFec3=totalFec3))

}

showDensityMetrics <- function(pL=NULL, prefix="", recPoint=20, manV=c(MAN_NC, MAN_DC, MAN_MPA), printMetrics=TRUE) {
  if (is.null(pL)) {
	fileStr = combine(prefix, "ts")
    pD = loadMatrixFromCSVFile(fileStr)
    pL = calcDensityMetrics(pD, prefix, manV=manV)
  }
  if (printMetrics) {
    metricV = c("Metric", "Ext.prop.", paste("Rec. at", recPoint), paste("Catch by", recPoint), "Equil. dens.", "Equil. catch")
    metricD = as.data.frame(metricV)
    metricD$NC  = c(manV[1], pL$propExtinct1, pL$recProp1, pL$discCatch1, pL$equilSSB1, pL$equilCatch1)
    metricD$DC  = c(manV[2], pL$propExtinct2, pL$recProp2, pL$discCatch2, pL$equilSSB2, pL$equilCatch2)
    metricD$MPA = c(manV[3], pL$propExtinct3, pL$recProp3, pL$discCatch3, pL$equilSSB3, pL$equilCatch3)
  
    myP("Scenario:", prefix)
    print(metricD)
  }
  return(pL)
}

runSizeMovieModel <- function(cW=CLASS_WIDTH, hProp=0.1, mpaProp=MPA_PROP, numReplic=5, maxT=300, mpaSize=MPA_SIZE, isAllee=TRUE, addCat=TRUE, prefix="", dropSmallest=TRUE) {
  pL = getParamList(hProp=hProp, mpaProp=mpaProp, mpaSize=mpaSize, isAllee=isAllee, maxT=maxT, numReplic=numReplic, classWidth=cW, addCat=addCat);

  resL  = runModelFromList(pL, showPlots=FALSE, saveToFile=FALSE, prefix=prefix);
  meanD = calcSpatialDensity(resL$means, resL$res, pL);
  meanD$ssbDensity = meanD$meanSSB / (TOTAL_AREA*10000);
  meanD$catchHa    = meanD$meanCatch / (TOTAL_AREA);

  # draw SSB density plot
  thePlot = ggplot(meanD[meanD$t>=(maxT/2 - 5),], aes(x=t)) + 
            geom_line(aes(y=resDens, linetype="Reserve")) +
			geom_line(aes(y=fishDens, linetype="Fishery"));
			#+ scale_linetype_manual(values=c("Reserve"="solid";
  thePlot = thePlot + labs(title=paste(prefix, "SSB density"), x="Time", y="Density (#/m^2)");
  myPPlot(thePlot);
  
  # get mean size data within reserves and outside
  # dim: rep, size, time, block
  sizeDataA = resL$res$compSizeDataA
  blockEffortV = getBlockEffortVector(pL$mpaSize, pL$mpaProp, pL$hProp);
  fishV = blockEffortV>0;
  resV  = !fishV; 

  # sum across relevant blocks and replicates
  # dim: size, time
  if (numReplic>1) sumColV = c(2,3)
  else sumColV = c(1,2); # deal with dimension loss for numReplic=1
  resSizeM  = apply(sizeDataA[,,,resV], sumColV, sum) / (sum(resV)*numReplic*10000);
  fishSizeM = apply(sizeDataA[,,,fishV], sumColV, sum) / (sum(fishV)*numReplic*10000);

  # reform matrices as melted data
  resSizeD  = melt(resSizeM, varnames=c("size", "t"), value.name="density");
  fishSizeD = melt(fishSizeM, varnames=c("size", "t"), value.name="density");

  # dropping smallest size class
  
  if (dropSmallest) {
    resSizeD  = resSizeD[resSizeD$size>1,]; 
    fishSizeD = fishSizeD[fishSizeD$size>1,]; 
  }	
  # filling in size value
  resSizeD$size  = (resSizeD$size-1)*CLASS_WIDTH + CLASS_WIDTH/2;
  fishSizeD$size = (fishSizeD$size-1)*CLASS_WIDTH + CLASS_WIDTH/2;;
  resSizeD$Type  = "Reserve";
  fishSizeD$Type = "Fishery";

  # merge together
  allSizeD = rbind(resSizeD, fishSizeD);
  
  # create cat vector
  eventA = resL$res$eventA;
  catV   = apply(eventA, 3, mean);
  
  fileStr = "sizeDens";
  if (prefix!="") fileStr = combine(prefix, fileStr);
  if (maxT>50) {
    catT = maxT/2;
    tV   = c(catT-10, catT+30);
  } else tV = NULL;
  createSizeMovie(allSizeD, meanD, catV, tV, fileStr=fileStr, prefix=prefix);
}

########################
##  MSY Calculation
########################

contMSY <- function(paramStr=PARAM_K, numReplic=20, maxT=100, hInc=0.01, minH=0, maxH=0.25, contV=NULL) {
  if (paramStr==PARAM_AGG) {
    # Minimum 0.5, or population extinct
	# Include -1, meaning no Allee effect
    if (is.null(contV)) contV = c(seq(0.5, 2.0, by=0.1),-1) # aggA mult
	manV = c(MAN_NC, MAN_MPA)
  } else if (paramStr==PARAM_K) {
    if (is.null(contV)) contV = seq(0.2, 2.0, by=0.1) # blockK mult 
	manV = c(MAN_NC, MAN_MPA)
  } else if (paramStr==PARAM_DISP) {
    if (is.null(contV)) contV = seq(0.25, 2.0, by=0.25) # dispersal mult 
    manV = c(MAN_MPA)  # run MPA only
  } else if (paramStr==PARAM_MPA_W) {
    if (is.null(contV)) contV = c(1, 2, 5, 10, 15, 20) # MPA width 
    manV = c(MAN_MPA)  # run MPA only
  } else if (paramStr==PARAM_MPA_P) {
    if (is.null(contV)) contV = seq(0, 0.8, by=0.1) # MPA proportion 
    manV = c(MAN_MPA)  # run MPA only
  }

  contAgg = AGG_A
  contK   = DEF_K
  contDisp = 1.0
  contMPASize = MPA_SIZE
  contMPAProp = MPA_PROP
  for (contVal in contV) {
    if (paramStr==PARAM_AGG) contAgg = AGG_A*contVal
	else if (paramStr==PARAM_K) contK = DEF_K*contVal
	else if (paramStr==PARAM_DISP) contDisp = contVal
	else if (paramStr==PARAM_MPA_W) contMPASize = contVal
	else if (paramStr==PARAM_MPA_P) contMPAProp = contVal
	
	prefix = combine(paramStr, contVal)
	getMSYChart(numReplic, maxT, hInc, minH, maxH, manV=manV,
	            blockK=contK, aggA=contAgg, dispMult=contDisp, mpaSize=contMPASize, mpaProp=contMPAProp,
				prefix=prefix, saveToFile=TRUE, showPlot=FALSE, savePlot=FALSE)
  }
  
  # Aggregate files
  mergeMSYFiles(prefix=paramStr, contV=contV)
}

getMSYChart <- function(numReplic=20, maxT=100, hInc=0.01, minH=0, maxH=0.25, 
						blockK=DEF_K, aggA=AGG_A, dispMult=1.0, mpaSize=MPA_SIZE,
						manV=c(MAN_NC, MAN_MPA), mpaProp=MPA_PROP,
						saveToFile=TRUE, prefix="", showPlot=FALSE, savePlot=FALSE) {
  runNC  = MAN_NC %in% manV
  runMPA = MAN_MPA %in% manV

  # assume no Allee if negative
  isAllee = aggA>0
  
  # assumes no catastrophe, so shorter time to equilibrium
  if (runNC) pL1 = getParamList(mpaProp=0, mpaSize=mpaSize, blockK=blockK, aggA=aggA, isAllee=isAllee, maxT=maxT, numReplic=numReplic, addCat=FALSE)
  if (runMPA) pL2 = getParamList(mpaProp=mpaProp, mpaSize=mpaSize, blockK=blockK, aggA=aggA, isAllee=isAllee, maxT=maxT, numReplic=numReplic, addCat=FALSE)
  
  set.seed(RAND_SEED)
  
  # vary dispersal here
  dispSDM  = getRandomDispSDMatrix(maxT, numReplic, rType="gamma", shape=DISP_GAMMA_SHAPE*dispMult)
  randRecM = getRandomRecruitMatrix(maxT, numReplic, rType="log-normal", shape=REC_LN_SDLOG);
  
  if (runNC) res1D = testHValues(pL1, drawPlot=FALSE, hInc=hInc, minH=minH, maxH=maxH, dispSDM=dispSDM, randRecM=randRecM, label=MAN_NC, prefix=prefix)
  if (runMPA) res2D = testHValues(pL2, drawPlot=FALSE, hInc=hInc, minH=minH, maxH=maxH, dispSDM=dispSDM, randRecM=randRecM, label=MAN_MPA, prefix=prefix)
  
  if (!runNC) msyD = res2D
  else if (!runMPA) msyD = res1D
  else msyD  = rbind(res1D, res2D)
  if (showPlot || savePlot) {
    thePlot = ggplot(msyD, aes(x=h)) + theme_bw() + 
              geom_line(aes(y=msy, group=name, color=name, linetype=name), size=2) +
              geom_line(aes(y=abun, group=name, color=name, linetype=name), size=2) +
	          labs(x="Harvest proportion", y="Equilibrium catch (kg)") + 
  		 	  theme(legend.text=element_text(size=rel(1.25)), 
	              legend.title=element_text(size=0),
        	      legend.position=c(.90, .8)) +
              theme(axis.text.x=element_text(size=rel(1.25)), 
                   axis.text.y=element_text(size=rel(1.25)),
 	               axis.title.x=element_text(size=rel(1.25)),
 			       axis.title.y=element_text(size=rel(1.25)))
    if (showPlot) myPPlot(thePlot)
  }
  if (saveToFile) {
	fileName = combine(prefix, "msy")
	saveMatrixToCSVFile(fileName, msyD)
   	if (savePlot) {
	  fileName = paste("output/", combine(prefix, "msy.png"), sep="")
      myPTiff(thePlot, fileName, h=5, w=7)      
	}
  }
}

testHValues <- function(paramL, drawPlot=FALSE, aveOver=40, hInc=0.05, minH=0, maxH=0.25, dispSDM=NULL, randRecM=NULL, label="", prefix="") {
  hV = seq(minH,maxH,hInc)
  cV = c()
  aV = c()
  lastT = paramL$maxT
  aveV = (lastT-aveOver+1):lastT
  if (is.null(dispSDM)) dispSDM = getRandomDispSDMatrix(lastT, paramL$numReplic, rType="gamma")
  if (is.null(randRecM)) randRecM = getRandomRecruitMatrix(lastT, paramL$numReplic, rType="log-normal", shape=REC_LN_SDLOG)
  for (h in hV) {
    hParamL         = paramL
	hParamL$hProp   = h
#	hParamL$hPreSim = h
    meanD = runModelFromList(hParamL, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine(combine(prefix, "test"),h))$means
	# get final catch value (averaged over final period)
	cV = c(cV, mean(meanD$meanCatch[aveV]))
	# get final abundance value
	aV = c(aV, mean(meanD$meanSSB[aveV]))
  }

  # draw total catch plot
  if (drawPlot) {
	catchD  = data.frame(h=hV, catch=cV)
	thePlot = ggplot(catchD, aes(x=h)) + 
				geom_line(aes(y=catch))
	thePlot = thePlot + labs(title="MSY Catch", x="h", y="Kg")
	myPPlot(thePlot)
  }
  return(data.frame(h=hV, msy=cV, abun=aV, name=label))
}

recurseFindMSY <- function(paramL, oldH, hInc=0.01, oldCatch=-1, dispSDM=NULL, randRecM=NULL, aveOver=40) {
  if (is.null(dispSDM)) dispSDM = getRandomDispSDMatrix(paramL$maxT, paramL$numReplic, rType="gamma");
  if (is.null(randRecM)) randRecM = getRandomRecruitMatrix(paramL$maxT, paramL$numReplic, rType="log-normal", shape=REC_LN_SDLOG);

  lastT     = paramL$maxT;
  aveV      = (lastT-aveOver+1):lastT;
  firstTest = (oldCatch==-1)
  
  if (firstTest) {
    myP("Generating starting value for h=", oldH)
    hParamL         = paramL;
    hParamL$hProp   = oldH;
    hParamL$hPreSim = oldH;
    meanD = runModelFromList(hParamL, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine("test",oldH))$means;
    oldCatch = mean(meanD$meanCatch[aveV]);
  } 
    
  newH = oldH + hInc;
  if ((newH<0) || (newH>1)) return(list(msyH=oldH, msyCatch=oldCatch))

  hParamL         = paramL;
  hParamL$hProp   = newH;
  hParamL$hPreSim = newH;
  myP("Running test catch for h=", newH)
  meanD = runModelFromList(hParamL, showPlots=FALSE, saveToFile=FALSE, dispSDM=dispSDM, randRecM=randRecM, prefix=combine("test",newH))$means;
  newCatch = mean(meanD$meanCatch[aveV]);
  
  if (newCatch>oldCatch) {
    # recurse on next value (either up 1 or down 1)
	# make sure to pass on hInc since it could be negative now
	myP("New catch greater. Recursing...")
    return(recurseFindMSY(paramL, oldH=newH, hInc=hInc, oldCatch=newCatch, dispSDM=dispSDM, randRecM=randRecM))
  } else if (newCatch<=oldCatch) {
    # don't increase any more
	if (firstTest) {
  	  # if first, then go down.  Otherwise, we came from down
	  myP("New catch smaller. Decrementing...")
	  return(recurseFindMSY(paramL, oldH=oldH, hInc=-0.01, oldCatch=oldCatch, dispSDM=dispSDM, randRecM=randRecM))
	} else {
	  # not the first, so we've already looked above and come from or looked at below
	  myP("New catch smaller. Stopping...")
      return(list(msyH=oldH, msyCatch=oldCatch))
    }
  }
}

mergeMSYFiles <- function(prefix=PARAM_AGG, contV=seq(0.5, 2.0, by=0.1)) {
  ncM   = NULL
  mpaM  = NULL
  ncPM  = NULL
  mpaPM = NULL

  for (i in 1:length(contV)) {
    cVal        = contV[i]
    contFileStr = combine(combine(prefix, cVal), "msy")
#    dirStr      = paste("msy/", prefix, sep=""), "files", sep=" ")
    msyFileD    = loadMatrixFromCSVFile(contFileStr)
	
	# load NC and MPA values
	manV   = unique(msyFileD$name)
	hasNC  = MAN_NC %in% manV
	hasMPA = MAN_MPA %in% manV
    
    # Get values
	hV   = unique(msyFileD$h)
	if (hasNC) {
	  ncV = msyFileD[msyFileD$name==MAN_NC,]$msy
	  
      if (is.null(ncM)) {
        ncM       = matrix(nrow=length(hV)+1, ncol=length(contV)+1)
        ncPM      = matrix(nrow=length(hV), ncol=length(contV))
  	    ncM[-1,1] = hV
	    ncM[1,-1] = contV
	  }
      ncM[-1,i+1] = ncV
	  ncPM[,i]    = ncV / max(ncV)
	}
	if (hasMPA) {
	  mpaV = msyFileD[msyFileD$name==MAN_MPA,]$msy

      if (is.null(mpaM)) {
        mpaM       = matrix(nrow=length(hV)+1, ncol=length(contV)+1)
        mpaPM      = matrix(nrow=length(hV), ncol=length(contV))
  	    mpaM[-1,1] = hV
  	    mpaM[1,-1] = contV
	  }
	  mpaM[-1,i+1] = mpaV
	  mpaPM[,i]    = mpaV / max(mpaV)
    }
  }	  
  
  if (hasNC) {
    ncD  = melt(ncM[-1,-1], varnames=c("h", prefix), value.name="nc_yield")
    ncPD = melt(ncPM, varnames=c("h", prefix), value.name="nc_prop_msy")
  }
  
  if (hasMPA) {
    mpaD  = melt(mpaM[-1,-1], varnames=c("h", prefix), value.name="mpa_yield")
    mpaPD = melt(mpaPM, varnames=c("h", prefix), value.name="mpa_prop_msy")
  }
  if (!hasNC) msyD = merge(mpaD, mpaPD)
  else if (!hasMPA) msyD = merge(ncD, ncPD)
  else {
    msyD   = merge(ncD, mpaD)
    propD  = merge(ncPD, mpaPD)
    msyD   = merge(msyD, propD)
  }
  msyD$h         = hV[msyD$h]
  msyD[[prefix]] = contV[msyD[[prefix]]]
	 
  if (hasNC) saveMatrixToCSVFile(combine(prefix, "NC_msy"), ncM) 
  if (hasMPA) saveMatrixToCSVFile(combine(prefix, "MPA_msy"), mpaM) 
  saveMatrixToCSVFile(combine(prefix, "allDF_msy"), msyD) 
}

# Function which determines a harvest level for given strategy and proportion of MSY
# Changed to default to proportion of fMSY instead of proportion of MSY itself
getMSYValue <- function(paramStr="Base", paramVal=1.0, manType=MAN_NC, propMSY=1.0, calcPreciseMSY=TRUE, useFMSY=TRUE) {
  # select appropriate MSY file
  fileStr = combine(combine(paramStr, manType), "msy")
  msyD    = loadMatrixFromCSVFile(fileStr, dirName="msy")
  
  # All negative values map to -1
  if (paramVal<0) paramVal = -1
  
  if (paramStr!="Base") {
    # only consider the 'h' column and the appropriate parameter column
	# XXX '-' not allowed in column names, replaced by '.'
	# This could cause an error, but seems ok right now (real fractions are "0.1")
	if (paramVal==-1) paramVal = ".1"
    paramIndex = which(colnames(msyD)==combine(paramStr, paramVal))
    msyD = msyD[,c(1,paramIndex)] 
  }

  msyD$prop = msyD[,2] / max(msyD[,2])
  
  if (useFMSY) {
    # just use proportion of fMSY h value
	maxIndex = which(msyD$prop==1.0)
	return(msyD$h[maxIndex]*propMSY)
  }
  
  if (propMSY==1.0) {
    maxIndex = which(msyD$prop==1.0)
	return(msyD$h[maxIndex])
  }
  
  # search through to find the transition
  transIndex = -1
  for (i in 1:(nrow(msyD)-1)) {
    if (msyD$prop[i]<=propMSY) {
      if (msyD$prop[i+1]>propMSY) {
	    transIndex = i
		break
	  }
    }	
  }  
  
  if (calcPreciseMSY) {
    lowH    = msyD$h[transIndex]
    highH   = msyD$h[transIndex+1]
	lowMSY  = msyD$prop[transIndex]
	highMSY = msyD$prop[transIndex+1]
	
	diffH    = highH - lowH
	diffMSY  = highMSY - lowMSY 
	propBtwn = (propMSY-lowMSY) / diffMSY
	preciseH = lowH + diffH*propBtwn
	return(preciseH) 
  } else return(msyD$h[transIndex])

  
}

#####################################
#           RUN MODEL
#####################################
# OA Effect V = c(mult for density, mult for larval survival, mult for mean growth)

getParamList <- function(mpaSize=5, mpaProp=0.2, hProp=0.1, hPresim=H_PRESIM, blockK=DEF_K,
                     maxT=T_MAX, numReplic=NUM_REPLICATIONS, classWidth=CLASS_WIDTH, closeThresh=DYNAMIC_THRESH,
					 isAllee=TRUE, aggA=AGG_A, randRec=FALSE, oaEffectV=rep(1,length(oaStrV)),
					 addCat=TRUE, catMort=DEF_CAT_MORT, noHarvPostCat=-1, numCat=1, catInt=10,
					 hSDProp=0, classicSurv=FALSE) {
  return(list(mpaSize=mpaSize, mpaProp=mpaProp, hProp=hProp, hPreSim=hPresim, blockK=blockK, closeThresh=closeThresh,
              maxT=maxT, numReplic=numReplic, classWidth=classWidth, isAllee=isAllee, aggA=aggA, randRec=randRec, 
			  oaEffectV=oaEffectV, addCat=addCat, catMort=catMort, noHarvPostCat=noHarvPostCat, numCat=numCat, catInt=catInt,
			  hSDProp=hSDProp, classicSurv=classicSurv));  
}

runModelFromList <- function(paramL, dispSDM=NULL, randRecM=NULL, runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE, showPlots=FALSE, saveToFile=TRUE, createMovie=FALSE, prefix="") {
  runModel(mpaSize=paramL$mpaSize, mpaProp=paramL$mpaProp, hProp=paramL$hProp, hPreSim=paramL$hPreSim, blockK=paramL$blockK, closeThresh=paramL$closeThresh,
           maxT=paramL$maxT, numReplic=paramL$numReplic, classWidth=paramL$classWidth, isAllee=paramL$isAllee, aggA=paramL$aggA, randRec=paramL$randRec,
		   oaEffectV=paramL$oaEffectV, addCat=paramL$addCat, catMort=paramL$catMort, noHarvPostCat=paramL$noHarvPostCat, numCat=paramL$numCat, catInt=paramL$catInt,
		   dispSDM=dispSDM, randRecM=randRecM, hSDProp=paramL$hSDProp, classicSurv=paramL$classicSurv, 
		   runPreSim=runPreSim, showPreSim=showPreSim, runPreSimOnly=runPreSimOnly, showPlots=showPlots, saveToFile=saveToFile, createMovie=createMovie, prefix=prefix);
}

runModel <- function(mpaSize=5, mpaProp=0.2, hProp=0.05, hPreSim=H_PRESIM, blockK=DEF_K,
                     maxT=T_MAX, numReplic=NUM_REPLICATIONS, classWidth=CLASS_WIDTH,
					 isAllee=TRUE, aggA=AGG_A, randRec=FALSE, oaEffectV=rep(1,length(oaStrV)),
					 addCat=TRUE, catMort=DEF_CAT_MORT, noHarvPostCat=-1, closeThresh=DYNAMIC_THRESH, numCat=1, catInt=10,
					 hSDProp=0, dispSDM=NULL, randRecM=NULL, classicSurv=FALSE,
					 runPreSim=TRUE, showPreSim=FALSE, runPreSimOnly=FALSE,
					 showPlots=TRUE, saveToFile=FALSE, createMovie=FALSE, prefix="") {
  initializeGlobalValues(classWidth)

  if (runPreSim) {
    myP("Running preliminary simulation...")
    myP("Replicates:", numReplic, "Time:", T_PRESIM)
	myP("Harvest rate (after 2/3 time):", hPreSim)
	preSimAbunV = runPreSimulation(numReplic, hPreSim, (showPreSim || runPreSimOnly), prefix, blockK, 
	                               isAllee, aggA, classicSurv=classicSurv)
  } else {
    myP("Skipping preliminary simulation...")
	myP("Setting initial abundance to", INIT_ABUN_NO_PRESIM,"per size...")
	preSimAbunV = rep(INIT_ABUN_NO_PRESIM, NUM_CLASSES)
  }
  
  if (!runPreSimOnly) {
    myP("Running simulation...")
    myP("Replicates:", numReplic, "Time:", maxT)
    myP("Proportion of coastline protected:", mpaProp)
    myP("Reserve block size:", mpaSize)
    myP("Harvest proportion:", hProp, "with SD proportion", hSDProp)
	if (addCat) {
	  if (numCat==1) myP("Simulating", catMort, "mortality catastrophe...")
	  else myP("Simulating", numCat, catMort, "mortality catastrophes", catInt, "years apart...")
	}
    modelResL = runSimulation(mpaSize, mpaProp, hProp, preSimAbunV, maxT, numReplic, blockK, isAllee, aggA, randRec, oaEffectV,
							  addCat, catMort, noHarvPostCat, closeThresh, numCat, catInt, hSDProp, dispSDM, randRecM, classicSurv)
  
    # save size structure matrix
#    compSizeDataA   = array(0, c(grpReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
    # sum across blocks
	compSizeRepA = apply(modelResL$compSizeDataA, c(1,2,3), sum)
    # ave across replicates	
	compSizeM    = apply(compSizeRepA, c(2,3), mean)
	saveMatrixToCSVFile(combine(prefix, "size"), compSizeM)
  
    # plot all graphs and save the matrices
    plotResults(modelResL, showPlots, saveToFile, createMovie, prefix)
	
	return(list(means=getMeanResults(modelResL, isAllee, prefix), res=modelResL))
  }
}


# Run replicated pre-simulation with no spatial structure
# Returns a vector of mean abundance for each size
runPreSimulation <- function(numReplic=NUM_REPLICATIONS, hProp=H_PRESIM, showPreSim=FALSE, prefix="", 
                             recK=5*10^7, isAllee=TRUE, aggA=AGG_A, createMovie=FALSE, dropYOY=TRUE, classicSurv=FALSE) {
  # XXX there is no randomness in the pre-simulation, so force only a single replicate for the moment
  numReplic = 1

  # assume 2/3 no harvest, then 1/3 harvest
  # Note: this is temporal effort, not spatial
  effortPreV = c(rep(0,2*T_PRESIM/3), rep(hProp,T_PRESIM/3)); # Harvest effort during pre-simulation
  densPreA   = array(0, c(NUM_CLASSES,T_PRESIM,numReplic)); # Pre-simulation density
  # changed to include all classes, even those with zero catch
  catchPreA  = array(0, c(NUM_CLASSES,T_PRESIM,numReplic)); # Pre-simulation catch

  # Set initial abundance across replicates
  densPreA[,1,] = rep(INIT_ABUN_PRESIM, NUM_CLASSES);           

  # Adjust by 1/2 class_width to get lower limit
  lowerSizeV = gMeanSizeV - CLASS_WIDTH/2;
  P = getTransitionMatrix(lowerSizeV);        # Growth matrix
  S = getSurvivalMatrix(lowerSizeV, P, classicSurv);
  if (classicSurv) PSMatrix = P%*%S  # in classic, survival is a single value based on initial length
  else PSMatrix = P*S;  # in modern version, survival is based on both initial and final length
  
  # put these in global to compare with other P,S
  assign("newP", P, envir=.GlobalEnv);
  assign("newS", S, envir=.GlobalEnv);
  
  # H must be calculated for each effort value
  uniqEffortV = unique(effortPreV);
  # load a vector with each unique H matrix
  harvA = array(dim=c(NUM_CLASSES, NUM_CLASSES, length(uniqEffortV)));
  for (i in 1:length(uniqEffortV)) harvA[, , i] = getHarvestMatrix(uniqEffortV[i]);

  # loop through replicates
  for (i in 1:numReplic) {
    myP("Running replicate", i, "...");
    # loop through presim time
    for (tp in 1:(T_PRESIM-1)) {
	  if (tp%%10==0) myP("Time step t =", tp, "...");
	  # pull H from precalculated matrices
      H = harvA[, , match(effortPreV[tp], uniqEffortV)];
	  # get max potential settlers
	  maxSettlers = sum(LARVAL_SURVIVAL*gEggsPerIndV*densPreA[,tp,i]);
	  # apply Allee effect, if any
	  if (isAllee) maxSettlers = maxSettlers*LARVAL_ALLEE_ADJ*getAlleeEffect(sum(gPropMatureV*densPreA[,tp,i]), aggA=aggA);
	  recruitV = getRecruitV(maxSettlers, classicSurv, recK);
	  # changed to use PSMatrix
	  densPreA[,tp+1,i] = PSMatrix%*%H%*%densPreA[,tp,i] + recruitV;
#	  if (tp==T_PRESIM-1) print(PSMatrix%*%H%*%densPreA[,tp,i]);
	  # add inverse of harvest survival to catch
  	  catchPreA[,tp,i]  = (diag(1,NUM_CLASSES)-H)%*%densPreA[,tp,i];
    }
    # should this be run every time? I don't think so
    H = harvA[, , match(effortPreV[T_PRESIM], uniqEffortV)];

    catchPreA[, T_PRESIM, i] = (diag(1,NUM_CLASSES)-H)%*%densPreA[,T_PRESIM,i];
  }
  # take mean for each size class across replicates
  preSimAbunV = apply(densPreA,c(1,2),mean)[,T_PRESIM];
  
  if (showPreSim) {
    myP("Summarizing presim data...")
#	myWin();
#    plotMatrixForKernel(gMeanSizeV, yV=NULL, P, title="Transition");
#	myWin();
#    plotMatrixForKernel(gMeanSizeV, yV=NULL, S, title="Survival");
    myP("Plotting presim results...")
    fileStr = "presimAbun";
    if (prefix!="") fileStr = combine(prefix, fileStr);
    titleStr = "Presim abundance";
	presimAbunM = apply(densPreA, c(2,3), sum)*BLOCK_AREA;
	# transpose to match shape of normal results
    plotAndSaveMatrix(t(presimAbunM), title=titleStr, yLab="Count", fileStr=fileStr);
    fileStr = "presimSizeDistrib";
    if (prefix!="") fileStr = combine(prefix, fileStr);
    titleStr = "Presim size density (no YOY)";
	if (dropYOY) sizeV = c(0, preSimAbunV[-1])
	else sizeV = preSimAbunV;
	# turn into density /m2
	preSimSizeDistribM = matrix((sizeV/10000), nrow=1, ncol=length(sizeV));
    plotAndSaveMatrix(preSimSizeDistribM, xV=gMeanSizeV, title=titleStr, xLab="Size", yLab="Proportion", fileStr=fileStr);
	
	if (createMovie) {
  	  # make a size-distribution movie
      sizeM = apply(densPreA, c(1,2), mean);
	  if (dropYOY) sizeM[1,] = 0;
	  # using "mass" and "block" as a name to simplify shared movie code
      sizeD = melt(sizeM, varnames=c("block", "t"), value.name="mass");
      sizeD$Type = "All";
	  sizeD$block = sizeD$block * CLASS_WIDTH/2 + MIN_SIZE;
	  sizeD$mass = sizeD$mass / 10000;

	  totalD = data.frame(t=1:T_PRESIM, mass=apply(presimAbunM, 1, mean), ssb=0, catch=0);
 	  fileStr = "preSimSize";
      if (prefix!="") fileStr = combine(prefix, fileStr);
      createSpatialMovie(sizeD, totalD, rep(1,T_PRESIM), yLab="#/m2", fileStr=fileStr, prefix=prefix);
	} else {
	  # plot sizes over time
	  # average over replicates
	  sizeM = apply(densPreA, c(1,2), mean);
	  if (dropYOY) sizeM[1,] = 0;
	  sizeD = melt(sizeM, varnames=c("size", "t"), value.name="density");
	  sizeD$size = sizeD$size * CLASS_WIDTH/2 + MIN_SIZE;
	  sizeD$density = sizeD$density / 10000;
	  thePlot = ggplot(data=sizeD, aes(x=t, y=density)) + 
	            geom_line(aes(group=size, color=size)) +
				labs(title="Size density over time", x="Time", y="#/m2");
	  myPPlot(thePlot);
	}
  }
  
  return(preSimAbunV)
}

###••••••••••••••••••••• SIMULATION MODEL•••••••••••••••••••••••••••••••••••••••••••••••

# function that returns matrix of totalAbunM, totalCatchM, totalCatchMassM                                                          
# given the spatial arrangement
# arguments are:
# I.   mpaSize (size of individual reserves [block])
# II.  protection proportion (from 0.01 to 1.0), 
# III. harvest proportion (from 0 to 1)
# Returns list with three elements:
#  -totalAbunM			Total numerical abundance
#  -totalAbunSSBM		Total number of mature individuals
#  -totalCatchM         Total catch biomass
#  ALSO RETURNS SPATIAL VALUES
#  Dimensions: row=replicate, col=time
runSimulation <- function(mpaSize, mpaProp, hProp, initAbunV, maxT=T_MAX, numReplic=NUM_REPLICATIONS, 
                          recK=DEF_K, isAllee=TRUE, aggA=AGG_A, randRec=FALSE, oaEffectV=rep(1,length(oaStrV)),
						  addCat=TRUE, catMort=DEF_CAT_MORT, noHarvPostCat=-1, closeThresh=DYNAMIC_THRESH, numCat=1, catInt=10,
						  hSDProp=0, dispSDM=NULL, randRecM=NULL, classicSurv=FALSE) { 
						  
  # break apart into groups
  # XXX currently not implemented (i.e. MAX set to 99999)
  totalReplicates = numReplic
  numGroups = totalReplicates %/% MAX_REPLICATE_GROUP                        # number of max size groups
  if (totalReplicates %% MAX_REPLICATE_GROUP > 0) numGroups = numGroups + 1  # add one for any remainder
  
  replicatesRemaining = totalReplicates
  
  # Get effort vector across blocks (h=0 for reserves)
  # Determines reserve structure
  blockEffortV   = getBlockEffortVector(mpaSize, mpaProp, hProp)
  uncertainHarv  = hSDProp>0
  dynamicClosure = (noHarvPostCat==DYNAMIC_DELAY)

  # if OA effect on growth, redo recruit distribution vector
  if (oaEffectV[OA_GROWTH_RATE_INDEX]!=1) gRecruitDistribV = getRecruitDistributionVector(FALSE, oaEffectV[OA_GROWTH_RATE_INDEX]);
  # if OA effect on maturity or fecundity, redo eggsPerIndV
  if ((oaEffectV[OA_SIZE_MAT_INDEX]!=1) || (oaEffectV[OA_FEC_INDEX]!=1)) {
    MAT_SIZE = MAT_SIZE*oaEffectV[OA_SIZE_MAT_INDEX]
    EGGS_A   = EGGS_A*oaEffectV[OA_FEC_INDEX]
    gPropMatureV = round(MAT_A/(1+exp(-(gMeanSizeV-MAT_SIZE)/MAT_B)),1)        
    assign("gPropMatureV", gPropMatureV, envir=.GlobalEnv)
    # re-initialize gEggsPerIndV
    gEggsPerIndV = 0.5*EGGS_A*wForL(gMeanSizeV)*gPropMatureV
    assign("gEggsPerIndV", gEggsPerIndV, envir=.GlobalEnv)
  }

  # Arrays for model dynamics
  # IMPORTANT: density is #/ha, not #/m^2
  densA     = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))    # Initialize Dens_SEXPL 
  catchA    = array(0, c(NUM_CLASSES,maxT,NUM_BLOCKS))    # Initialize Catch_SEXPL 
  eggsProdV = rep(0, NUM_BLOCKS)                          # Eggs produced vector
  eggsArrV  = rep(0, NUM_BLOCKS)                          # Eggs arrived vector
  eggDispM  = matrix(0, nrow=NUM_BLOCKS, ncol=NUM_BLOCKS) # Egg dispersal matrix

  # Adjust by 1/2 class_width to get lower limit
  lowerSizeV = gMeanSizeV - CLASS_WIDTH/2;
  # Potentially add OA growth and maxSize effects via multiplier
  P = getTransitionMatrix(lowerSizeV, oaGrowthMult=oaEffectV[OA_GROWTH_RATE_INDEX], oaMaxMult=oaEffectV[OA_MAX_SIZE_INDEX]);     # Growth matrix
  S = getSurvivalMatrix(lowerSizeV, P, classicSurv);

  # Potentially add OA survival effect via multiplier
  # Alter the juvenile, young adult, and mature adult portions by the OA multipliers
  S[,gMeanSizeV<50] = S[,gMeanSizeV<50]*oaEffectV[OA_JUV_SURV_INDEX]
  S[,(gMeanSizeV>=50 & gMeanSizeV<MAT_SIZE)] = S[,(gMeanSizeV>=50 & gMeanSizeV<MAT_SIZE)]*oaEffectV[OA_YA_SURV_INDEX]
  S[,gMeanSizeV>=MAT_SIZE] = S[,gMeanSizeV>=MAT_SIZE]*oaEffectV[OA_MA_SURV_INDEX]

  if (classicSurv) PSMatrix = P%*%S  # survival matrix based only on initial length
  else PSMatrix = P*S;  # survival matrix based on initial and final length

  # H must be calculated for each effort value
  uniqEffortV = unique(blockEffortV)
  # load an array with each unique H matrix
  if (!uncertainHarv) {
    harvA = array(dim=c(NUM_CLASSES, NUM_CLASSES, length(uniqEffortV)))
    for (i in 1:length(uniqEffortV)) harvA[, , i] = getHarvestMatrix(uniqEffortV[i])
  }

  # run in groups
  for (g in 1:numGroups) {
    if (replicatesRemaining>MAX_REPLICATE_GROUP) grpReplic = MAX_REPLICATE_GROUP
	else grpReplic = replicatesRemaining
    replicatesRemaining = replicatesRemaining - grpReplic
	
    # get random dispersal SD value matrix
    # shape: time x replicates
    if (is.null(dispSDM)) dispSDM = getRandomDispSDMatrix(maxT, grpReplic, rType="gamma")
    if (!is.null(randRecM)) randRec = TRUE
    else randRecM = getRandomRecruitMatrix(maxT, grpReplic, rType="log-normal")
  
    # Matrices for summed results
    # Row = replicate
    # Col = time
    totalAbunM    = matrix(0, nrow=grpReplic, ncol=maxT)  # Total abundance (individuals)
    totalAbunSSBM = matrix(0, nrow=grpReplic, ncol=maxT)  # Total abundance (breeders)
    totalCatchM   = matrix(0, nrow=grpReplic, ncol=maxT)  # Total catch (weight in tons)
  
    # Arrays for spatial results (mass)
    # dim = rep, time, blocks
    spatialMassA    = array(0, c(grpReplic,maxT,NUM_BLOCKS))  # Spatial mass (individuals)
    spatialSSBA     = array(0, c(grpReplic,maxT,NUM_BLOCKS))  # Spatial mass (breeders)
    spatialCatchA   = array(0, c(grpReplic,maxT,NUM_BLOCKS))  # Spatial catch (weight in tons)
    spatialSSBDensA = array(0, c(grpReplic,maxT,NUM_BLOCKS))  # Spatial density (breeders)
    compSizeDataA   = array(0, c(grpReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
  

    # get event array
    eventA = getEventSurvivalArray(maxT, grpReplic, catMort, numCat, catInt);
    if (!addCat) eventA[,,,] = 1;  # eliminate all catastrophes
  
    # Run simulation
    closureLengthV = c()  # total closure
    harvestableV   = gMeanSizeV>MIN_LAND_SIZE
    for (rep in 1:grpReplic) { # Start loop on replicates
      repClosureV = c()
      myP("Starting replicate", rep, "...")
      densA[,1,] = initAbunV

	  noHarvUntil  = -1  # timer to allow cessation of harvest for X years after catastrophe
	  harvThresh   = -1  # threshold to resume harvest
	  closureStart = -1  # dynamic closure start
      for (t in 1:(maxT-1)) {  # Start loop on time
	    if (t%%10==0) myP("Time step t =", t, "...")
	    # Get random dispersal SD value
	    dispSD = getDispSD(dispSDM[t,rep])
	                                 
        # Calculate egg dispersal matrix for this time step									 
        for (j in 1:(NUM_BLOCKS/2)) { 
	      eggDispM[,j] = propDispToBlock(c((1-j):(NUM_BLOCKS-j)), dispSD) + 
		             rev(propDispToBlock(c(j:(NUM_BLOCKS+(j-1))), dispSD))
	    }
        # adds circular nature to coastline
        eggDispM[1:NUM_BLOCKS,(NUM_BLOCKS/2+1):NUM_BLOCKS] = eggDispM[NUM_BLOCKS:1,(NUM_BLOCKS/2):1]

	    # apply Allee effect to each block, if any
	    # Adds OA fertilization effect via overall effect or density multiplier
	    if (isAllee) effDensM  = getEffectiveDensityMatrix(densA[,t,], aggA, oaFertMult=oaEffectV[OA_FERT_SUCCESS_INDEX], 
	                                                     oaDensityMult=oaEffectV[OA_ALLEE_SHIFT_INDEX])
	    else effDensM = densA[,t,]
        eggsProdV = apply(sweep(effDensM,gEggsPerIndV,MARGIN=1,'*'), 2, sum) # vectorized 
        eggsArrV  = apply(sweep(eggDispM,MARGIN=2,eggsProdV,'*'), 1, sum)    # vectorized
	    if (isAllee) eggsArrV = eggsArrV*LARVAL_ALLEE_ADJ

	    # XXX hard-coding the spatial aspect of the catastrophe here
	    # Rewrite later to be more flexible
	    if (addCat) {
	      # from time of catastrophe
		  wasCat = (mean(eventA[,1,t,rep])<1)
	      if (wasCat) {
		    if (dynamicClosure) {
		      # set threshold based on prior time steps
		      harvThresh = mean(densA[harvestableV,((t-10):t),]) * closeThresh
			  # keep track of closure length
			  closureStart = t+1
		    } else noHarvUntil = t + noHarvPostCat  # constant delay (none if -1)
		  }
	    } else wasCat = FALSE

        # Post-catastrophe, assess if dynamic closure continues
	    if (!wasCat&&(closureStart>0)) {
	      curDens = mean(densA[harvestableV,t,])
		  if (curDens>=harvThresh) {
		    # add closure to list
		    repClosureV = c(repClosureV, t-closureStart)
		    # end the closure
		    noHarvUntil  = -1
		    closureStart = -1
		  } else noHarvUntil = t + 1  # continue for one step
	    }
 
	    # Apply growth, survival, and harvest and add in new recruits
        for (block in 1:NUM_BLOCKS) {
	      evPS = PSMatrix
	      # apply event survival (defaults to 1 everywhere)
		  evPS = evPS*eventA[,block,t,rep]
		
		  # pull H from unique vector
		  if ((t<noHarvUntil) || (blockEffortV[block]==0)) blockH = diag(NUM_CLASSES)
		  else {
            # have to generate on the fly if uncertain
		    if (uncertainHarv) blockH = getHarvestMatrix(blockEffortV[block], hSDProp, verbose=FALSE)
		    else blockH = harvA[, , match(blockEffortV[block], uniqEffortV)]
		  }
		
		  # apply larval survival, including possible OA effect
		  numSettlers = LARVAL_SURVIVAL*eggsArrV[block]*oaEffectV[OA_LARV_SURV_INDEX]
		  # add random recruitment pre-DD
		  if (randRec) numSettlers = numSettlers*randRecM[t,rep]
		  recruitV = getRecruitV(numSettlers, classicSurv, recK)
#		  if (randRec) recruitV = recruitV*randRecM[t,rep]
          # adjust if random recruitment survival
          densA[,t+1,block] = ((evPS%*%blockH))%*%densA[,t,block] + recruitV
          # include all sizes, even non-harvested
          catchA[,t,block]  = (diag(1,NUM_CLASSES)-blockH)%*%densA[,t,block]
        }
      } # End loop on time
    
	  # Add total closure time to replicate value
	  if (length(repClosureV)>0) closureLengthV = c(closureLengthV, sum(repClosureV))
	 
      # shouldn't run every time? (moved outside time loop)
      for (block in 1:NUM_BLOCKS) {
	    if (uncertainHarv) blockH = getHarvestMatrix(blockEffortV[block], hSDProp, verbose=FALSE)
	    else blockH = harvA[, , match(blockEffortV[block], uniqEffortV)]
        catchA[,maxT,block] = (diag(1,NUM_CLASSES)-blockH)%*%densA[,maxT,block]
	  }

	  # Weight of catch [kg]
      abunWeightA  = sweep(densA, wForL(gMeanSizeV)/1000, MARGIN=1, '*') 
      ssbWeightA   = sweep(densA[,,]*gPropMatureV, wForL(gMeanSizeV)/1000, MARGIN=1, '*') 
      catchWeightA = sweep(catchA, wForL(gMeanSizeV)*W_NOSHELL/1000, MARGIN=1, '*') 

	  # mass
	  # row, col now time, block
	  spatialMassA[rep,,]  = apply(abunWeightA, c(2,3), sum)
	  spatialSSBA[rep,,]   = apply(ssbWeightA, c(2,3), sum)
      blockCatchMassM      = apply(catchWeightA, c(2,3), sum)
	  spatialCatchA[rep,,] = blockCatchMassM
	
	  # still density here
	  # row, col now time, block
      blockDensM      = apply(densA, c(2,3), sum)
      blockDensSSBM   = apply(densA[,,]*gPropMatureV, c(2,3), sum)

      totalAbunM[rep,]    = apply(blockDensM, 1, sum)*BLOCK_AREA
      totalAbunSSBM[rep,] = apply(blockDensSSBM, 1, sum)*BLOCK_AREA
      totalCatchM[rep,]   = apply(blockCatchMassM, 1, sum)*BLOCK_AREA
	  spatialSSBDensA[rep,,] = blockDensSSBM # not summed across blocks
	  compSizeDataA[rep,,,]  = densA # not summed across blocks or sizes
    } # End loop on replicates
  } # End loop on groups
  
  return(list(totalAbunM=totalAbunM, totalAbunSSBM=totalAbunSSBM, totalCatchM=totalCatchM,
			  spatialMassA=spatialMassA, spatialSSBA=spatialSSBA, spatialCatchA=spatialCatchA,
			  spatialSSBDensA=spatialSSBDensA, compSizeDataA=compSizeDataA, eventA=eventA, closureLengthV=closureLengthV))
}

plotResults <- function(modelResL, showPlots=TRUE, saveToFile=FALSE, createMovie=TRUE, prefix="") {
  # Dimensions: row=replicate, col=time
  totalAbunM    = modelResL$totalAbunM
  totalAbunSSBM = modelResL$totalAbunSSBM
  totalCatchM   = modelResL$totalCatchM

  numReplic = nrow(totalAbunM)
  # create mean/sd/95% rows for each
  myP("Summarizing replicate data...")
  myP("Plotting results...")
  fileStr = "totalAbun"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total abundance, Rep:", numReplic))
  plotAndSaveMatrix(totalAbunM, title=titleStr, xLab="Time", yLab="Count", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)

  fileStr = "totalAbunSSB"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total SSB, Rep:", numReplic))
  plotAndSaveMatrix(totalAbunSSBM, title=titleStr, xLab="Time", yLab="Count", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)

  fileStr = "totalCatchMass";
  if (prefix!="") fileStr = combine(prefix, fileStr)
  titleStr = paste(prefix, paste("Total catch, Rep:", numReplic))
  plotAndSaveMatrix(totalCatchM, title=titleStr, xLab="Time", yLab="Kg", fileStr=fileStr, 
                    showPlot=showPlots, savePlotToFile=saveToFile)

  # leave now if we're not saving a movie
  if (!createMovie) return()
  
  # create movies of spatial mass patterns
  spatialMassA    = modelResL$spatialMassA
  spatialSSBA     = modelResL$spatialSSBA
  spatialCatchA   = modelResL$spatialCatchA
  spatialSSBDensA = modelResL$spatialSSBDensA
  
  # average across replicates
  # row,col time, block
  spatialMassM    = apply(spatialMassA, c(2,3), mean)
  spatialSSBM     = apply(spatialSSBA, c(2,3), mean)
  spatialCatchM   = apply(spatialCatchA, c(2,3), mean)
  spatialSSBDensM = apply(spatialSSBDensA, c(2,3), mean)
  
  # reform matrices as melted data
  massD  = melt(spatialMassM, varnames=c("t", "block"), value.name="mass")
  ssbD   = melt(spatialSSBM, varnames=c("t", "block"), value.name="mass")
  catchD = melt(spatialCatchM, varnames=c("t", "block"), value.name="mass")
  massD$Type  = "All"
  ssbD$Type   = "SSB"
  catchD$Type = "Catch"

  # merge together
  allMassD = rbind(massD, ssbD, catchD)
  
  # keep density separate
  # (keeping name as "mass" to simplify shared movie code
  ssbDensD = melt(spatialSSBDensM, varnames=c("t", "block"), value.name="mass")
  ssbDensD$Type = "SSB"
  # convert from /ha to /m2
  ssbDensD$mass = ssbDensD$mass / 10000

  # assemble total
  totalD       = data.frame(t=1:nrow(spatialMassM))
  totalD$mass  = apply(spatialMassM, 1, sum)
  totalD$ssb   = apply(spatialSSBM, 1, sum)
  totalD$catch = apply(spatialCatchM, 1, sum)
  
  # create cat vector
  eventA = modelResL$eventA
  catV   = apply(eventA, 3, mean)
  
  fileStr = "spatialMass"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  createSpatialMovie(allMassD, totalD, catV, fileStr=fileStr, prefix=prefix)

  fileStr = "spatialDens"
  if (prefix!="") fileStr = combine(prefix, fileStr)
  createSpatialMovie(ssbDensD, totalD, catV, yLab="#/m2", fileStr=fileStr, prefix=prefix)

}

# Summarizes the result matrix, plots summary, and saves both matrices to file.
plotAndSaveMatrix <- function(resM, title="", xV=NULL, xLab="", yLab="", fileStr, 
                              sepExtinct=TRUE, showPlot=TRUE, savePlotToFile=FALSE) {
  if (is.null(xV)) xV = 1:ncol(resM)
  resD       = data.frame(t=xV)
  resD$mean  = apply(resM, 2, mean, na.rm=TRUE)
  resD$sd    = apply(resM, 2, sd, na.rm=TRUE)
  resD$q_025 = apply(resM, 2, quantile, 0.025, na.rm=TRUE)
  resD$q_975 = apply(resM, 2, quantile, 0.975, na.rm=TRUE)
  
  if (sepExtinct) {
    extinctV   = resM[,ncol(resM)] < 10
	numExtinct = sum(extinctV)
	if (numExtinct<nrow(resM)) {
      # survivors
  	  if (numExtinct==(nrow(resM)-1)) {
        resD$surv_mean  = resM[!extinctV,]
        resD$surv_sd    = NA
        resD$surv_q_025 = NA
        resD$surv_q_975 = NA
	  } else {
        resD$surv_mean  = apply(resM[!extinctV,], 2, mean, na.rm=TRUE)
        resD$surv_sd    = apply(resM[!extinctV,], 2, sd, na.rm=TRUE)
        resD$surv_q_025 = apply(resM[!extinctV,], 2, quantile, 0.025, na.rm=TRUE)
        resD$surv_q_975 = apply(resM[!extinctV,], 2, quantile, 0.975, na.rm=TRUE)
	  }
	}
	
	if (numExtinct>0) {
      # extinct
	  if (numExtinct==1) {
        resD$ext_mean  = resM[extinctV,]
        resD$ext_sd    = NA
        resD$ext_q_025 = NA
        resD$ext_q_975 = NA
	  } else {
        resD$ext_mean  = apply(resM[extinctV,], 2, mean, na.rm=TRUE)
        resD$ext_sd    = apply(resM[extinctV,], 2, sd, na.rm=TRUE)
        resD$ext_q_025 = apply(resM[extinctV,], 2, quantile, 0.025, na.rm=TRUE)
        resD$ext_q_975 = apply(resM[extinctV,], 2, quantile, 0.975, na.rm=TRUE)
	  }
	}
  }
	
  # draw total abundance plot
  thePlot = ggplot(resD, aes(x=t)) + 
            geom_line(aes(y=mean)) +
			geom_line(aes(y=q_025), linetype="dashed") +
			geom_line(aes(y=q_975), linetype="dashed");
  thePlot = thePlot + labs(title=title, x=xLab, y=yLab);
  if (savePlotToFile) {
    myP("Saving plot to file...");
  	fileName = paste("output/", paste(fileStr, ".png", sep=""), sep="");
    myPTiff(thePlot, fileName, h=5, w=7);    
  } else if (showPlot) myPPlot(thePlot);
  
  # save the raw results and the summaries  
  saveMatrixToCSVFile(fileName, resM)
  fileName = combine(fileName, "summary")
  saveMatrixToCSVFile(fileName, resD)
}

# Summarizes the result matrix, plots summary, and saves both matrices to file.
createSpatialMovie <- function(allMassD, totalD, catV, xLab="Block", yLab="Kg", fileStr, prefix="", delFrames=TRUE) {
  myP("Making movie...");
  maxVal     = max(allMassD$mass);
  maxT       = max(allMassD$t);
  maxTot     = max(totalD$mass);
  # quick loop to drop all 1 points and set others=0
  for (i in 1:maxT) catV[i] = ifelse(catV[i]<1, 0, -100);
  catD      = data.frame(t=1:maxT, cat=catV);
  
  # function to draw each frame and save as PNG
  drawTimeStep <- function(ts) {
    if (ts%%10==0) myP("Drawing frame t=", ts, "...");
	tsStr = as.character(ts);
	if (ts<10) tsStr = paste("0", tsStr, sep="");
	if (ts<100) tsStr = paste("0", tsStr, sep="");
	fileName = paste(combine(prefix, combine("timestep", tsStr)), ".png", sep="");
	png(file=fileName, width=1200, height=800);
	tsD    = allMassD[allMassD$t==ts,];
	tsPlot = ggplot(data=tsD, aes(x=block, y=mass));
	tsPlot = tsPlot + geom_line(aes(group=Type, color=Type), size=2);
	tsPlot = tsPlot + scale_y_continuous(limits=c(0, maxVal)) +
	                  scale_color_manual(values=c("All"="blue", "SSB"="green", "Catch"="red"), breaks=c("All", "SSB", "Catch"));
	tsPlot = tsPlot + labs(title=combine(prefix, paste("t=",ts,sep="")), x=xLab, y=yLab) +
	                  theme(plot.title=element_text(size=rel(5)),
		  				    legend.text=element_text(size=rel(4)), 
							legend.title=element_text(size=rel(4)),
							axis.text.x=element_text(size=rel(5)), 
	                        axis.text.y=element_text(size=rel(5)),
	   				        axis.title.x=element_text(size=rel(5)),
							axis.title.y=element_text(size=rel(5)));
    # draw sub plot
	print(tsPlot);
	vp = viewport(width=0.15, height=0.15, x=0.9, y=0.75);
	# create inset temp plot with cats
	totTSD = totalD[totalD$t<=ts,];
	totPlot = ggplot(data=totTSD, aes(x=t)) + 
	           geom_line(aes(y=mass), size=1.5, color="blue") + # mass
	           geom_line(aes(y=ssb), size=1.5, color="green") + # ssb
	           geom_line(aes(y=catch), size=1.5, color="red") + # catch
	           geom_point(data=catD, aes(x=t, y=cat), size=3, color="yellow"); # catastrophes
	totPlot = totPlot + labs(title="Total", x="time", y="kg") +
				          scale_x_continuous(limits=c(0, maxT)) +
				          scale_y_continuous(limits=c(0, maxTot)) +
	   				      theme(axis.title.x=element_text(size=rel(1)),
							    axis.title.y=element_text(size=rel(1)));
	  
    print(totPlot, vp=vp);
	dev.off();
  }

  # loop function for use in saveGIF()
  loopTimeSteps <- function() {
	lapply(1:maxT, function(t) { drawTimeStep(t); });
  }
    
  # saveGIF PATH info not working quite right, do by hand
  loopTimeSteps();
  # hardcode PATH info
  # NOTE: incorporates ALL PNG files
  #       Must delete earlier files before running again
  system(paste('"C:\\Program Files\\ImageMagick-6.9.1-Q16\\convert.exe" -delay 15 *.png', 
	           paste("C:\\Research\\Code\\abalone\\output\\", paste(fileStr, ".gif", sep=""), sep="")));
  if (delFrames) file.remove(list.files(pattern=".png"));
}

createSizeMovie <- function(allSizeD, totalD, catV, tV=NULL, xLab="Size", yLab="Density (#/m2)", fileStr, prefix="", delFrames=TRUE) {
  myP("Making movie...");
  maxVal = max(allSizeD$density);
  maxT   = max(allSizeD$t);
  if (is.null(tV)) tV = c(1, maxT);
  maxTot = max(totalD$resDens, totalD$fishDens);
  # quick loop to drop all 1 points and set others=0
  for (i in 1:maxT) catV[i] = ifelse(catV[i]<1, 0, -100);
  catD = data.frame(t=1:maxT, cat=catV);
  
  # function to draw each frame and save as PNG
  drawTimeStep <- function(ts) {
    if (ts%%10==0) myP("Drawing frame t=", ts, "...");
	tsStr = as.character(ts);
	if (ts<10) tsStr = paste("0", tsStr, sep="");
	if (ts<100) tsStr = paste("0", tsStr, sep="");
	fileName = paste(combine(prefix, combine("timestep", tsStr)), ".png", sep="");
	png(file=fileName, width=1200, height=800);
	tsD    = allSizeD[allSizeD$t==ts,];
	tsPlot = ggplot(data=tsD, aes(x=size, y=density));
	tsPlot = tsPlot + geom_line(aes(group=Type, color=Type), size=2);
	tsPlot = tsPlot + scale_y_continuous(limits=c(0, maxVal)) +
	                  scale_color_manual(values=c("Reserve"="blue", "Fishery"="red"), breaks=c("Reserve", "Fishery"));
	tsPlot = tsPlot + labs(title=combine(prefix, paste("t=",ts,sep="")), x=xLab, y=yLab) +
	                  theme(plot.title=element_text(size=rel(5)),
		  				    legend.text=element_text(size=rel(4)), 
							legend.title=element_text(size=rel(4)),
							axis.text.x=element_text(size=rel(5)), 
	                        axis.text.y=element_text(size=rel(5)),
	   				        axis.title.x=element_text(size=rel(5)),
							axis.title.y=element_text(size=rel(5)));
    # draw sub plot
	print(tsPlot);
	vp = viewport(width=0.15, height=0.15, x=0.9, y=0.75);
	# create inset temp plot with cats
	totTSD = totalD[totalD$t<=ts,];
	totPlot = ggplot(data=totTSD, aes(x=t)) + 
	           geom_line(aes(y=resDens), size=1.5, color="blue") + # mass
	           geom_line(aes(y=fishDens), size=1.5, color="red") + # catch
	           geom_point(data=catD, aes(x=t, y=cat), size=3, color="yellow"); # catastrophes
	totPlot = totPlot + labs(title="Total", x="time", y="kg") +
				          scale_x_continuous(limits=c(0, maxT)) +
				          scale_y_continuous(limits=c(0, maxTot)) +
	   				      theme(axis.title.x=element_text(size=rel(1)),
							    axis.title.y=element_text(size=rel(1)));
	  
    print(totPlot, vp=vp);
	dev.off();
  }

  # loop function for use in saveGIF()
  loopTimeSteps <- function(startT=1, endT=maxT) {
	lapply(startT:endT, function(t) { drawTimeStep(t); });
  }
    
  # saveGIF PATH info not working quite right, do by hand
  loopTimeSteps(startT=tV[1], endT=tV[2]);
  # hardcode PATH info
  # NOTE: incorporates ALL PNG files
  #       Must delete earlier files before running again
  system(paste('"C:\\Program Files\\ImageMagick-6.9.1-Q16\\convert.exe" -delay 30 *.png', 
	           paste("C:\\Research\\Code\\abalone\\output\\", paste(fileStr, ".gif", sep=""), sep="")));
  if (delFrames) file.remove(list.files(pattern=".png"));
}

# returns a data frame with the mean time series
getMeanResults <- function(modelResL, isAllee, prefix) {
  # Dimensions: row=replicate, col=time
  totalAbunM    = modelResL$totalAbunM
  totalAbunSSBM = modelResL$totalAbunSSBM
  totalCatchM   = modelResL$totalCatchM
  
  maxT          = ncol(totalAbunM)
  meanResD      = data.frame(t=1:maxT)
  meanResD$meanAbun  = apply(totalAbunM, 2, mean, na.rm=TRUE)
  meanResD$meanSSB   = apply(totalAbunSSBM, 2, mean, na.rm=TRUE)
  meanResD$meanCatch = apply(totalCatchM, 2, mean, na.rm=TRUE)
  meanResD$cumCatch  = 0
  # add cumulative catch
  for (cT in (maxT/2):maxT) {
	  lastCumCatch = meanResD[meanResD$t==(cT-1),]$cumCatch 
	  meanResD[meanResD$t==cT,]$cumCatch = lastCumCatch + meanResD[meanResD$t==cT,]$meanCatch
  }
  
  # function to find mean size
  meanSizeF <- function(xV, ssbOnly=FALSE) {
    if (ssbOnly) xV = xV*gPropMatureV

    totalSizeV = xV*gMeanSizeV
	return(sum(totalSizeV)/sum(xV))
  }
  # function to find fecundity
  fecF <- function(xV, total=FALSE) {
    eggsV = xV*gEggsPerIndV
    
	if (total) return(sum(eggsV))
	else return(sum(eggsV)/sum(xV))
  }

  # determine which have collapsed
  extinctV = totalCatchM[,maxT] < (COLLAPSE_THRESHOLD*totalCatchM[,(maxT/2)-1])

  # calculate mean size, mean SSB size, mean fecundity, total fecundity
  # array(c(grpReplic,NUM_CLASSES,maxT,NUM_BLOCKS)) # complete size/space data
  repSizeDataA         = apply(modelResL$compSizeDataA, c(1,2,3), sum, na.rm=TRUE) # sum across all blocks

  # replicate-specific values
  # for each rep, metrics at equil: abun, abun SSB, catch, size, size SSB, mean fec, total fec
  repEquilD = data.frame(rep=1:nrow(totalCatchM))
  AVE_OVER  = 10
  aveV      = (maxT-AVE_OVER+1):maxT
  repEquilD$collapse = extinctV
  repEquilD$abun     = apply(totalAbunM[,aveV], 1, mean, na.rm=TRUE)   
  repEquilD$abunSSB  = apply(totalAbunSSBM[,aveV], 1, mean, na.rm=TRUE)   
  repEquilD$catch    = apply(totalCatchM[,aveV], 1, mean, na.rm=TRUE)
  equilRepSizeM      = apply(repSizeDataA[,,aveV], c(1,2), mean, na.rm=TRUE)
  repEquilD$size     = apply(equilRepSizeM, 1, meanSizeF)  
  repEquilD$sizeSSB  = apply(equilRepSizeM, 1, meanSizeF, ssbOnly=TRUE)  
  repEquilD$meanFec  = apply(equilRepSizeM, 1, fecF)  
  repEquilD$totalFec = apply(equilRepSizeM, 1, fecF, total=TRUE)  
  fileName = combine(prefix, "repData")
  saveMatrixToCSVFile(fileName, repEquilD)

  # metrics across replicates
  # note that the averages are slightly different for size and fec. than when averaging the by-replicate values
  # average of mean sizes is different from mean size of average
  if (sum(extinctV)>0) {
    # average across extinct replicates
    if (sum(extinctV)==1) meanSizeDataExtM = repSizeDataA[extinctV,,] # don't average
	else meanSizeDataExtM = apply(repSizeDataA[extinctV,,], c(2,3), mean, na.rm=TRUE) 
	
    meanResD$meanSizeExt    = apply(meanSizeDataExtM, 2, meanSizeF)
    meanResD$meanSizeSSBExt = apply(meanSizeDataExtM, 2, meanSizeF, ssbOnly=TRUE)
    meanResD$meanFecExt     = apply(meanSizeDataExtM, 2, fecF)
    meanResD$totalFecExt    = apply(meanSizeDataExtM, 2, fecF, total=TRUE)
  } else {
    meanResD$meanSizeExt    = rep(-1, maxT)
    meanResD$meanSizeSSBExt = rep(-1, maxT)
    meanResD$meanFecExt     = rep(-1, maxT)
    meanResD$totalFecExt    = rep(-1, maxT)
  }
  if (sum(!extinctV)>0) {
    # average across non-extinct replicates
    if (sum(!extinctV)==1) meanSizeDataSurvM = repSizeDataA[!extinctV,,] # don't average
	else meanSizeDataSurvM = apply(repSizeDataA[!extinctV,,], c(2,3), mean, na.rm=TRUE) 

    meanResD$meanSizeSurv    = apply(meanSizeDataSurvM, 2, meanSizeF)
    meanResD$meanSizeSSBSurv = apply(meanSizeDataSurvM, 2, meanSizeF, ssbOnly=TRUE)
    meanResD$meanFecSurv     = apply(meanSizeDataSurvM, 2, fecF)
    meanResD$totalFecSurv    = apply(meanSizeDataSurvM, 2, fecF, total=TRUE)
  } else {
    meanResD$meanSizeSurv    = rep(-1, maxT)
    meanResD$meanSizeSSBSurv = rep(-1, maxT)
    meanResD$meanFecSurv     = rep(-1, maxT)
    meanResD$totalFecSurv    = rep(-1, maxT)
  }
  
  meanResD$name      = prefix
  meanResD$isAllee   = isAllee

  return(meanResD)
}

# adds mean spatial density in reserve and outside
# creates separate means for time series which go extinct
calcSpatialDensity <- function(theD, resL, paramL, sepExtinct=TRUE) {
  numReplic       = paramL$numReplic
  # get location of reserves
  blockEffortV    = getBlockEffortVector(paramL$mpaSize, paramL$mpaProp, paramL$hProp)
  fishV           = blockEffortV>0
  resV            = !fishV
  spatialSSBDensA = resL$spatialSSBDensA
  totalCatchM     = resL$totalCatchM
  # sum across relevant blocks and replicates
  if (numReplic>1) sumCol = 2  
  else sumCol = 1 # deal with dimension loss for numReplic=1
  
  theD$resDens  = apply(spatialSSBDensA[,,resV], sumCol, sum) / (sum(resV)*numReplic*10000)
  theD$fishDens = apply(spatialSSBDensA[,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000)

  maxT         = ncol(spatialSSBDensA)
  # Collapse criterion: <90% pre-catastrophe catch (Worm et al. 2006)
  # XXX assumes 1st catastrophe at maxT/2
  extinctV = totalCatchM[,maxT] < (COLLAPSE_THRESHOLD*totalCatchM[,(maxT/2)-1])
  
  # Abundance
  #  replicFinalV = apply(spatialSSBDensA[,maxT,], 1, sum);
  #  initialV     = apply(spatialSSBDensA[,(maxT/2)-1,], 1, sum);
  #  extinctV     = replicFinalV<(COLLAPSE_THRESHOLD*initialV);
  propExtinct  = sum(extinctV) / numReplic
  myP("Proportion extinct:", propExtinct)
  
  if (sepExtinct) {
    # separate out extinct time series from non-
    if (propExtinct>0) {
	  if (sum(extinctV)==1) sumCol = 1
	  else sumCol = 2
      theD$resDensExt      = apply(spatialSSBDensA[extinctV,,resV], sumCol, sum) / (sum(resV)*numReplic*10000*propExtinct)
      theD$fishDensExt     = apply(spatialSSBDensA[extinctV,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000*propExtinct)
      theD$resDensExtSD    = apply(spatialSSBDensA[extinctV,,resV], sumCol, sd) / (10000)
      theD$fishDensExtSD   = apply(spatialSSBDensA[extinctV,,fishV], sumCol, sd) / (10000)
      theD$resDensExtQ025  = apply(spatialSSBDensA[extinctV,,resV], sumCol, quantile, 0.025) / (10000)
      theD$fishDensExtQ025 = apply(spatialSSBDensA[extinctV,,fishV], sumCol, quantile, 0.025) / (10000)
      theD$resDensExtQ975  = apply(spatialSSBDensA[extinctV,,resV], sumCol, quantile, 0.975) / (10000)
      theD$fishDensExtQ975 = apply(spatialSSBDensA[extinctV,,fishV], sumCol, quantile, 0.975) / (10000)
	  if (sumCol==1) theD$catchExt = totalCatchM[extinctV,]
	  else theD$catchExt = apply(totalCatchM[extinctV,], sumCol, sum) / (numReplic*propExtinct)
    } else {
      myP("No extinct replicates.")
	  theD$resDensExt      = rep(-1, maxT)
	  theD$fishDensExt     = rep(-1, maxT)
      theD$resDensExtSD    = rep(-1, maxT)
      theD$fishDensExtSD   = rep(-1, maxT)
      theD$resDensExtQ025  = rep(-1, maxT)
      theD$fishDensExtQ025 = rep(-1, maxT)
      theD$resDensExtQ975  = rep(-1, maxT)
      theD$fishDensExtQ975 = rep(-1, maxT)
	  theD$catchExt        = rep(-1, maxT)
    }
    if (propExtinct<1) {
	  if (sum(!extinctV)==1) sumCol = 1
	  else sumCol = 2;
      theD$resDensSurv      = apply(spatialSSBDensA[!extinctV,,resV], sumCol, sum) / (sum(resV)*numReplic*10000*(1-propExtinct))
      theD$fishDensSurv     = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, sum) / (sum(fishV)*numReplic*10000*(1-propExtinct))
      theD$resDensSurvSD    = apply(spatialSSBDensA[!extinctV,,resV], sumCol, sd) / (10000)
      theD$fishDensSurvSD   = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, sd) / (10000)
      theD$resDensSurvQ025  = apply(spatialSSBDensA[!extinctV,,resV], sumCol, quantile, 0.025) / (10000)
      theD$fishDensSurvQ025 = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, quantile, 0.025) / (10000)
      theD$resDensSurvQ975  = apply(spatialSSBDensA[!extinctV,,resV], sumCol, quantile, 0.975) / (10000)
      theD$fishDensSurvQ975 = apply(spatialSSBDensA[!extinctV,,fishV], sumCol, quantile, 0.975) / (10000)
	  if (sumCol==1) theD$catchSurv = totalCatchM[!extinctV,]
	  else theD$catchSurv = apply(totalCatchM[!extinctV,], sumCol, sum) / (numReplic*(1-propExtinct))
    } else {
      myP("No non-extinct replicates.")
	  theD$resDensSurv      = rep(-1, maxT)
	  theD$fishDensSurv     = rep(-1, maxT)
      theD$resDensSurvSD    = rep(-1, maxT)
      theD$fishDensSurvSD   = rep(-1, maxT)
      theD$resDensSurvQ025  = rep(-1, maxT)
      theD$fishDensSurvQ025 = rep(-1, maxT)
      theD$resDensSurvQ975  = rep(-1, maxT)
      theD$fishDensSurvQ975 = rep(-1, maxT)
	  theD$catchSurv        = rep(-1, maxT)
    }
	
	cumCatchSurv = rep(0, (maxT/2)-1)
	cumCatchExt  = rep(0, (maxT/2)-1)
	for (t in (maxT/2):maxT) {
      cumCatchSurv = c(cumCatchSurv, cumCatchSurv[t-1] + theD$catchSurv[t])
      cumCatchExt  = c(cumCatchExt, cumCatchExt[t-1] + theD$catchExt[t])
	}
	theD$cumCatchSurv = cumCatchSurv
	theD$cumCatchExt  = cumCatchExt
  }
  
  # add extinction prop to data for ease of access
  theD$propExtinct = propExtinct
  return(list(densD=theD, propExtinct=propExtinct))
}

#============================================================================# 
#  KERNEL FUNCTIONS
#============================================================================# 
ab_survX <- function(x) { 
  # s = e^(a-b*log(weight))
  SURV_A = 0.6346126
  SURV_B = 0.317388
  mort   = exp(SURV_A-SURV_B*log(wForL(x))) ### ALL NAT DATA
  return(exp(-mort))
} 
ab_harvX <- function(x, hProp, minSize=MIN_LAND_SIZE, maxSize=9999) { 
  # XXX ignoring maxsize for now
  ifelse((x>=minSize), 1-hProp, 1); 
} 
#ab_growYX <- function(y,x,meanC=1, meanXC=0.7, sd=0.3) { dnorm(y, mean=meanC+meanXC*x, sd=sd); }
#ab_fecYX <- function(y,x, fec=0.75, meanC=0.8, meanXC=0.2, sd=0.35) { fec*dnorm(y, mean=meanC+meanXC*x, sd=sd); }
#ab_kernelYX <- function(y,x) { ab_survX(x)*(ab_growYX(y,x)+ab_fecYX(y,x)); }

getHarvestMatrix <- function(hProp, hSDProp=0, minSize=MIN_LAND_SIZE, maxSize=9999, verbose=TRUE) {
  if (verbose) {
    myP("Creating harvest matrix for h =", hProp);
    if (hSDProp>0) myP("Harvest is uncertain with SD proportion =", hSDProp);
  }

  # calculate uncertain hProp
  # no error if h=0
  if ((hProp==0) || (hSDProp==0)) hVal = hProp
  else {
    # draw from rand, using SD value as a proportion of the mean
	# draw from here so all size classes have same random harvest
    hVal = rnorm(1, mean=hProp, sd=hSDProp*hProp);
	if (hVal<0) hVal = 0
	else if (hVal>1) hVal = 1;
  }
  
  H = matrix(0,NUM_CLASSES,NUM_CLASSES);  # Harvest matrix (survival from harvest)
  diag(H) <- ab_harvX(gMeanSizeV, hVal);
  return(H);
}

getEventSurvivalArray <- function(maxT, numReplic=NUM_REPLICATIONS, catMort=DEF_CAT_MORT, numCat=1, catInt=10) {
  # 75% mortality seen at Isla Natividad (Micheli et al. 2012)
  eventA = array(1, dim=c(NUM_CLASSES, NUM_BLOCKS, maxT, numReplic));
  # very simple universal death at maxT/2
  # survival, so subtract catMort
  initCatT = maxT/2;
  for (i in 1:numCat) {
    catT = initCatT + (i-1)*catInt;
    eventA[,,catT,] = 1 - catMort;
  }
  return(eventA);
}

getAlleeEffect <- function(densHa, aggA=AGG_A, alleeType=ALLEE_PROB) {
  # IMPORTANT! Density comes in as #/ha, not #/m2
  densM2 = densHa / 10000;

  # age_agg = a*d + b
  aveAggSize = aggA*densM2 + AGG_B
  
  if (alleeType==ALLEE_PROB) {
    # Probabilistic Allee:
    # P_mixed(x) = 1 - 0.5^(x-1)
    alleeMult = 1 - 0.5^(aveAggSize-1)
  } else if (alleeType==ALLEE_LIN) {
    # Linear Allee:
	# Assumes 0.2 = 80% 
	# A = 4*dens if dens<0.25, else =1-1
	if (densM2<0.25) alleeMult = 4*densM2
	else alleeMult = 1
  } else if (alleeType==ALLEE_EXP) {
    # Exponential Allee:
	# Assumes 0.2 = 80% 
    # A = exp(-(3.38-agg)*2) if agg<3.38
	if (aveAggSize<3.38) alleeMult = exp(-2*(3.38-aveAggSize))
	else alleeMult = 1
  } else alleeMult = 1  # no Allee effect
  return(alleeMult);
}

# Calculates the effective density of each block utilizing the Allee effect
# Adjusts for possible OA effect on fertilization success
getEffectiveDensityMatrix <- function(localDensM, aggA=AGG_A, oaFertMult=1, oaDensityMult=1) {
  # assumes row,col = class,block
  for (i in 1:ncol(localDensM)) {
    # calculate total density of mature individuals
	# apply OA multiplier
    totalMature = sum(gPropMatureV*localDensM[,i])*oaDensityMult;
	# apply Allee effect, if any
    # IMPORTANT! Density here is in #/ha, not #/m2
	# will convert inside
    localDensM[,i] = localDensM[,i]*getAlleeEffect(totalMature, aggA)*oaFertMult;
  }
  return(localDensM);
}

# Creates a size distribution vector for year-1 recruits 
# Must include possible OA effect
getRecruitDistributionVector <- function(isConst=FALSE, oaMult=1) {
  gRecruitDistribV = rep(0, NUM_CLASSES);
  if (isConst) {
    recruitInd = 1;
    for (i in 2:NUM_CLASSES) {
      if (gMeanSizeV[i]>LARVAL_CONST_SIZE) break
      else recruitInd = i;
    }
    gRecruitDistribV[recruitInd] = 1;
  } else {
    # Using Bardos(2005) and starting with larvae all of 500um, determine share in each size class after 1 year
	lowerSizesV = gMeanSizeV - CLASS_WIDTH/2;
	# Piggy-backing on the transition matrix code here
	gRecruitDistribV = getTransitionMatrix(lowerSizesV, oaMult, recruitVectorOnly=TRUE);
  }
  return(gRecruitDistribV);
}

# Making recruit vector
# Applies density dependence on incoming settlers
getRecruitV <- function(numSettlers, classicSurv=FALSE, recK=DEF_K) { # totSettlers is the number of arrived settlers    
  # recruits = a*settlers*e^(-settlers/b)
  REC_A = 0.01; 
  
#  REC_B = 5*10^7;
  newRecruits = REC_A*numSettlers*(exp(-numSettlers/recK));
  
  # multiply by global constant distribution vector
  # for classic, all recruits to first size class
  if (classicSurv) recruitV = c(newRecruits, rep(0, NUM_CLASSES-1))
  else recruitV = newRecruits*gRecruitDistribV;
  return(recruitV);
} 

# Creates the per-block effort vector which "creates" the MPAs
getBlockEffortVector <- function(mpaSize, mpaProp, hProp) {
  # Determine block structure of reserve network
  numBlockProt = round(mpaProp*NUM_BLOCKS, 0);           # Total number of blocks protected
  numReserve   = round(numBlockProt/mpaSize, 0); # Total number of reserves
  numBlockFish = NUM_BLOCKS - numBlockProt;               # Total number of fished blocks

  x  = floor(numBlockFish/numReserve);
  y  = ceiling(numBlockFish/numReserve);
  nx = numBlockFish - x*numReserve;
  ny = numReserve - nx;

  # Reserves exist as h=0 blocks
  # generate spatial effort vector across blocks
  if (is.na(nx)|is.na(ny)) { 
    # No reserves
    blockEffortV = rep(hProp, NUM_BLOCKS); 
  } else if (nx==ny) { 
    blockEffortV = rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), ny); 
  } else if (nx>ny) { 
    blockEffortV = c(rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), ny),
	                 rep(c(rep(0,mpaSize),rep(hProp,y)), (nx-ny))); 
  } else if (nx<ny) { 
    blockEffortV = c(rep(c(rep(0,mpaSize),rep(hProp,x),rep(0,mpaSize),rep(hProp,y)), nx), 
	                 rep(c(rep(0,mpaSize),rep(hProp,x)), (ny-nx)));
  }

  return(blockEffortV);
}

# Random matrix of dispersal SD values across time and replicates
# Currently only supports gamma distribution
getRandomDispSDMatrix <- function(maxT, nReplic, rType="gamma", shape=DISP_GAMMA_SHAPE) {
  # generate random value vector
  if (rType=="gamma") randV = rgamma(maxT*nReplic, shape=shape, rate=DISP_GAMMA_RATE)
  else {
    myP("WARNING: Distribution", rType, "unsupported for getRandomDispSDMatrix.");
	randV = rep(0, maxT*nReplic);
  }
  
  # create and return matrix for time and replicates
  return(matrix(randV, nrow=maxT, ncol=nReplic));
}

# Random matrix of recruitment SD values across time and replicates
# Values vary from 
getRandomRecruitMatrix <- function(maxT, nReplic, rType="log-normal", shape=REC_LN_SDLOG) {
  # generate random value vector
  if (rType=="gamma") {
    myP("Creating random recruit matrix with gamma distribution...");
    randV = rgamma(maxT*nReplic, shape=shape, rate=REC_GAMMA_RATE);
  } else if (rType=="log-normal") {
    myP("Creating random recruit matrix with log-normal distribution...");
    randV   = rlnorm(maxT*nReplic, meanlog=REC_LN_MEANLOG, sdlog=shape);
  }	else if (rType=="good-bad") {
    myP("Creating random recruit matrix with good-bad years...");
    randV   = runif(maxT*nReplic);
	threshV = randV > REC_GB_THRESHOLD;
	randV[threshV]  = shape;
	randV[!threshV] = REC_GB_BAD;
  } else {
    myP("WARNING: Distribution", rType, "unsupported for getRandomRecruitMatrix.");
	randV = rep(0, maxT*nReplic);
  }
  
  # create and return matrix for time and replicates
  return(matrix(randV, nrow=maxT, ncol=nReplic));
}

# Creates transition matrix by integrating fitted gamma function across each size class
# Takes a size class vector of lower bounds, not midpoints
getTransitionMatrix <- function(sizesV, oaGrowthMult=1, oaMaxMult=1, recruitVectorOnly=FALSE) {
  # growth parameters from "Fitting_Bardos" MLE fit
  # load("Fitting_Bardos"); estimate5_3;
  
  #XXX test cynthia's growth function (MEPS catton)
  B_GI    = 0.5635249;
  B_SIGMA = 55.9560109;
  B_LN    = 150.3854074;
  B_GAMMA = 1.7190492;
  B_BETA  = 1.4780344;

  dT    = 1;
  m     = 3;
  n     = 3;

  # apply OA effect to growth
  B_GI = B_GI*oaGrowthMult
  # apply OA effect to max
  B_LN = B_LN*oaMaxMult
  
  # local functions to generate gamma values
  # calculate Linf (Eq. A2 in Appendix)
  LinfF <- function(dL, L1) { return(((L1+dL)*L1^(-exp(-B_GI*dT)))^(1/(1-exp(-B_GI*dT)))); }
  # nF and dF are numerator/denominator functions for lambdaF and roF
  numF <- function(L1) { return(B_LN/(1+(B_BETA*L1/B_LN)^m)); }
  denF <- function(L1) { return(B_SIGMA/(1+(B_GAMMA*L1/B_LN)^n)); }
  # gamma rate (Eq. A3 in Appendix)
  lambdaF <- function(L1) { return(numF(L1)/(denF(L1))^2); }
  # gamma shape (Eq. A4 in Appendix)
  roF <- function(L1) { return((numF(L1))^2/(denF(L1))^2); }

  # function to generate the matrix itself
  # NOTE: Needs sizes vector with one extra large size (i.e. beyond the limit)
  createTransMatrix <- function(sizesPlusV) {
    myP("Creating transition matrix...");
    numSizesM1 = length(sizesPlusV) - 1;
    P = matrix(0, nrow=numSizesM1, ncol=numSizesM1);

    # integrate k across size class
    for (i in 1:numSizesM1) {
	  for (j in 1:numSizesM1) {
	    # integrand function
	    integrand <- function(L1) { 
	      return(pgamma(LinfF(sizesPlusV[i+1]-L1,L1)-L1, rate=lambdaF(L1), shape=roF(L1)) - 
		         pgamma(LinfF(sizesPlusV[i]-L1,L1)-L1, rate=lambdaF(L1), shape=roF(L1)));
	    }
        P[i, j] = round((integrate(integrand,lower=sizesPlusV[j],upper=sizesPlusV[j+1]))$value/(sizesPlusV[j+1]-sizesPlusV[j]), 4);
  	  }
    }
    P[numSizesM1, numSizesM1] = 1;
	
	# make sure each column adds to 1
	for (col in 1:ncol(P)) P[,col] = P[,col] / sum(P[,col]);
	
    return(P);
  }

  if (recruitVectorOnly) {
    # Assume recruits start at L0 and grow for one year
    L0 = 1;
    # Using Bardos(2005) and starting with larvae all of 500um, determine share in each size class after 1 year
    recruitV = rep(0, length(sizesV));
    # convert lower sizes vector to upper and add L0
	upperSizesPlusV = c(L0, sizesV + CLASS_WIDTH);
    # integrate growth across size class
	# note that we assume that anything below MIN_SIZE is bumped up to size class 1
	for (i in 1:length(sizesV)) {
	  recruitV[i] = pgamma(LinfF(upperSizesPlusV[i+1]-L0, L0)-L0, rate=lambdaF(L0), shape=roF(L0)) - 
	                pgamma(LinfF(upperSizesPlusV[i]-L0, L0)-L0, rate=lambdaF(L0), shape=roF(L0));
	}
#	recruitV = recruitV / sum(recruitV);
	return(recruitV);
  } else {
    # create sizes vector with one extra large size
    sizesPlusV = c(sizesV, sizesV[length(sizesV)]+CLASS_WIDTH);
    # call and return
    return(createTransMatrix(sizesPlusV));
  }
}

# Creates transition matrix by integrating fitted gamma function across each size class
# Takes a size class vector of lower bounds, not midpoints
getSurvivalMatrix <- function(sizesV, pM, classicSurv=FALSE) {
  # This is more complicated now
  # Need to integrate across expected growth, or survival is way too low (due to low mean size)
  # Depends a lot on size-class
  
  if (classicSurv) {
    # generate simple length-based diagonal using class midpoint
	S = matrix(0, nrow=length(gMeanSizeV), ncol=length(gMeanSizeV));
	diag(S) <- ab_survX(gMeanSizeV);
    return(S);
  }
  
  # function to generate the matrix itself
  # NOTE: Needs sizes vector with one extra large size (i.e. beyond the limit)
  createSurvMatrix <- function(sizesPlusV) {
    myP("Creating survival matrix...");
    numSizesM1 = length(sizesPlusV) - 1;
    S = matrix(0, nrow=numSizesM1, ncol=numSizesM1);

    # integrate across size class
    for (i in 1:numSizesM1) {
	  for (j in 1:numSizesM1) {
	    if (i>=j) { # no negative growth
    	  # integrand function
		  if (i==j) midPoint1 = sizesPlusV[j]
		  else midPoint1 = sizesPlusV[j] + CLASS_WIDTH/2;
	      integrand <- function(L1) { 
		    intSum = 0;
		    dL     = L1 - midPoint1;
			if (dL<0) return(log(ab_survX(midPoint1)));
		    D_K    = 100;
		    kdL    = dL/D_K;
		    for (k in 1:D_K) intSum = intSum + log(ab_survX(midPoint1+k*kdL));
	        return(intSum/D_K);
	      }
		  ans = integrate(integrand,lower=sizesPlusV[i],upper=sizesPlusV[i+1])$value/(sizesPlusV[i+1]-sizesPlusV[i]);
#		  if (j==1) myP(ans, exp(ans));
          S[i, j] = exp(ans);
#		  S[i,j] = exp((log(ab_survX(midPoint1)) + log(ab_survX(midPoint2)))/2);
        }
  	  }
    }

    return(S)
  }

  # create sizes vector with one extra large size
  sizesPlusV = c(sizesV, sizesV[length(sizesV)]+CLASS_WIDTH)
  # call and return
  return(createSurvMatrix(sizesPlusV))
}

initializeGlobalValues()

################################
#      PLOTTING RESULTS
################################
plotGrowthLikelihood <- function(pM, lV=NULL) {
  # pMatrix: row = Lt+1, col=Lt
  if (is.null(lV)) lV = gMeanSizeV;
  growthM = matrix(0, nrow=length(lV)^2, ncol=3);
  for (i in 1:nrow(pM)) {
    Lt_1 = lV[i];
	for (j in 1:ncol(pM)) {
	  Lt = lV[j];
	  dL = Lt_1 - Lt;
	  prop = pM[i,j];
	  prop = prop / max(pM[,j]);
	  growthM[(i-1)*ncol(pM)+j,] = c(Lt, dL, prop);
	}
  }
  growthD = data.frame(l=growthM[,1], dL=growthM[,2], prop=growthM[,3]);
  growthD = growthD[growthD$prop>0.00001,];
  thePlot = ggplot(growthD, aes(x=l, y=dL)) + theme_bw() +
            geom_point(aes(color=prop), size=3) +
	        labs(x="Length (mm)", y="Growth (mm)") + 
  			theme(legend.text=element_text(size=rel(1.25)), 
	              legend.title=element_text(size=0),
        	      legend.position=c(.90, .8)) +
            theme(axis.text.x=element_text(size=rel(1.25)), 
                   axis.text.y=element_text(size=rel(1.25)),
 	               axis.title.x=element_text(size=rel(1.25)),
 			       axis.title.y=element_text(size=rel(1.25)));

  myPPlot(thePlot);
}

plotCohort <- function(numInd, pM, sM, maxT=10, showHist=FALSE) {
  indM = matrix(0, nrow=numInd, ncol=maxT);
  indM[,1] = 8; # start at ave size class
  cumPM    = pM;
  totalRec = 0;
  for (i in 2:nrow(pM)) cumPM[i,] = cumPM[i,] + cumPM[i-1,];
  for (t in 2:maxT) {
    growRandV = runif(numInd);
    survRandV = runif(numInd);
    for (i in 1:numInd) {
      indL = indM[i,t-1];
	  if (indL>0) { # alive
		if (survRandV[i]>sM[indL,indL]) indM[i,t] = 0 # dies
		else {
		  # grows
		  newIndL   = sum(cumPM[,indL]<growRandV[i])+1;
		  if (t==maxT) totalRec  = totalRec + gEggsPerIndV[newIndL];
#		  if (indL==1) myP(growRandV[i], indL, newIndL, gMeanSizeV[indL], gMeanSizeV[newIndL]);
		  indM[i,t] = newIndL;
		}
	  }
    }    
  }
  
  # count # of >0 values in each column
  survV = apply(apply(indM, 2, function(x) { ifelse(x>0,1,0);}), 2, sum) / numInd;
  survD = data.frame(t=1:maxT, surv=survV);

  indD  = melt(indM, c("ind", "t"), value.name="length");
  indD  = indD[indD$length>0,]; # drop all dead individuals
  # adjust from index
  indD$length = indD$length * CLASS_WIDTH + MIN_SIZE - CLASS_WIDTH/2;
  
  thePlot = ggplot(indD, aes(x=t, y=length)) + theme_bw() + # geom_line(aes(group=ind)) +
            geom_smooth(size=2, color="black") + labs(title="Growth", x="Time", y="Length") +
	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)));
  if (showHist)	myPPlot(thePlot, plotDim=c(3,1))
  else myPPlot(thePlot, plotDim=c(2,1));
  thePlot = ggplot(survD, aes(x=t, y=surv)) + theme_bw() +
            geom_line(size=2, color="black") + labs(title="Survival", x="Time", y="Prop. alive") +
	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)));

  myPPlot(thePlot, plotPos=c(2,1), newWindow=FALSE);
  if (showHist) {
    thePlot = ggplot(indD[indD$t==maxT,], aes(x=length)) + theme_bw() +
              geom_histogram() + labs(x="Length", y="Count");
    myPPlot(thePlot, plotPos=c(3,1), newWindow=FALSE);
  }
  myP("Total survivors:", survV[maxT]*numInd);
  myP("Total mean per-capita eggs output in final year:", myF(totalRec/numInd));  
  
           
}

plotRecruitmentOptions <- function(lnVal=REC_LN_SDLOG, gamVal=REC_GAMMA_SHAPE, gbVal=REC_GB_GOOD) {
  numVal = 100000;
  
  rand1V = getRandomRecruitMatrix(numVal, 1, "log-normal", lnVal);
  rand2V = getRandomRecruitMatrix(numVal, 1, "gamma", gamVal);
  rand3V = getRandomRecruitMatrix(numVal, 1, "good-bad", gbVal);
  
  # drop outliers
  rand1V = dropOutliers(rand1V, c(0.025, 0.975));
  rand2V = dropOutliers(rand2V, c(0.025, 0.975));
#  rand3V = dropOutliers(rand3V, c(0.025, 0.975));
  
  myP("Ratio of largest to smallest:");
  myP("Log-normal:", max(rand1V)/min(rand1V));
  myP("Gamma:", max(rand2V)/min(rand2V));
  myP("Good-bad:", max(rand3V)/min(rand3V));
  
  theD    = data.frame(x=rand1V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Log-normal") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)));

  myPPlot(thePlot, plotDim=c(3,1), newWindow=TRUE);
  theD    = data.frame(x=rand2V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Gamma") 
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)));

  myPPlot(thePlot, plotPos=c(2,1), newWindow=FALSE);
  theD    = data.frame(x=rand3V);
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram() + labs(title="Good-bad") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)));
  myPPlot(thePlot, plotPos=c(3,1), newWindow=FALSE);
}

plotStochasticElements <- function(numVal=100000, lnVal=REC_LN_SDLOG, gamVal=DISP_GAMMA_SHAPE) {

  rand1V = getRandomRecruitMatrix(numVal, 1, "log-normal", lnVal)
  rand2V = getRandomDispSDMatrix(numVal, 1, "gamma", gamVal)
  
  # drop outliers
  rand1V = dropOutliers(rand1V, c(0.025, 0.975))
  rand2V = dropOutliers(rand2V, c(0.025, 0.975))
  
  myP("Ratio of largest to smallest:")
  myP("Log-normal:", max(rand1V)/min(rand1V))
  myP("Gamma:", max(rand2V)/min(rand2V))
  
  myMode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  myP("Mean/mode recruit:", mean(rand1V), myMode(rand1V))
  myP("Mean/mode disp:", mean(rand2V), myMode(rand2V))
    
  theD    = data.frame(x=rand1V)
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram(binwidth=(max(rand1V)-min(rand1V))/100) + theme_bw() +
            labs(x="Recruitment strength", y="Count") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))

  myPPlot(thePlot)
  theD    = data.frame(x=rand2V)
  thePlot = ggplot(theD, aes(x=x)) + geom_histogram(binwidth=(max(rand2V)-min(rand2V))/100)  + theme_bw() +
            labs(x="Mean dispersal distance (m)", y="Count") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))

  myPPlot(thePlot)
}

plotConstVectors <- function() {
  theD    = data.frame(x=gMeanSizeV, eggs=gEggsPerIndV, pMat=gPropMatureV, rec=getRecruitDistributionVector())
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=eggs)) +
            labs(x="Size", y="Eggs") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))

  myPPlot(thePlot, plotDim=c(2,2))
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=pMat)) +
            labs(x="Size", y="Prop. Maturity") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(1,2))
  
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=rec)) +
            labs(x="Size", y="Recruits") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(2,1))

  densityV = seq(0,1,by=0.1)
  theD = data.frame(x=densityV, allee=getAlleeEffect(densityV*10000))
  thePlot = ggplot(theD, aes(x=x)) + geom_line(aes(y=allee)) +
            labs(x="Density", y="Allee effect") +
  	        theme(axis.text.x=element_text(size=rel(1.25)), 
                  axis.text.y=element_text(size=rel(1.25)),
 	              axis.title.x=element_text(size=rel(1.25)),
 			      axis.title.y=element_text(size=rel(1.25)))
  myPPlot(thePlot, newWindow=FALSE, plotPos=c(2,2))
}