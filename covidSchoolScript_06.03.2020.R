# covidSchoolScript_06.03.2020.R

# changes from covidSchoolScript_05.30.2020.R
# Added half-day scenarios

# not commented very much since I sort of threw this together. May comment more thoroughly later.
# Basics: 3 levels of mixing *times* three levels of R0 to generate 9 scenarios.

library(data.table)
library(ggplot2)
library(reshape2)
library(SimInf)
library(dplyr)
library(Matrix)
library(zoo)
library(viridis)

setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("covid19_SimInf_06.02.2020.R")
source("eventFunctions5.16.2020.R")

schoolPop <- 450
maxTSpan <- 105
startDate <- "2020-09-01"

# low, moderate, high R0 parameters
R0ILow <- 1.2
R0IMod <- 1.8
R0IHigh <- 2.4

R0ULow <- .3
R0UMod <- .4
R0UHigh <- .5

# low mixing parameters - classrooms are mostly isolated
schoolGradesLow <- rep(0:5,each=3) # kindergarten is 0
numClassroomsLow <- length(schoolGradesLow)
igtRateLow = 1/28               # Reciprocal of expected waiting time for in-transfer event
igtNodeNumLow = 1             # Average number/proportion of nodes that transfer at each transfer event
igtMinPropLow = .1             # Minimum proportion of node population that transfers in each event
igtMaxPropLow = .2             # Maximum proportion of node population that transfers in each event
ogtRateLow =  1/105            # Reciprocal of expected waiting time for out-transfer event
ogtNodeNumLow = 1/numClassroomsLow           # Average number/proportion of nodes that transfer at each transfer event
ogtMinPropLow = .05           # Minimum proportion of node population that transfers in each event
ogtMaxPropLow = .1           # Maximum proportion of node population that transfers in each event

# moderate mixing parameters - assumes all classrooms within a grade mix evenly, with some mixing between grades
numGradesMod <- 6            # number of grades
igtRateMod = 1/5               # Reciprocal of expected waiting time for in-transfer event
igtNodeNumMod = 1             # Average number/proportion of nodes that transfer at each transfer event
igtMinPropMod = .1             # Minimum proportion of node population that transfers in each event
igtMaxPropMod = .2             # Maximum proportion of node population that transfers in each event

# high mixing parameters - assumes older kids and younger kids mix evenly with frequent mixing between older and younger
schoolAgeGroupsHigh <- 2
igtRateHigh = 1/2               # Reciprocal of expected waiting time for in-transfer event
igtNodeNumHigh = 1             # Average number/proportion of nodes that transfer at each transfer event
igtMinPropHigh = .2             # Minimum proportion of node population that transfers in each event
igtMaxPropHigh = .3             # Maximum proportion of node population that transfers in each event

lowMixingLowR0 <- covidWrapper(
                               # parameters for lowMixingLowR0
                               folderPath = NULL, 
                               simID=paste0("covidSchool",Sys.Date(),"lowMixingLowR0"), 
                               simDate = startDate, 
                               maxT = maxTSpan, 
                               trialPop = schoolPop, 
                               I0Pop = 0, 
                               N = numClassroomsLow, 
                               nodeGroupList = schoolGradesLow, 
                               numTrials = 500, 
                               unifPop = TRUE, 
                               R0I = R0ILow, 
                               R0U = R0ULow, 
                               R0P=0, 
                               parachuteRate = 1/35,
                               
                               inGroupTransferRate = igtRateLow,               # Reciprocal of expected waiting time for in-transfer event
                               inGroupTransferNodeNum = igtNodeNumLow,             # Average number/proportion of nodes that transfer at each transfer event
                               inGroupTransferMinProp = igtMinPropLow,             # Minimum proportion of node population that transfers in each event
                               inGroupTransferMaxProp = igtMaxPropLow,             # Maximum proportion of node population that transfers in each event
                               outGroupTransferRate = ogtRateLow,             # Reciprocal of expected waiting time for out-transfer event
                               outGroupTransferNodeNum = ogtNodeNumLow,           # Average number/proportion of nodes that transfer at each transfer event
                               outGroupTransferMinProp = ogtMinPropLow,           # Minimum proportion of node population that transfers in each event
                               outGroupTransferMaxProp = ogtMaxPropLow,           # Maximum proportion of node population that transfers in each event
                               titleString = "Elementary School", 
                               
                               # (included just for safety)
                               RPhysicalDistancing = R0ILow+R0ULow, 
                               phi0 = R0ILow+R0ULow, 
                               RTarget1 = R0ILow+R0ULow, 
                               RTarget2 = R0ILow+R0ULow,   
                               
                               # other parameters to remove policy response from model
                               switchOffPolicies = 1, 
                               switchOffDay = startDate, 
                               kbDay3 = startDate, 
                               R0Spread = 0, 
                               maxPrev1 = .99, 
                               maxPrev2 = .99, 
                               pdDecay = -1, 
                               upDelay = 1,                             # Number of days after prevalence passes threshold until minor/major intervention
                               downDelay = 1,                           # Number of days after prevalence drops below threshold until intervention lifted
                               phiMoveUp = 1,                          # Rate at which phi increases when interventions are imposed
                               phiMoveDown = 1,                        # Rate at which phi decreases when interventions are lifted
                               
                               # other parameters
                               monitoringSuccess = .5,
                               isoRate = .5,  # Assume pretty good contact tracing
                               postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
                               hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
                               hospRatePost = 0,
                               hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
                               allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
                               percentile = .5,                       # percentile for graphing "worst case
                               plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

lowMixingModR0 <- covidWrapper(
  # parameters for lowMixingLowR0
  folderPath = NULL, 
  simID=paste0("covidSchool",Sys.Date(),"lowMixingModR0"), 
  simDate = startDate, 
  maxT = maxTSpan, 
  trialPop = schoolPop, 
  I0Pop = 0, 
  N = numClassroomsLow, 
  nodeGroupList = schoolGradesLow, 
  numTrials = 500, 
  unifPop = TRUE, 
  R0I = R0IMod,
  R0U = R0UMod, 
  R0P=0,
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateLow,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumLow,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropLow,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropLow,             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = ogtRateLow,             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = ogtNodeNumLow,           # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = ogtMinPropLow,           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = ogtMaxPropLow,           # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School", 
  
  # (included just for safety)
  RPhysicalDistancing = R0IMod+R0UMod, 
  phi0 = R0IMod+R0UMod, 
  RTarget1 = R0IMod+R0UMod, 
  RTarget2 = R0IMod+R0UMod,   
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1, 
  switchOffDay = startDate, 
  kbDay3 = startDate, 
  R0Spread = 0, 
  maxPrev1 = .99, 
  maxPrev2 = .99, 
  pdDecay = -1, 
  upDelay = 1,                             # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                           # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                          # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                        # Rate at which phi decreases when interventions are lifted
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

lowMixingHighR0 <- covidWrapper(
  # parameters for lowMixingLowR0
  folderPath = NULL, 
  simID=paste0("covidSchool",Sys.Date(),"lowMixingHighR0"), 
  simDate = startDate, 
  maxT = maxTSpan, 
  trialPop = schoolPop, 
  I0Pop = 0, 
  N = numClassroomsLow, 
  nodeGroupList = schoolGradesLow, 
  numTrials = 500, 
  unifPop = TRUE, 
  R0I = R0IHigh, 
  R0U = R0UHigh, 
  R0P=0, 
  parachuteRate = 1/35, 
  
  inGroupTransferRate = igtRateLow,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumLow,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropLow,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropLow,             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = ogtRateLow,             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = ogtNodeNumLow,           # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = ogtMinPropLow,           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = ogtMaxPropLow,           # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School", 
  
  # (included just for safety)
  RPhysicalDistancing = R0IHigh+R0UHigh, 
  phi0 = R0IHigh+R0UHigh, 
  RTarget1 = R0IHigh+R0UHigh, 
  RTarget2 = R0IHigh+R0UHigh,   
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1, 
  switchOffDay = startDate, 
  kbDay3 = startDate, 
  R0Spread = 0, 
  maxPrev1 = .99, 
  maxPrev2 = .99, 
  pdDecay = -1, 
  upDelay = 1,                             # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                           # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                          # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                        # Rate at which phi decreases when interventions are lifted
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

modMixingLowR0 <- covidWrapper(
  # parameters for modMixingLowR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"modMixingLowR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = numGradesMod,
  unifPop = TRUE,
  numTrials = 500,
  R0I = R0ILow,
  R0U = R0ULow,
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateMod,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumMod,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropMod,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropMod,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0ILow+R0ULow,
  phi0 = R0ILow+R0ULow,
  RTarget1 = R0ILow+R0ULow,
  RTarget2 = R0ILow+R0ULow,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

modMixingModR0 <- covidWrapper(
  # parameters for modMixingLowR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"modMixingModR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = numGradesMod,
  unifPop = TRUE,
  numTrials = 500,
  R0I = R0IMod,
  R0U = R0UMod,
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateMod,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumMod,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropMod,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropMod,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0IMod+R0UMod,
  phi0 = R0IMod+R0UMod,
  RTarget1 = R0IMod+R0UMod,
  RTarget2 = R0IMod+R0UMod,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

modMixingHighR0 <- covidWrapper(
  # parameters for modMixingLowR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"modMixingHighR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = numGradesMod,
  unifPop = TRUE,
  numTrials = 500,
  R0I = R0IHigh,
  R0U = R0UHigh, 
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateMod,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumMod,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropMod,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropMod,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0IHigh+R0UHigh,
  phi0 = R0IHigh+R0UHigh,
  RTarget1 = R0IHigh+R0UHigh,
  RTarget2 = R0IHigh+R0UHigh,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

highMixingLowR0 <- covidWrapper(
  # parameters for highMixingLowR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"highMixingLowR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = schoolAgeGroupsHigh,
  numTrials = 500,
  R0I = R0ILow,
  R0U = R0ULow,
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateHigh,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumHigh,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropHigh,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropHigh,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0ILow+R0ULow,
  phi0 = R0ILow+R0ULow,
  RTarget1 = R0ILow+R0ULow,
  RTarget2 = R0ILow+R0ULow,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

highMixingModR0 <- covidWrapper(
  # parameters for highMixingModR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"highMixingModR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = schoolAgeGroupsHigh,
  numTrials = 500,
  R0I = R0IMod,
  R0U = R0UMod,
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateHigh,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumHigh,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropHigh,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropHigh,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0IMod+R0UMod,
  phi0 = R0IMod+R0UMod,
  RTarget1 = R0IMod+R0UMod,
  RTarget2 = R0IMod+R0UMod,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)

highMixingHighR0 <- covidWrapper(
  # parameters for highMixingLowR0
  folderPath = NULL,
  simID=paste0("covidSchool",Sys.Date(),"highMixingHighR0"),
  simDate = startDate,
  maxT = maxTSpan,
  trialPop = schoolPop,
  I0Pop = 0,
  N = schoolAgeGroupsHigh,
  numTrials = 500,
  R0I = R0IHigh,
  R0U = R0UHigh,
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateHigh,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumHigh,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropHigh,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropHigh,             # Maximum proportion of node population that transfers in each event
  titleString = "Elementary School",
  
  # (included just for safety)
  RPhysicalDistancing = R0IHigh+R0UHigh,
  phi0 = R0IHigh+R0UHigh,
  RTarget1 = R0IHigh+R0UHigh,
  RTarget2 = R0IHigh+R0UHigh,  
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1,
  switchOffDay = startDate,
  kbDay3 = startDate,
  R0Spread = 0,
  maxPrev1 = .99,
  maxPrev2 = .99,
  pdDecay = -1,
  upDelay = 1,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                       # Rate at which phi decreases when interventions are lifted
  
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                       # percentile for graphing "worst case
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)


#######################
#######################
#######################

# half day parameters
# half day transferEventsDF

halfDaySchoolPop <- 540 # some duplication to round out classroom size to 15
numTrials <- 500
tspan <- 1:105
morningNodeGroupList <- rep(0:5,each=3)
afternoonNodeGroupList <- rep(6:11,each=3)
N <- length(morningNodeGroupList)+length(afternoonNodeGroupList)
NnumTrials <- N*numTrials
numGroups <- length(unique(morningNodeGroupList))+length(unique(afternoonNodeGroupList))

compartments <- c("S","E","I","uI","R","Im","pI","H","cumI","M")

E <- cbind(diag(length(compartments)),c(rep(1,length(compartments)-4),rep(0,4)))
dimnames(E) <- list(compartments,c(1:ncol(E)))

nodeTrialMat <- data.table(node=c(1:(NnumTrials)),
                           nodeGroup = rep(c(morningNodeGroupList,afternoonNodeGroupList),numTrials),
                           trial=rep(1:numTrials,each=N))

nodeAndGroupList <- split(nodeTrialMat,nodeTrialMat$nodeGroup)
nodeAndGroupList <- lapply(nodeAndGroupList,function(x) split(x$node,x$trial))

S0 <- round(rep(halfDaySchoolPop/N,NnumTrials),0)

inGroupTransferRate = rep(igtRateLow,numGroups)
inGroupTransferNodeNum = rep(igtNodeNumLow,numGroups)
inGroupTransferMinProp = rep(igtMinPropLow,numGroups)
inGroupTransferMaxProp = rep(igtMaxPropLow,numGroups)
outGroupTransferRate = rep(ogtRateLow,numGroups)
outGroupTransferNodeNum = rep(ogtNodeNumLow,numGroups)
outGroupTransferMinProp = rep(ogtMinPropLow,numGroups)
outGroupTransferMaxProp = rep(ogtMaxPropLow,numGroups)

##### Morning section #####

morningNumInTransferList <- lapply(1:numGroups,function(x) 
  lapply(1:numTrials, function(y) rpois(1,inGroupTransferRate[x]*inGroupTransferNodeNum[x]*max(tspan))))

###### IN-GROUP TRANSFER EVENTS ######
morningInTransferEventsList <- lapply(1:numGroups, transferFunction,
                                      numTList = morningNumInTransferList,
                                      tS = tspan,
                                      nT = numTrials,
                                      nList = nodeAndGroupList,
                                      minProp = inGroupTransferMinProp,
                                      maxProp = inGroupTransferMaxProp,
                                      Emat = E,
                                      nTM = nodeTrialMat,
                                      u0S = S0,
                                      outLogic = FALSE
)
morningInTransferEvents <- bind_rows(morningInTransferEventsList[lengths(morningInTransferEventsList) != 0])

###### OUT-GROUP TRANSFER EVENTS ######
morningNumOutTransferList <- lapply(1:numGroups,function(x) 
  lapply(1:numTrials, function(y) rpois(1,outGroupTransferRate[x]*outGroupTransferNodeNum[x]*max(tspan))))

morningOutTransferEventsList <- lapply(1:numGroups, transferFunction,
                                       numTList = morningNumOutTransferList,
                                       tS = tspan,
                                       nT = numTrials,
                                       nList = nodeAndGroupList,
                                       minProp = outGroupTransferMinProp,
                                       maxProp = outGroupTransferMaxProp,
                                       Emat = E,
                                       nTM = nodeTrialMat,
                                       u0S = S0,
                                       outLogic = TRUE
)
morningOutTransferEvents <- bind_rows(morningOutTransferEventsList[lengths(morningOutTransferEventsList) != 0])
morningTransferEvents <- rbind(morningInTransferEvents,morningOutTransferEvents)

rm(morningOutTransferEventsList)
rm(morningOutTransferEvents)
rm(morningInTransferEventsList)
rm(morningInTransferEvents)

##### Afternoon section #####

afternoonNumInTransferList <- lapply(1:numGroups,function(x) 
  lapply(1:numTrials, function(y) rpois(1,inGroupTransferRate[x]*inGroupTransferNodeNum[x]*max(tspan))))

###### IN-GROUP TRANSFER EVENTS ######
afternoonInTransferEventsList <- lapply(1:numGroups, transferFunction,
                                        numTList = afternoonNumInTransferList,
                                        tS = tspan,
                                        nT = numTrials,
                                        nList = nodeAndGroupList,
                                        minProp = inGroupTransferMinProp,
                                        maxProp = inGroupTransferMaxProp,
                                        Emat = E,
                                        nTM = nodeTrialMat,
                                        u0S = S0,
                                        outLogic = FALSE
)
afternoonInTransferEvents <- bind_rows(afternoonInTransferEventsList[lengths(afternoonInTransferEventsList) != 0])

###### OUT-GROUP TRANSFER EVENTS ######
afternoonNumOutTransferList <- lapply(1:numGroups,function(x) 
  lapply(1:numTrials, function(y) rpois(1,outGroupTransferRate[x]*outGroupTransferNodeNum[x]*max(tspan))))

afternoonOutTransferEventsList <- lapply(1:numGroups, transferFunction,
                                         numTList = afternoonNumOutTransferList,
                                         tS = tspan,
                                         nT = numTrials,
                                         nList = nodeAndGroupList,
                                         minProp = outGroupTransferMinProp,
                                         maxProp = outGroupTransferMaxProp,
                                         Emat = E,
                                         nTM = nodeTrialMat,
                                         u0S = S0,
                                         outLogic = TRUE
)
afternoonOutTransferEvents <- bind_rows(afternoonOutTransferEventsList[lengths(afternoonOutTransferEventsList) != 0])
afternoonTransferEvents <- rbind(afternoonInTransferEvents,afternoonOutTransferEvents)

rm(afternoonOutTransferEventsList)
rm(afternoonOutTransferEvents)
rm(afternoonInTransferEventsList)
rm(afternoonInTransferEvents)

transferEvents <- rbind(morningTransferEvents,afternoonTransferEvents)
rm(morningTransferEvents)
rm(afternoonTransferEvents)

transferEvents <- transferEvents[sample(nrow(transferEvents),ceiling(nrow(transferEvents)/4)),]

##########################
##########################
##########################

fulldaySchool <- covidWrapper(
  # parameters for lowMixingLowR0
  folderPath = NULL, 
  simID=paste0("covidSchool",Sys.Date(),"fulldaySchool"), 
  simDate = startDate, 
  maxT = maxTSpan, 
  trialPop = halfDaySchoolPop, 
  I0Pop = 0, 
  N = numClassroomsLow, 
  nodeGroupList = schoolGradesLow, 
  numTrials = 500, 
  unifPop = TRUE, 
  R0I = R0IMod, 
  R0U = R0UMod, 
  R0P=0, 
  parachuteRate = 1/35,
  
  inGroupTransferRate = igtRateLow,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumLow,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropLow,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropLow,             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = ogtRateLow,             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = ogtNodeNumLow,           # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = ogtMinPropLow,           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = ogtMaxPropLow,           # Maximum proportion of node population that transfers in each event
  titleString = "Full Day Elementary School", 
  
  # (included just for safety)
  RPhysicalDistancing = R0IMod+R0UMod, 
  phi0 = R0IMod+R0UMod, 
  RTarget1 = R0IMod+R0UMod, 
  RTarget2 = R0IMod+R0UMod,   
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1, 
  switchOffDay = startDate, 
  kbDay3 = startDate, 
  R0Spread = 0, 
  maxPrev1 = .99, 
  maxPrev2 = .99, 
  pdDecay = -1, 
  upDelay = 1,                             # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                           # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                          # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                        # Rate at which phi decreases when interventions are lifted
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                        # Plot 90th percentile of point-in-time simulations
  plotRandomTrajs = 5,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  lString = "Median"
)


halfdaySchool <- covidWrapper(
  # parameters for lowMixingLowR0
  folderPath = NULL, 
  simID=paste0("covidSchool",Sys.Date(),"halfDaySchool"), 
  simDate = startDate, 
  maxT = maxTSpan, 
  trialPop = halfDaySchoolPop, 
  I0Pop = 0, 
  N = numClassroomsLow*2, 
  nodeGroupList = schoolGradesLow, 
  numTrials = 500, 
  unifPop = TRUE, 
  R0I = (1+R0IMod)/2, 
  R0U = (1+R0UMod)/2, 
  R0P=0, 
  parachuteRate = 1/35,
  
  ##### THIS RIGHT HERE #####
  transferEventsDF = transferEvents,
  
  inGroupTransferRate = igtRateLow,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = igtNodeNumLow,             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = igtMinPropLow,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = igtMaxPropLow,             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = ogtRateLow,             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = ogtNodeNumLow,           # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = ogtMinPropLow,           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = ogtMaxPropLow,           # Maximum proportion of node population that transfers in each event
  
  titleString = "Half day Elementary School", 
  
  # (included just for safety)
  RPhysicalDistancing = (1+R0IMod+R0UMod)/2, 
  phi0 = (1+R0IMod+R0UMod)/2, 
  RTarget1 = (1+R0IMod+R0UMod)/2, 
  RTarget2 = (1+R0IMod+R0UMod)/2,   
  
  # other parameters to remove policy response from model
  switchOffPolicies = 1, 
  switchOffDay = startDate, 
  kbDay3 = startDate, 
  R0Spread = 0, 
  maxPrev1 = .99, 
  maxPrev2 = .99, 
  pdDecay = -1, 
  upDelay = 1,                             # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 1,                           # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = 1,                          # Rate at which phi increases when interventions are imposed
  phiMoveDown = 1,                        # Rate at which phi decreases when interventions are lifted
  
  # other parameters
  monitoringSuccess = .5,
  isoRate = .5,  # Assume pretty good contact tracing
  postInfectiousPeriod = 1/10,                        # Reciprocal of length of isolation period
  hospRateExp = 0,                         # Proportion of infectious/isolated that are hospitalized
  hospRatePost = 0,
  hospDeathRate = 0,                     # Hospitalized fatality rate NOT CASE FATALITY RATE
  allTraj = TRUE,                         # Logical if all simulation trajectories are plotted or just median and spread
  percentile = .5,                        # Plot 90th percentile of point-in-time simulations
  plotRandomTrajs = 5,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  lString = "Median"
)
