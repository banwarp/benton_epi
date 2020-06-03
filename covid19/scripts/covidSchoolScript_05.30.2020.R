# covidSchoolScript_05.30.2020.R

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
  plotRandomTrajs = 5                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
)
