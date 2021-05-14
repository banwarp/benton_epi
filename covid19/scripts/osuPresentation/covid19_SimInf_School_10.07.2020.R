# covid19_SimInf_School_10.07.2020.R

# copyright Peter Banwarth Benton County Health Department 2020

# Using package SimInf to model stochastic disease spread
# https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
# https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
# https://github.com/stewid/SimInf

# Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
# Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

# Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
# International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

# changes from covid19_SimInf_School_09.16.2020.R
# Added option to test cohorts when a case is identified; tracks cumulative isolations: cumIso, whenever cumIso increases, tests whole cohort
# Added subtitle to plots

# changes from covid19_SimInf_School_08.08.2020.R
# changed confidence interval to 80% two sided
# added option for user-defined axes for graphs
# added plot for days in school/quarantine/etc.
# added tabulation of number of days in school etc.

# changes from covid19_SimInf_School_08.06.2020.R
# changed initial infections to total per trial; added probability for whether infections exist
# Increased schooltime R0 to be 2 - noSchoolFactor

# changes from covid19_SimInf_08.04.2020.R
# Changed from hourly to two periods: In school and out of school. In school is 1/2 of day; out of school is 1/2 of day.

# changes from covid19_SimInf_School_07.04.2020.R
# Brought back parachuting

# Updated parameters to most recent literature
# Sources:
# R0; asymptomatic transmission; presymptomatic transmission: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box
# Base R0 = 2.5
# Need to review the following:
# Asymptomatic/presymptomatic transmission fraction = .75/(1 + .75 + .1) = .4; denominator is relative infectiousness of symptomatic + presymptomatic + postsymptomatic
# Postsymptomatic tranmission fraction = .1/(1+.75+.1) = .05 (just a guess)
# Symptomatic transmission fraction = 1/1.85 = .55

# Child/Adult transmission differential: https://www.medrxiv.org/content/10.1101/2020.05.20.20108126v1
# Child is 44% of adult 
# Need more research on this in particular

# Probability of detection for each infectious state; Just a guess at this point
# pre/post/asymp = 0.1
# symp = .8

# Probability of symptomatic infection versus asymptomatic for middle/high school: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.29.2001352#html_fulltext
# 43% for middle/high school
# Probability of symptoms for staff: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box
# 65%
# Probability of symptomatic infection versus asymptomatic for elementary: without other data: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.29.2001352#html_fulltext
# 43% for elementary

# changes from covid19_SimInf_School06.20.2020.R
# Changed pts_fun to let user define the case number for which quarantine is started.
# Removed false positives - too complex for right now

# changes from covid19_SimInf_06.16.2020.R
# Major changes from the base script
# Model schools instead of mixed populations
# Added preSymptomatic, Symptomatic, post-symptomatic, and asymptomatic compartments.
# stationary and transitory compartments
# Ability to isolate individuals (maybe), classrooms (yes), grades (not yet - based on node groups), and schools (not yet - based on trials)
# Non zero, non-unit probability of detecting pre, symp, post, and asymp infections
# Renamed covidWrapper to covidSimulator
# See covid19_SimInf_06.16.2020.R for history of other changes

#########################################

covidSimulator <- function(
  ### Simulation parameters
  folderPath = NULL,                       # folder path for subroutines and output
  simID = "covid.mo.day.2020.School",      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                          # Number of trials in simulation
  
  ### population parameters
  N = 18,                                  # Number of classrooms
  nodeGroupList = NULL,                    # List of group IDs for node groups. Default is 1 group of nodes.
  Symp0nodeGroups = NULL,                  # Node groups where initial infectious are distributed
  Trans0nodes = NULL,                      # Nodes group where transitory individuals are distributed
  S0Pops = 450,                            # Initial population - stationary; total per trial. Infectious will be extracted from this total.
  preSymp0Pops = 0,                        # Expected initial number of presymptomatic - stationary; total per trial
  Symp0Pops = 0,                           # Expected initial number of symptomatic - stationary; total per trial
  postSymp0Pops = 0,                       # Expected initial number of post-symptomatic - stationary; total per trial
  aSymp0Pops = 0,                          # Expected initial number of asymptomatic - stationary; total per trial
  R0Pops = 0,                              # Initial number of recovered - stationary; count per node
  Im0Pops = 0,                             # Initial number of immune - stationary; count per node
  Iso0Pops = 0,                            # Initial number of isolated - stationary; count per nodey
  S0Popt = 10,                             # Initial population - transitory; total per trial. Infectious will be extracted from this total. 
  preSymp0Popt = 0,                        # Expected initial number of pre-symptomatic - transitory; total per trial
  Symp0Popt = 0,                           # Expected initial number of symptomatic - transitory; total per trial
  postSymp0Popt = 0,                       # Expected initial number of post-symptomatic - transitory; total per trial
  aSymp0Popt = 0,                          # Expected initial number of asymptomatic - transitory; total per trial
  R0Popt = 0,                              # Initial number of recovered - transitory; count per node
  Im0Popt = 0,                             # Initial number of immune - transitory; count per node
  Iso0Popt = 0,                            # Initial number of isolated - transitory; count per node
  
  ### gdata parameters
  R0Base = 3,                              # Baseline R0; compromize betwen 2.5 from CDC and 3.28 from OHA
  R0preFrac = .4,                          # Presymptomatic R0 fraction (what proportion of new cases come from presymptomatic)
  R0sympFrac = .55,                        # Symptomatic R0 fraction
  R0postFrac = .05,                        # Postsymptomatic R0 fraction
  R0asympFrac = .4,                        # Asymptomatic R0 fraction
  R0ComplianceFrac = 1,                    # Proportionate reduction in R0 with high compliance to masks/distancing
  R0VentilationFrac = 1,                   # Proportionate reduction in R0 with good ventilation
  studentTeacherDiff = c(.33,1,.5,1),      # Differential R0 for stationary-stationary, transitory-transitory, and stationary-transitory interactions
  preSympPeriod = .5,                      # Reciprocal of presymptomatic period
  symptomaticPeriod = .25,                 # Reciprocal of symptomatic period
  postSympPeriod = .1667,                  # Reciprocal of post-symptomatic period
  aSympPeriod = .1667,                     # Reciprocal of asymptomatic period
  isoPeriod = .1,                          # Reciprocal of isolation period
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = .01,                     # Reciprocal of temporary immunity period, after which R becomes Im or S
  
  ### ldata parameters
  R0Spread = 0,                            # Uniform variation in R0 across trials (measured as %change from R0)
  
  ### continuous initialized parameters (v0)
  symptomaticProp = c(.43,.65),            # Average proportion of infectious who are symptomatic; (stationary, transitory)
  symptomaticVariance = 0,                 # Uniform variance of symptomatic proportion
  schoolDay = 1,                           # Indicator for school day: 1 = school day, 0 = evening, list for each node or single entry
  schoolWeek = 1,                          # Indicator for school week: 1 = week, 0 = weekend, list for each node or single entry
  schoolTerm = 1,                          # Indicator for school term: 1 = school term, 0 = break, list for each node or single entry
  noQuarantine = 1,                        # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
  weekDays = list(c(0,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### pts_fun parameters
  cosAmp = 0.125,                          # Amplitude of seasonal variation in beta
  quarantineDaysClassroom = 14,            # Length in days of quarantine
  quarantineDaysGrade = 14,                # Length in days of quarantine
  quarantineDaysSchool = 28,               # Length in days of quarantine
  preDetectSuccess = 0.1,                  # Probability of detecting a pre-symptomatic individual
  sympDetectSuccess = 0.8,                 # Probability of detecting a symptomatic individual
  postDetectSuccess = 0.1,                 # Probability of detecting a post-symptomatic individual
  asympDetectSuccess = 0.1,                # Probability of detecting an asymptomatic individual
  detectionProbabilityVar = 0,             # Proportionate variable decrease in detection probability - only decreases probability, never increases
  classroomThreshold = 2,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = NULL,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 16,                    # Number of active detected cases before school quarantined
  countClassrooms = TRUE,                  # Whether to count classrooms on quarantine when counting active detected cases
  nightFactor = .05,                       # Factor to reduce beta when nighttime
  noSchoolFactor = .2,                     # Factor to reduce beta when outside of school
  quarantineFactor = .1,                   # Factor to reduce beta when quarantine
  cohortTesting = 1,                       # Binary logical: test cohort if a case is identified; 1 = TRUE, 0 = FALSE
  testSensitivity = .75,                   # Test sensitivity (true positives); for simplicity assumes no false positives
  
  ### *day* specified infection events parameters - list to allow multiple events, use list structure to accommodate stationary/transitory 
  infectionEventsDF = NULL,                # data.frame with specified infection events
  infectionCount = c(),                    # Number of infections caused by the infection spreader
  infectionNodes = c(),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c(),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(),                       # Infection group type: 1=stationary or 2=transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraEventsDF = NULL,                     # pre-built parachuter events data frame if the events will be preset before the script is run
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  # parachuteRate = c(1/10,1/25),                  # Reciprocal of expected waiting time in days for a parachute event
  parachuteRate = c(0,0),                  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  transferEventsDF = NULL,                 # pre-built *hourly* transfer events data frame if the events will be preset before the script is run
  inGroupTransferRate = rep(0,4),          # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = rep(0,4),       # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = rep(0,4),       # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = rep(0,4),       # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),    # Maximum proportion of node population that transfers in each event
  
  ### plot parameters
  rollM = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  allTraj = FALSE,                          # Logical if all simulation trajectories are plotted or just median and spread
  plotRandomTrajs = 5,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  percentile = .5,                         # Percentile of point-in-time simulations to plot, default to median
  confIntv = .8,                          # two-sided confidence interval for plotting spread
  plotGroups = NULL,                       # Which node groups to plot. NULL plots all
  sumGroups = TRUE,                        # Whether to sum the node groups.
  dateBreaks = "1 month",                  # Plot parameter, x-axis displays months
  dataTitle = "Infections",                # Title of data that is plotted
  miny = NULL,                             # optional minimum for y-axis; NULL uses built-in
  maxy = NULL,                             # optional maximum for y-axis; NULL uses built-in
  titleString = "Generic Title",           # Title of plot
  subtitleString = "Subtitle",             # Subtitle of plot
  xString = "Date",                        # Title of x-axis
  yString = "Frequency",                   # Title of y-axis
  lString = "Median",                      # Title of legend
  cString = NULL                           # Plot caption
  ) {
  
  ######## Begin function #######
  ######## Begin function #######
  ######## Begin function #######
  ######## Begin function #######
  ######## Begin function #######
  
  # set.seed(123)
  
  library(ggplot2)
  library(reshape2)
  library(SimInf)
  library(data.table)
  library(dplyr)
  library(simFrame)
  library(Matrix)
  library(zoo)
  library(viridis)
  library(combinat)
  
  if(!is.null(folderPath)) {setwd(folderPath)}
  
  if(is.null(folderPath)) {setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")}
  source("transitionsScriptSchool_10.07.2020.R")
  source("eventFunctionsSchool_08.06.2020.R")
  source("simInfPlottingFunction_09.16.2020.R")
  source("pts_funScript_school_10.07.2020.R")
  if(is.null(folderPath)) {setwd("../")}

  # saving parameters
  if(is.null(folderPath)) {setwd("./parameters")}
  sink(paste0("parms",simID,".txt"))
  print("All parameters:")
  print(mget(names(formals())))
  sink()
  if(is.null(folderPath)) {setwd("../")}

  ####### Setting additional parameters #######
  
  # Start of simulation
  startofSimDay <- as.numeric(as.Date(simDate)) - as.numeric(as.Date("2020-01-01"))

  # Time span in 2-period days
  tspan <- 1:(maxT*2)
  
  # additional plotting parameters
  fileName <- paste0("data",simID,".txt")
  plotName <- paste0("plot",simID)
  if(is.null(cString)) {
    cString <- plotName
  }
  
  # Trial and node parameters
  if(is.null(nodeGroupList)) {nodeGroupList <- rep(1,N)} # default is one nodeGroup
  numGroups <- length(unique(nodeGroupList))
  NnumTrials <- N*numTrials
  nodeTrialMat <- data.table(node=c(1:(NnumTrials)),
                             nodeGroup = rep(nodeGroupList,numTrials),
                             trial=rep(1:numTrials,each=N))
  NList <- split(nodeTrialMat$node,nodeTrialMat$trial)
  nodeAndGroupList <- split(nodeTrialMat,nodeTrialMat$nodeGroup)
  nodeAndGroupList <- lapply(nodeAndGroupList,function(x) split(x$node,x$trial))
  if(is.null(Symp0nodeGroups)) {Symp0nodeGroups <- unique(nodeGroupList)}
  if(is.null(Trans0nodes)) {Trans0nodes <- 1:N}
  
  # Model parameters
  # compartments
  compartments <- c("Ss",        # Susceptible stationary
                    "St",        # Susceptible transitory
                    "preSymps",  # Pre-symptomatic stationary stationary
                    "preSympt",  # Pre-symptomatic transitory
                    "Symps",     # symptomatic Infectious stationary
                    "Sympt",     # symptomatic Infectious transitory
                    "postSymps", # post-symptomatic stationary
                    "postSympt", # post-symptomatic transitory
                    "aSymps",    # asymptomatic Infectious stationary
                    "aSympt",    # asymptomatic Infectious transitory
                    "Rs",        # Recovered stationary
                    "Rt",        # Recovered transitory
                    "Ims",       # Immune stationary
                    "Imt",       # Immune transitory
                    "Isos",      # Isolated stationary
                    "Isot",      # Isolated transitory
                    "cumIs",     # cumulative stationary Infections (used for final tabulation)
                    "cumIt",     # cumulative transitory Infections (used for final tabulation)
                    "cumIso"     # cumulative isolations for cohort testing
                    )
  
  # number of compartments that don't participate in stationary transfers and in transitory transfers
  xCompt <- 5 # Isos, Isot, cumIs, cumIt, cumIso
  
  # More model parameters
  # global parameters
  gdata = data.frame(
    betaPre = R0Base*R0preFrac*preSympPeriod/2,         # Transmission rate among presymptomatic
    betaSymp = R0Base*R0sympFrac*symptomaticPeriod/2,   # Transmission rate among symptomatic
    betaPost = R0Base*R0postFrac*postSympPeriod/2,      # Transmission rate among post symptomatic
    betaA = R0Base*R0asympFrac*aSympPeriod/2,           # Transmission rate among asymptomatic
    R0ComplianceFrac = R0ComplianceFrac,                # Proportionate reduction in R0 for high/low mask/distancing compliance
    R0VentilationFrac = R0VentilationFrac,              # Proportionate reduction in R0 for high/low ventilation
    preSympPeriod = preSympPeriod/2,                    # Reciprocal of pre symptomatic period
    symptomaticPeriod = symptomaticPeriod/2,            # Reciprocal of infectious period
    postSympPeriod = postSympPeriod/2,                  # Reciprocal of post symptomatic period
    aSympPeriod = aSympPeriod/2,                        # Reciprocal of asymptomatic period
    isoPeriod = isoPeriod/2,                            # Reciprocal of isolation period
    reSuscepRate = reSuscepRate,                        # Proportion of recovered who become re-susceptible
    tempImmPeriod = tempImmPeriod/2,                    # Reciprocal of recovered period before re-susceptibility
    ssD = studentTeacherDiff[1],                        # Differential on R0 for stationary-stationary interactions
    stD = studentTeacherDiff[2],                        # Differential on R0 for stationary-transitory interactions
    tsD = studentTeacherDiff[3],                        # Differential on R0 for transitory-stationary interactiosn
    ttD = studentTeacherDiff[4]                         # Differential on R0 for transitory-transitory interactions
  )
  
  # local parameters
  # firstGroupNode: the first node in each group
  # numberGroupNode: the number of nodes in each group
  # firstTrialnode: the first node in each trial
  # betaRandomizer: random factor for beta
  ldata <- data.frame(
    firstGroupNode = unlist(lapply(which(!duplicated(nodeTrialMat[,c("nodeGroup","trial")])), function(x) rep(x-1,length(nodeGroupList[which(nodeGroupList==nodeTrialMat$nodeGroup[x])])))),
    numberGroupNode = rep(unlist(lapply(1:numGroups,function(x) rep(length(which(nodeGroupList==x)),length(which(nodeGroupList==x))))),numTrials),
    firstTrialNode = sort(rep.int(0:(numTrials-1)*N,N),method="quick"),
    betaRandomizer = rep(runif(numTrials,min=1-R0Spread,max=1+R0Spread),each=N)
  )
  
  # rep weekDays and breakDays if needed and split into different lists
  if(length(weekDays) == 0) {weekDays <- list(c(-1,-1))}
  if(length(breakDays) == 0) {breakDays <- list(c(-1,-1))}
  
  # if breakDays are dates
  if(is.character(unlist(breakDays))) {
    breakDays <- lapply(breakDays,function(x) lapply(x, function(y) as.numeric(as.Date(y)) - as.numeric(as.Date("2020-01-01"))))
  }
  
  maxWeekDays <- max(unlist(lapply(weekDays,length)))
  maxBreakDays <- max(unlist(lapply(breakDays,length)))
  
  if(length(weekDays) == 1) {weekDays <- rep(weekDays,N)}
  if(length(breakDays) == 1) {breakDays <- rep(breakDays,N)}
  
  weekDayDF <- data.frame(matrix(data=-1,nrow=N,ncol=maxWeekDays))
  for(i in 1:length(weekDays)) {
    weekDayDF[i,c(1:length(weekDays[[i]]))] <- weekDays[[i]]
  }
  names(weekDayDF) <- paste0("wd",1:maxWeekDays)
  
  breakDayDF <- data.frame(matrix(data=-1,nrow=N,ncol=maxBreakDays))
  for(i in 1:length(breakDays)) {
    breakDayDF[i,c(1:length(breakDays[[i]]))] <- breakDays[[i]]
  }
  names(breakDayDF) <- paste0("bd",1:maxBreakDays)
  
  ldata <- cbind(ldata,weekDayDF,breakDayDF)
 
  # Transitions
  transitions <- transitionsScript()
  
  # Initial state
  # Uniformly distribute stationary population between nodes
  S0s <- round(rep(S0Pops/N,NnumTrials),0)
  
  # Distribute transitory population according to Trans0nodeGroups
  S0t <- rep(0,N)
  if(S0Popt <= length(Trans0nodes)) {
    if(length(Trans0nodes > 1)) {
      S0t[sample(Trans0nodes,S0Popt)] <- 1
    } else {
      S0t[Trans0nodes] <- 1
    }
  } else {
    repNum <- floor(S0Popt/length(Trans0nodes))
    repRem <- S0Popt %% length(Trans0nodes)
    S0t[Trans0nodes] <- repNum
    S0t[sample(Trans0nodes,repRem)] <- 1+repNum
  }
  S0t <- rep(S0t,numTrials)  

  # loading infectious individuals
  
  # function to randomize placement of infectious individuals
  initialInfectiousFunction <- function(z, infPop) {
    y <- rep(0,N)
    numInf <- rpois(1,infPop)
    if(numInf > 0) {
      infT <- table(sample(N,numInf,replace=TRUE))
      y[as.numeric(names(infT))] <- infT  
    }
    return(y)
  }
  
  # selecting individuals from each node; switching them from susceptible to infectious
  # listing infectious pops for function
  inf0PopsList <- c(preSymp0Pops,Symp0Pops,postSymp0Pops,aSymp0Pops)
  inf0PoptList <- c(preSymp0Popt,Symp0Popt,postSymp0Popt,aSymp0Popt)
  # initializing list of infectious pops
  inf0s <- vector(mode = "list", length = 4)
  inf0t <- vector(mode = "list", length = 4)
  # once for each infectious group: pre, S, post, a
  for(i in 1:4) {
    inf0s[[i]] <- unlist(lapply(1:numTrials, initialInfectiousFunction, infPop = inf0PopsList[i]))
    S0s <- S0s - inf0s[[i]]
    inf0t[[i]] <- unlist(lapply(1:numTrials, initialInfectiousFunction, infPop = inf0PoptList[i]))
    S0t <- S0s - inf0t[[i]]
  }
  
  # Initial state matrix
  u0 <- data.frame(
    Ss = S0s,
    St = S0t,
    preSymps = inf0s[[1]],
    preSympt = inf0t[[1]],
    Symps = inf0s[[2]],
    Sympt = inf0t[[2]],
    postSymps = inf0s[[3]],
    postSympt = inf0t[[3]],
    aSymps = inf0s[[4]],
    aSympt = inf0t[[4]],
    Rs = rep(R0Pops, NnumTrials),
    Rt = rep(R0Popt, NnumTrials),
    Ims = rep(Im0Pops, NnumTrials),
    Imt = rep(Im0Popt, NnumTrials),
    Isos = rep(Iso0Pops, NnumTrials),
    Isot = rep(Iso0Popt, NnumTrials),
    cumIs = rep(Symp0Pops, NnumTrials),
    cumIt = rep(Symp0Popt, NnumTrials),
    cumIso = rep(0,NnumTrials)
  )
  
  # Initialized continuous variables
  v0 <- data.frame(
    outOfSchool = rep(1,NnumTrials),                           # out-of-school factor for beta
    season = rep(1,NnumTrials),                                # seasonality factor for beta
    symptomaticProps = rep(symptomaticProp[1],NnumTrials),     # proportion of infections that are symptomatic; stationary
    symptomaticPropt = rep(symptomaticProp[2],NnumTrials),     # proportion of infections that are symptomatic; transitory
    schoolDay = rep(schoolDay,NnumTrials),                     # Indicator for school day: 1 = in school, 0 = out of school
    schoolWeek = rep(schoolWeek,NnumTrials),                   # Indicator for school week: 1 = week, 0 = weekend
    schoolTerm = rep(schoolTerm,NnumTrials),                   # Indicator for school term: 1 = school term, 0 = break
    preDetectSuccess = rep(preDetectSuccess,NnumTrials),       # Probability of detecting a pre-symptomatic individual
    sympDetectSuccess = rep(sympDetectSuccess,NnumTrials),     # Probability of detecting a symptomatic individual
    postDetectSuccess = rep(postDetectSuccess,NnumTrials),     # Probability of detecting a post-symptomatic individual
    asympDetectSuccess = rep(asympDetectSuccess,NnumTrials),   # Probability of detecting an asymptomatic individual
    noQuarantineC = rep(noQuarantine,NnumTrials),              # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
    quarTimerC = rep(0,NnumTrials),                            # Timer for length of quarantine - classroom
    quarCounterC = rep(0,NnumTrials),                          # Counter for number of quarantine events - classroom
    noQuarantineG = rep(noQuarantine,NnumTrials),              # Indicator for grade level quarantine: 1 = no quarantine, 0 = quarantine
    quarTimerG = rep(0,NnumTrials),                            # Timer for length of quarantine - grade
    quarCounterG = rep(0,NnumTrials),                          # Counter for number of quarantine events - grade
    noQuarantineS = rep(noQuarantine,NnumTrials),              # Indicator for school level quarantine: 1 = no quarantine, 0 = quarantine
    quarTimerS = rep(0,NnumTrials),                            # Timer for length of quarantine - school
    quarCounterS = rep(0,NnumTrials),                          # Counter for number of quarantine events - school
    isoC = rep(0,NnumTrials),                                  # Counter for number of isolated in classroom
    isoG = rep(0,NnumTrials),                                  # Counter for number of isolated in grade
    isoS = rep(0,NnumTrials),                                  # Counter for number of isolated in school
    cumIsoCount = rep(0,NnumTrials)                            # Tracks cumulative isolations for cohort testing
  )

  # building pts_fun
  # Safeties for thresholds
  if(is.null(classroomThreshold)) {classroomThreshold <- S0Pops+S0Popt+1}
  if(is.null(gradeThreshold)) {gradeThreshold <- S0Pops+S0Popt+1}
  if(is.null(schoolThreshold)) {schoolThreshold <- S0Pops+S0Popt+1}
  if(gradeThreshold <= classroomThreshold) {classroomThreshold <- S0Pops+S0Popt+1} # Overrides classroom quarantine if grade has equal or lower threshold
  if(schoolThreshold <= min(classroomThreshold,gradeThreshold)) {
    classroomThreshold <- S0Pops+S0Popt+1 # Overrides classroom quarantine if school has equal or lower threshold
    gradeThreshold <- S0Pops+S0Popt+1 # Overrides grade quarantine if school has equal or lower threshold
  }

  # Convert countClassrooms to numeric
  if(countClassrooms) {
    countClassrooms <- 1
  } else {countClassrooms <- 0}
  
  # Equalize quarantine days if grade and school are null
  if(is.null(quarantineDaysGrade)) {quarantineDaysGrade <- quarantineDaysClassroom}
  if(is.null(quarantineDaysSchool)) {quarantineDaysSchool <- quarantineDaysClassroom}
  
  pts_fun <- pts_funScriptSchool(
    startDay = startofSimDay,               # start of simulation
    enn = N,                                # number of nodes in a trial
    numComp = length(compartments),         # number of compartments
    ncolV0 = ncol(v0),                      # number of continuous variables
    cosA = cosAmp,                          # Amplitude of seasonal variation in beta
    sympProp0s = symptomaticProp[1],        # Baseline symptomatic proportion that gets re-randomized each time step; stationary
    sympProp0t = symptomaticProp[2],        # Baseline symptomatic proportion that gets re-randomized each time step; transitory
    sympVar = symptomaticVariance,          # Uniform variance of symptomatic proportion
    preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
    sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
    postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
    asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
    detVar = detectionProbabilityVar,       # Proportionate variable decrease in probability
    cThreshold = classroomThreshold,        # Threshold for classroom quarantine
    gThreshold = gradeThreshold,            # Threshold for grade quarantine
    sThreshold = schoolThreshold,           # Threshold for school quarantine
    cClassrooms = countClassrooms,          # Logical if to include classrooms that are on quarantine when deciding if to quarantine grades/school
    qDaysC = quarantineDaysClassroom,       # Length in days of quarantine
    qDaysG = quarantineDaysGrade,           # Length in days of quarantine
    qDaysS = quarantineDaysSchool,          # Length in days of quarantine
    wDays = maxWeekDays,                    # number of school week start/stops
    bDays = maxBreakDays,                   # number of school break start/stops
    noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
    night = nightFactor,                    # proportionate reduction in beta at night
    quarFactor = quarantineFactor,          # proportionate reduction in beta due to quarantine
    cohortTest = cohortTesting,             # Binary variable for testing cohort when case is identified
    testSens <- testSensitivity             # Test sensitivity; assumes no false positives
  )
  
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  
  # create events: occasional infections enter into the nodes; occasional transfers between the nodes
  
  # Events matrix
  # leaves off Isos, Isot, cumIs, and cumIt from last column
  E <- cbind(diag(length(compartments)),
             c(rep(c(1,0),(length(compartments)-xCompt)/2),rep(0,xCompt)),
             c(rep(c(0,1),(length(compartments)-xCompt)/2),rep(0,xCompt)),
             c(rep(1,length(compartments)-xCompt),rep(0,xCompt)))
  dimnames(E) <- list(compartments,c(1:ncol(E)))

  # Shift matrix
  Nmat <- matrix(c(2,2,rep(0,length(compartments)-2)),ncol=1)
  dimnames(Nmat) <- list(compartments,1)
  
  ###### Infection EVENTS #####
  if(is.null(infectionEventsDF)) {
    
    # if infectionNodeGroups not defined, spread infections evenly across node groups
    if(is.null(infectionNodeGroups)) {infectionNodeGroups <- lapply(infectionCount,function(x) unique(nodeGroupList))}
    
    if(length(infectionCount) > 0){
      eligibleinfectionNodes <- lapply(infectionNodeGroups,function(x) nodeTrialMat[which(nodeTrialMat$nodeGroup %in% x & nodeTrialMat$trial==1),"node"]$node)
      
      infectionNodeList <- lapply(1:length(infectionNodes),
                                  function(x) if(length(eligibleinfectionNodes[[x]])>1)
                                  {sort(sample(eligibleinfectionNodes[[x]],min(length(eligibleinfectionNodes[[x]]),infectionNodes[x])))}
                                  else{eligibleinfectionNodes[[x]]})
      
      infectionEventList <- lapply(1:length(infectionCount), infectionFunction,
                                   IList=infectionCount,
                                   nodeList=infectionNodeList,
                                   dayList=infectionDate,
                                   spreadList=infectionSpread,
                                   cohortList = infectionST,
                                   startDay=startofSimDay,
                                   nT = numTrials,
                                   enn = N)
      infectionEvents <- bind_rows(infectionEventList[lengths(infectionEventList) != 0])
    }
  } else {infectionEvents <- infectionEventsDF}
  
  ###### PARACHUTE EVENTS ######
  # Follows a poisson process, shifts one or more individuals from susceptible to presymptomatic
  # Use fractions to avoid trying to shift 1+ individuals from a pop of 0
  
  if(is.null(paraEventsDF)) {
    
    # Allowable node groups: all are allowed if none are selected
    if(is.null(parachuteNodeGroups)) {parachuteNodeGroups <- list(unique(nodeGroupList),unique(nodeGroupList))}
    
    parachuteEvents <- data.frame(NULL)
    paraCompartments <- c("Ss","St")
    paraPops <- c(S0Pops/N,S0Popt)
    
    # 2 iterations: stationary parachuters; transitory parachuters
    # compensate for two-period by dividing/multiplying by 2 and subtracting 1 to set period at mod 0
    for(i in 1:2) {
      # number of events in each trial
      numParachuteList <- lapply(1:numTrials,function(x) rpois(1,parachuteRate[i]*max(tspan)/2))
      parachuteNTM <- nodeTrialMat[nodeTrialMat$nodeGroup %in% parachuteNodeGroups[[i]],]
      pNList <- split(parachuteNTM$node,parachuteNTM$trial)

      if(max(unlist(numParachuteList)) > 0){
        # parachute events data frame
        parachuteEV <- data.frame(
          event = "intTrans",
          time = unlist(lapply(numParachuteList, function(x) 2*ceiling(rbeta(x,paraMu[i],paraSig[i])*max(tspan)/2) - 1)),
          node = unlist(lapply(c(1:numTrials),function(x) sample(pNList[[x]],numParachuteList[[x]],replace=TRUE))),
          dest = 0,
          n = 0,
          proportion = min(parachuteNum[i]/paraPops[i],1),
          select = which(compartments == paraCompartments[i]),
          shift = 1
        )
        # Adding to the cumulative infections
        parachuteCum <- parachuteEV
        parachuteCum$event <- "enter"
        parachuteCum$select <- which(compartments == paraCompartments[i])+16
        parachuteCum$shift <- 0
        parachuteEV <- rbind(parachuteEV,parachuteCum)
        parachuteEvents <- rbind(parachuteEvents,parachuteEV)
      }
    }
  } else {parachuteEvents <- paraEventsDF}
  
  ###### TRANSFER EVENTS ######
  if(is.null(transferEventsDF)) {

    transferEvents <- data.frame(NULL)
    
    # 4 iterations: stationary in school; stationary out of school; transitory in school; transitory out of school
    for(i in 1:4) {
      
      iGTR <- inGroupTransferRate[[i]]
      iGTNN <- inGroupTransferNodeNum[[i]]
      iGTMinP <- inGroupTransferMinProp[[i]]
      iGTMaxP <- inGroupTransferMaxProp[[i]]
      oGTR <- outGroupTransferRate[[i]]
      oGTNN <- outGroupTransferNodeNum[[i]]
      oGTMinP <- outGroupTransferMinProp[[i]]
      oGTMaxP <- outGroupTransferMaxProp[[i]]
      cGTR <- crossGroupTransferRate[[i]]
      cGTNN <- crossGroupTransferNodeNum[[i]]
      cGTMinP <- crossGroupTransferMinProp[[i]]
      cGTMaxP <- crossGroupTransferMaxProp[[i]]
      
      # Making list of transfer parameters if they are the same across groups
      if(numGroups > 1 & length(iGTR) == 1)
      {
        iGTR = rep(iGTR,numGroups)
        iGTNN = rep(iGTNN,numGroups)
        iGTMinP = rep(iGTMinP,numGroups)
        iGTMaxP = rep(iGTMaxP,numGroups)
        oGTR = rep(oGTR,numGroups)
        oGTNN = rep(oGTNN,numGroups)
        oGTMinP = rep(oGTMinP,numGroups)
        oGTMaxP = rep(oGTMaxP,numGroups)
        cGTR = rep(cGTR,numGroups)
        cGTNN = rep(cGTNN,numGroups)
        cGTMinP = rep(cGTMinP,numGroups)
        cGTMaxP = rep(cGTMaxP,numGroups)
      }
      
      ###### IN-GROUP TRANSFER EVENTS ######
      if(iGTNN[1]<1) {iGTNN <- ceiling(iGTNN*N)} # converts proportion to number
      
      # Number of transfer events
      numInTransferList <- lapply(1:numGroups,function(x) 
        lapply(1:numTrials, function(y) rpois(1,iGTR[x]*iGTNN[x]*max(tspan))))
      
      inTransferEventsList <- lapply(1:numGroups, transferFunction,
                                     popType = i,
                                     numTList = numInTransferList,
                                     tS = tspan,
                                     nT = numTrials,
                                     nList = nodeAndGroupList,
                                     minProp = iGTMinP,
                                     maxProp = iGTMaxP,
                                     selectCol = length(compartments)+ceiling(i/2),
                                     nTM = nodeTrialMat,
                                     u0S = u0$S,
                                     outLogic = "in"
      )
      inTransferEvents <- bind_rows(inTransferEventsList[lengths(inTransferEventsList) != 0])
      
      ###### OUT-GROUP TRANSFER EVENTS ######
      if(numGroups > 1) {
        if(oGTNN[1]<1) {oGTNN<- ceiling(oGTNN*N)} # converts proportion to number
        
        # Number of transfer events
        numOutTransferList <- lapply(1:numGroups,function(x) 
          lapply(1:numTrials, function(y) rpois(1,oGTR[x]*oGTNN[x]*max(tspan))))
        
        outTransferEventsList <- lapply(1:numGroups, transferFunction,
                                        numTList = numOutTransferList,
                                        tS = tspan,
                                        nT = numTrials,
                                        nList = nodeAndGroupList,
                                        minProp = oGTMinP,
                                        maxProp = oGTMaxP,
                                        selectCol = length(compartments)+ceiling(i/2),
                                        nTM = nodeTrialMat,
                                        u0S = u0$S,
                                        outLogic = "out"
        )
        outTransferEvents <- bind_rows(outTransferEventsList[lengths(outTransferEventsList) != 0])
        transferEvents <- rbind(transferEvents,inTransferEvents,outTransferEvents)
        rm(outTransferEventsList)
        rm(outTransferEvents)
        
        ###### CROSS-GROUP TRANSFER EVENTS ######
        if(cGTNN[1]<1) {cGTNN<- ceiling(cGTNN*N)} # converts proportion to number
          
        numCrossTransferList <- lapply(1:numGroups,function(x) 
          lapply(1:numTrials, function(y) rpois(1,cGTR[x]*cGTNN[x]*max(tspan))))
          
        crossTransferEventsList <- lapply(1:numGroups, transferFunction,
                                          numTList = numCrossTransferList,
                                          tS = tspan,
                                          nT = numTrials,
                                          nList = nodeAndGroupList,
                                          minProp = cGTMinP,
                                          maxProp = cGTMaxP,
                                          selectCol = length(compartments)+ceiling(i/2),
                                          nTM = nodeTrialMat,
                                          u0S = u0$S,
                                          outLogic = "cross"
        )
        crossTransferEvents <- bind_rows(crossTransferEventsList[lengths(crossTransferEventsList) != 0])
        transferEvents <- rbind(transferEvents,crossTransferEvents)
        rm(crossTransferEventsList)
        rm(crossTransferEvents)
          
      } else {transferEvents <- rbind(transferEvents,inTransferEvents)}
      rm(inTransferEventsList)
      rm(inTransferEvents)
    } # end of i in 1:4
  } else {transferEvents <- transferEventsDF}
  
  
    
  ###### ALL EVENTS ######
  allEvents <- data.frame(NULL)
  if(length(infectionCount | !is.null(infectionEventsDF))>0) { allEvents <- rbind(allEvents,infectionEvents)}
  if(max(parachuteRate) > 0) { allEvents <- rbind(allEvents,parachuteEvents)}
  if(max(inGroupTransferRate,outGroupTransferRate)>0 | !is.null(transferEventsDF)) { allEvents <- rbind(allEvents,transferEvents) }
  if(nrow(allEvents) > 0) {allEvents <- allEvents[order(allEvents$time),]}
  
  ###### BUILDING MODEL ######
  ###### BUILDING MODEL ######
  ###### BUILDING MODEL ######
  ###### BUILDING MODEL ######
  ###### BUILDING MODEL ######
  
  if(nrow(allEvents)>0) {
    model <- mparse(transitions = transitions, compartments = compartments, events = allEvents, E=E,N=Nmat,
                    gdata = gdata,ldata=ldata, u0 = u0, v0=v0, tspan = tspan, pts_fun = pts_fun)
  } else {
    model <- mparse(transitions = transitions, compartments = compartments,
                    gdata = gdata,ldata=ldata, u0 = u0, v0=v0, tspan = tspan, pts_fun = pts_fun)
  }
  
  ###### RUNNING MODEL ######
  ###### RUNNING MODEL ######
  ###### RUNNING MODEL ######
  ###### RUNNING MODEL ######
  ###### RUNNING MODEL ######
  
  # set.seed(123)
  result <- run(model)
  
  ###### PLOTTING RESULT ######
  ###### PLOTTING RESULT ######
  ###### PLOTTING RESULT ######
  ###### PLOTTING RESULT ######
  ###### PLOTTING RESULT ######
  
  if(is.null(folderPath)) {setwd("./trajectories")}
  
  # drawing plot of active infectious
  trajPlotInfections <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("preSymps_preSympt_Symps_Sympt_postSymps_postSympt_aSymps_aSympt_Isos_Isot"),           # compartments that will be summed and plotted
    newI = FALSE,              # Logical to plot new infections
    groups = plotGroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                     # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    includeMedian = TRUE,                      # Whether to include the median
    confIntv = confIntv,                       # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    twoPeriod = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = dataTitle,                    # Title of data
    miny = miny,                             # optional minimum for y-axis; NULL uses built-in
    maxy = maxy,                             # optional maximum for y-axis; NULL uses built-in
    titleString = paste0("Active Infections; ",titleString),        # plot parameter: Title of plot
    subtitleString = subtitleString,           # plot parameter: subtitle of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )
  
  # drawing plot of cumulative infectious
  trajPlotCumI <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("cumIs_cumIt"),           # compartments that will be summed and plotted
    newI = FALSE,              # Logical to plot new infections
    groups = plotGroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                     # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    includeMedian = TRUE,                      # Whether to include the median
    confIntv = confIntv,                       # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    twoPeriod = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = "Cumulative infections",                    # Title of data
    titleString = paste0("Cumulative infections; ",titleString),        # plot parameter: Title of plot
    subtitleString = subtitleString,           # plot parameter: subtitle of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )
  
  trajPlotNewI <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("cumIs_cumIt"),           # compartments that will be summed and plotted
    newI = TRUE,              # Logical to plot new infections
    groups = plotGroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                     # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    includeMedian = TRUE,                      # Whether to include the median
    confIntv = confIntv,                       # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    twoPeriod = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = "New infections",                    # Title of data
    titleString = paste0("New infections; ",titleString),        # plot parameter: Title of plot
    subtitleString = subtitleString,           # plot parameter: subtitle of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )
  
  # schoolOpenPlot <- simInfPlottingFunction(
  #   result = result,                           # model result
  #   table = "V",                               # which table: U or V
  #   compts= c("noQuarantineS"),           # compartments that will be summed and plotted
  #   groups = plotGroups,                       # List of groups to aggregate and plot
  #   sumGroups = sumGroups,                     # Automatically sums over all groups
  #   uNames = NULL,                        # list of compartments
  #   vNames = names(v0),                             # list of variables
  #   rollM = rollM,                             # number of days for rolling mean
  #   allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
  #   plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
  #   percentile = percentile,                   # percentile for central (or other) tendency
  #   includeMedian = TRUE,                      # Whether to include the median
  #   confIntv = confIntv,                       # confidence interval for plotting spread
  #   nTM = nodeTrialMat,                        # node-Trial matrix
  #   tS = tspan,                                # length of simulation
  #   enn = N,                                   # number of nodes per trial
  #   nT = numTrials,                            # number of trials in simulation
  #   startDate = startofSimDay,                 # start date of simulation
  #   twoPeriod = TRUE,                            # time span is by hours
  #   dateBreaks = dateBreaks,          # plot parameter: Date axis format
  #   dataTitle = "School in attendance",                    # Title of data
  #   titleString = paste0("School, ",titleString),        # plot parameter: Title of plot
  #   xString = "Date",                          # plot parameter: Title of x axis
  #   yString = "School in attendance",          # plot parameter: Title of y axis
  #   lString = lString,                         # plot parameter: Title of legend
  #   cString = cString                          # plot parameter: Plot caption
  # )
  
  ggsave(paste0(plotName,"I.png"),trajPlotInfections,width=9,height=5,units="in")
  ggsave(paste0(plotName,"cumI.png"),trajPlotCumI,width=9,height=5,units="in")
  ggsave(paste0(plotName,"newI.png"),trajPlotNewI,width=9,height=5,units="in")
  # ggsave(paste0(plotName,"schoolOpen.png"),schoolOpenPlot,width=9,height=5,units="in")
  
  ### Number of days in school/in quarantine
  vNames <- names(v0)
  compts <- c("quarCounterC","quarCounterG","quarCounterS")
  vIndex <- which(vNames %in% compts)
  vIndices <- data.table(outer((0:(N*numTrials-1))*length(vNames),vIndex,FUN="+"))
  names(vIndices) <- compts
  dt <- data.table(apply(vIndices,2,function(x) result@V[x,]))
  dt[,trial:=rep.int(nodeTrialMat$trial,max(tspan))][,node:=rep.int(nodeTrialMat$node,max(tspan))][,nodeGroup:=rep.int(nodeTrialMat$nodeGroup,max(tspan))][,time:=sort(rep.int(tspan,N*numTrials),method="quick")]
  
  quarTable <- dt[,lapply(.SD,max),by=.(node,trial),.SDcols = c("quarCounterC","quarCounterG","quarCounterS")]
  quarTable[,classroom:=rep(1:N,numTrials)]
  # median number quarantines/closures over trial
  quarTableMedian <- quarTable[,lapply(.SD,median),by=.(classroom),.SDcols=c("quarCounterC","quarCounterG","quarCounterS")]
  quarTableMedian[,quarCdays:=pmin(maxT,quarTableMedian$quarCounterC*quarantineDaysClassroom)]
  quarTableMedian[,quarGdays:=pmin(maxT,quarTableMedian$quarCounterG*quarantineDaysGrade)]
  quarTableMedian[,quarSdays:=pmin(maxT,quarTableMedian$quarCounterS*quarantineDaysSchool)]

  # max number of quarantines/closures over trial
  quarTableMax <- quarTable[,lapply(.SD,max),by=.(classroom),.SDcols=c("quarCounterC","quarCounterG","quarCounterS")]
  quarTableMax[,quarCdays:=pmin(maxT,quarTableMax$quarCounterC*quarantineDaysClassroom)]
  quarTableMax[,quarGdays:=pmin(maxT,quarTableMax$quarCounterG*quarantineDaysGrade)]
  quarTableMax[,quarSdays:=pmin(maxT,quarTableMax$quarCounterS*quarantineDaysSchool)]
  
  sink(paste0("quarTableMedian.",simID,".txt"))
  print(data.frame(quarTableMedian))
  sink()
  sink(paste0("quarTableMax.",simID,".txt"))
  print(data.frame(quarTableMax))
  sink()

  return(result)
}