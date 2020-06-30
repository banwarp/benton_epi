# covid19_SimInf_06.20.2020.R

# copyright Peter Banwarth Benton County Health Department 2020

# Using package SimInf to model stochastic disease spread
# https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
# https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
# https://github.com/stewid/SimInf

# Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
# Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

# Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
# International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

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
  maxT = 270,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 18,                                  # Number of classrooms
  nodeGroupList = NULL,                    # List of group IDs for node groups. Default is 1 group of nodes.
  Symp0nodeGroups = NULL,                  # Node groups where initial infectious are distributed
  Trans0nodes = NULL,                      # Nodes group where transitory individuals are distributed
  S0Pops = 25*18,                          # Initial number of susceptible - stationary
  preSymp0Pops = 0,                        # Initial number of presymptomatic - stationary
  Symp0Pops = 0,                           # Initial number of symptomatic and asymptomatic - stationary
  postSymp0Pops = 0,                       # Initial number of post-symptomatic - stationary
  R0Pops = 0,                              # Initial number of recovered - stationary
  Im0Pops = 0,                             # Initial number of immune - stationary
  Iso0Pops = 0,                            # Initial number of isolated - stationary
  S0Popt = 10,                             # Initial number of susceptible - transitory
  preSymp0Popt = 0,                        # Initial number of pre-symptomatic - transitory
  Symp0Popt = 0,                           # Initial number of symptomatic and asymptomatic - transitory
  postSymp0Popt = 0,                       # Initial number of post-symptomatic - transitory
  R0Popt = 0,                              # Initial number of recovered - transitory
  Im0Popt = 0,                             # Initial number of immune - transitory
  Iso0Popt = 0,                            # Initial number of isolated - transitory
  
  ### gdata parameters
  R0pre = .5,                              # Basic reproduction number for presypmtomatic
  R0symp = 2,                              # Basic reproduction number for symptomatic
  R0post = .5,                             # Basic reproduction number for postsymptomatic
  R0asymp = .8,                            # Basic reproduction number for asymptomatic
  studentTeacherDiff = rep(1,3),           # Differential R0 for stationary-stationary, transitory-transitory, and stationary-transitory interactions
  preSympPeriod = 1/3,                     # Reciprocal of presymptomatic period
  symptomaticPeriod = 1/4,                 # Reciprocal of symptomatic period
  postSympPeriod = 1/6,                    # Reciprocal of post-symptomatic period
  aSympPeriod = 1/6,                       # Reciprocal of asymptomatic period
  isoPeriod = 1/10,                        # Reciprocal of isolation period
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = 1/100,                   # Reciprocal of temporary immunity period, after which R becomes Im or S
  
  ### ldata parameters
  R0Spread = .1,                           # Uniform variation in R0 across trials (measured as %change from R0)
  
  ### continuous initialized parameters (v0)
  symptomaticProp = .65,                   # Average proportion of infectious who are symptomatic
  symptomaticVariance = 0,                 # Uniform variance of symptomatic proportion
  schoolDay = 1,                           # Indicator for school day: 1 = school day, 0 = evening, list for each node or single entry
  schoolWeek = 1,                          # Indicator for school week: 1 = week, 0 = weekend, list for each node or single entry
  schoolTerm = 1,                          # Indicator for school term: 1 = school term, 0 = break, list for each node or single entry
  noQuarantine = 1,                        # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### pts_fun parameters
  quarantineType = "individual",           # Is quarantine implemented: individual, classroom, grade, or school
  cosAmp = 0.25,                           # Amplitude of seasonal variation in beta
  quarantineDays = 1,                      # Length in days of quarantine
  preDetectSuccess = 0,                    # Probability of observing a pre-symptomatic individual
  sympDetectSuccess = 0,                   # Probability of observing a symptomatic individual
  postDetectSuccess = 0,                   # Probability of observing a post-symptomatic individual
  asympDetectSuccess = 0,                  # Probability of observing an asymptomatic individual
  detectionProbabilityVar = 0,             # Proportionate variable decrease in detection probability - only decreases probability, never increases
  detectSuccessTimeMu = 0,                 # First parameter for function to increase probability of observing as time goes on; can be linear, sigmoid, or uniform
  detectSuccessTimeSig = 0,                # Second parameter for function to increase probability of observing as time goes on; can be linear, sigmoid, or uniform
  detectSuccessLowerBound = 1,             # Lower bound of probability of detection
  detectSuccessUpperBound = 1,             # Upper bound of probability of detection
  noSchoolFactor = .2,                     # Factor to reduce beta when outside of school
  nightFactor = .05,                       # Transmission rate at night; incorporates afternoon also
  quarantineFactor = .1,                   # Factor to reduce beta when quarantine
  
  ### *day* specified infection events parameters - list to allow multiple events
  infectionEventsDF = NULL,                # pre-built infection spreader events data frame if the events will be preset before the script is run
  infectionCount = c(),                    # Number of infections caused by the infection spreader
  infectionNodes = c(),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c(),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(),                   # Symmetric spread in days of infection spreader infections
  infectionCohort = c(),                   # Infection cohort: 1=stationary or 2=transitory
  
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
  plotCompList = "cumI",                   # List of compartments that will be plotted
  rollM = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  allTraj = TRUE,                          # Logical if all simulation trajectories are plotted or just median and spread
  plotRandomTrajs = 5,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  percentile = .5,                         # Percentile of point-in-time simulations to plot, default to median
  confIntv = .95,                          # two-sided confidence interval for plotting spread
  plotGroups = NULL,                       # Which node groups to plot. NULL plots all
  sumGroups = TRUE,                        # Whether to sum the node groups.
  dateBreaks = "1 month",                  # Plot parameter, x-axis displays months
  dataTitle = "Infections",                # Title of data that is plotted
  titleString = "Generic Title",           # Title of plot
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
  
  set.seed(123)
  
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
  
  # if(!is.null(folderPath)) {setwd(folderPath)}
  
  if(is.null(folderPath)) {setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")}
  source("transitionsScriptSchool_06.25.2020.R")
  source("eventFunctionsSchool_06.20.2020.R")
  source("simInfPlottingFunction_06.20.2020.R")
  source("pts_funScript_school_06.30.2020.R")
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

  # Time span in 24 hour days
  tspan <- 1:(maxT*24)
  
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
                    "cumIt"      # cumulative transitory Infections (used for final tabulation)
                    )
  
  # number of compartments that don't participate in stationary transfers and in transitory transfers
  xCompt <- 4 # Isos, Isot, cumIs, cumIt
  
  # More model parameters
  # global parameters
  gdata = data.frame(
    betaPre = R0pre*preSympPeriod/24,          # Transmission rate among presypmtomatic
    betaSymp = R0symp*symptomaticPeriod/24,    # Transmission rate among symptomatic
    betaPost = R0post*postSympPeriod/24,       # Transmission rate among post symptomatic
    betaA = R0asymp*aSympPeriod/24,            # Transmission rate among asymptomatic
    preSympPeriod = preSympPeriod/24,          # Reciprocal of pre symptomatic period
    symptomaticPeriod = symptomaticPeriod/24,  # Reciprocal of infectious period
    postSympPeriod = postSympPeriod/24,        # Reciprocal of post symptomatic period
    aSympPeriod = aSympPeriod/24,              # Reciprocal of asymptomatic period
    isoPeriod = isoPeriod/24,                  # Reciprocal of isolation period
    reSuscepRate = reSuscepRate,               # Proportion of recovered who become re-susceptible
    tempImmPeriod = tempImmPeriod/24,          # Reciprocal of recovered period before re-susceptibility
    ssD = studentTeacherDiff[1],               # Differential on R0 for stationary-stationary interactions
    ttD = studentTeacherDiff[2],               # Differential on R0 for transitory-transitory interactions
    stD = studentTeacherDiff[3]                # Differential on R0 for stationary-transitory interactions
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
  
  # rep dayTimes, weekDays, and breakDays if needed and split into different lists
  if(length(dayTimes) == 0) {dayTimes <- list(c(-1,-1))}
  if(length(weekDays) == 0) {weekDays <- list(c(-1,-1))}
  if(length(breakDays) == 0) {breakDays <- list(c(-1,-1))}
  
  maxDayTimes <- max(unlist(lapply(dayTimes,length)))
  maxWeekDays <- max(unlist(lapply(weekDays,length)))
  maxBreakDays <- max(unlist(lapply(breakDays,length)))
  
  if(length(dayTimes) == 1) {dayTimes <- rep(dayTimes,N)}
  if(length(weekDays) == 1) {weekDays <- rep(weekDays,N)}
  if(length(breakDays) == 1) {breakDays <- rep(breakDays,N)}
  
  dayTimeDF <- data.frame(matrix(data=-1,nrow=N,ncol=maxDayTimes))
  for(i in 1:length(dayTimes)) {
     dayTimeDF[i,c(1:length(dayTimes[[i]]))] <- dayTimes[[i]]
  }
  names(dayTimeDF) <- paste0("dt",1:maxDayTimes)
  
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
  
  ldata <- cbind(ldata,dayTimeDF,weekDayDF,breakDayDF)
  
 
  # Transitions
  ##### Individual quarantine
  if(quarantineType == "individual") {
    transitions <- transitionsIndividualIsolation()
  } else {  ##### All other quarantine types
    transitions <-  transitionsGroupQuarantine()
  }
  
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

  # Distribute initial symptomatic individuals
  if(length(Symp0Pops == 1)) {
    Symp0srepeat <- NnumTrials
  } else {
    Symp0srepeat <- numTrials
  }
  
  if(length(Symp0Popt == 1)) {
    Symp0trepeat <- NnumTrials
  } else {
    Symp0trepeat <- numTrials
  }
  
  # Initial state matrix
  u0 <- data.frame(
    Ss = S0s,
    St = S0t,
    preSymps = rep(preSymp0Pops, NnumTrials),
    preSympt = rep(preSymp0Popt, NnumTrials),
    Symps = rep(Symp0Pops*symptomaticProp, Symp0srepeat),
    Sympt = rep(Symp0Popt*symptomaticProp, Symp0trepeat),
    postSymps = rep(postSymp0Pops, NnumTrials),
    postSympt = rep(postSymp0Popt, NnumTrials),
    aSymps = rep(Symp0Pops*(1-symptomaticProp), Symp0srepeat),
    aSympt = rep(Symp0Popt*(1-symptomaticProp), Symp0trepeat),
    Rs = rep(R0Pops, NnumTrials),
    Rt = rep(R0Popt, NnumTrials),
    Ims = rep(Im0Pops, NnumTrials),
    Imt = rep(Im0Popt, NnumTrials),
    Isos = rep(Iso0Pops, NnumTrials),
    Isot = rep(Iso0Popt, NnumTrials),
    cumIs = rep(Symp0Pops, NnumTrials),
    cumIt = rep(Symp0Popt, NnumTrials)
  )
  
  # Initialized continuous variables
  v0 = data.frame(
    outOfSchool = rep(1,NnumTrials),                         # out-of-school factor for beta
    season = rep(1,NnumTrials),                              # seasonality factor for beta
    symptomaticProp = rep(symptomaticProp,NnumTrials),       # proportion of infections that are symptomatic
    schoolDay = rep(schoolDay,NnumTrials),                   # Indicator for school day: 1 = school day, 0 = evening 
    schoolWeek = rep(schoolWeek,NnumTrials),                 # Indicator for school week: 1 = week, 0 = weekend
    schoolTerm = rep(schoolTerm,NnumTrials),                 # Indicator for school term: 1 = school term, 0 = break
    noQuarantine = rep(noQuarantine,NnumTrials),             # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
    infectionTimer = rep(0,NnumTrials),                      # Timer for length that at least one infection is present in node
    quarTimer = rep(0,NnumTrials),                           # Timer for length of quarantine
    quarCounter = rep(0,NnumTrials),                         # Counter for number of quarantine events
    preDetectSuccess = rep(preDetectSuccess,NnumTrials),     # Probability of detecting a pre-symptomatic individual
    preDetectFailure = rep(1-preDetectSuccess,NnumTrials),   # Probability of failing to detect a pre-symptomatic individual
    sympDetectSuccess = rep(sympDetectSuccess,NnumTrials),   # Probability of detecting a symptomatic individual
    sympDetectFailure = rep(1-sympDetectSuccess,NnumTrials), # Probability of failing to detect a symptomatic individual
    postDetectSuccess = rep(postDetectSuccess,NnumTrials),   # Probability of detecting a post-symptomatic individual
    postDetectFailure = rep(1-postDetectSuccess,NnumTrials), # Probability of failing to detect a post-symptomatic individual
    asympDetectSuccess = rep(asympDetectSuccess,NnumTrials),  # Probability of detecting an asymptomatic individual
    asympDetectFailure = rep(1-asympDetectSuccess,NnumTrials) # Probability of failing to detect an asymptomatic individual
  )

  # building pts_fun
  if(quarantineType == "individual") {
    pts_fun <- pts_funScriptIndividualQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      detVar = detectionProbabilityVar,       # Proportionate variable decrease in probability
      dTimes = maxDayTimes,                   # number of school day time start/stops
      wDays = maxWeekDays,                    # number of school week start/stops
      bDays = maxBreakDays,                   # number of school break start/stops
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor                     # proportionate reduction in beta at night
    )
  } else if(quarantineType == "classroom") {
    pts_fun <- pts_funScriptClassQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      dsTimeMu = detectSuccessTimeMu,         # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeSig = detectSuccessTimeSig,       # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeLB = detectSuccessLowerBound,     # Lower bound on probability of detecting as time increases
      dsTimeUB = detectSuccessUpperBound,     # Upper bound on probability of detecting as time increases
      qDays = quarantineDays,                 # Length in days of quarantine
      dTimes = maxDayTimes,                   # number of school day time start/stops
      wDays = maxWeekDays,                    # number of school week start/stops
      bDays = maxBreakDays,                   # number of school break start/stops
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor,                    # proportionate reduction in beta at night
      quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
    )
  } else if(quarantineType == "grade") {
    pts_fun <- pts_funScriptGradeQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      dsTimeMu = detectSuccessTimeMu,         # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeSig = detectSuccessTimeSig,       # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeLB = detectSuccessLowerBound,     # Lower bound on probability of detecting as time increases
      dsTimeUB = detectSuccessUpperBound,     # Upper bound on probability of detecting as time increases
      qDays = quarantineDays,                 # Length in days of quarantine
      dTimes = maxDayTimes,                   # number of school day time start/stops
      wDays = maxWeekDays,                    # number of school week start/stops
      bDays = maxBreakDays,                   # number of school break start/stops
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor,                    # proportionate reduction in beta at night
      quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
    )
  } else if(quarantineType == "school") {
    pts_fun <- pts_funScriptSchoolQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      dsTimeMu = detectSuccessTimeMu,         # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeSig = detectSuccessTimeSig,       # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeLB = detectSuccessLowerBound,     # Lower bound on probability of detecting as time increases
      dsTimeUB = detectSuccessUpperBound,     # Upper bound on probability of detecting as time increases
      qDays = quarantineDays,                 # Length in days of quarantine
      dTimes = maxDayTimes,                   # number of school day time start/stops
      wDays = maxWeekDays,                    # number of school week start/stops
      bDays = maxBreakDays,                   # number of school break start/stops
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor,                    # proportionate reduction in beta at night
      quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
    )
  }
  
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
                                   cohortList = infectionCohort,
                                   startDay=startofSimDay,
                                   nT = numTrials,
                                   enn = N)
      infectionEvents <- bind_rows(infectionEventList[lengths(infectionEventList) != 0])
    }
  } else {infectionEvents <- infectionEventsDF}
  
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
  if(max(inGroupTransferRate,outGroupTransferRate)>0 | !is.null(transferEventsDF)) { allEvents <- rbind(allEvents,transferEvents) }
  
  
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
    compts= c("preSymps_preSympt_Symps_Sympt_postSymps_postSympt_aSymps_aSympt"),           # compartments that will be summed and plotted
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
    byHour = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = dataTitle,                    # Title of data
    titleString = paste0("Active infections in ",titleString),        # plot parameter: Title of plot
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
    byHour = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = "Cumulative infections",                    # Title of data
    titleString = paste0("Active infections in ",titleString),        # plot parameter: Title of plot
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
    byHour = TRUE,                            # time span is by hours
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    dataTitle = "Cumulative infections",                    # Title of data
    titleString = paste0("Active infections in ",titleString),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )
  
  ggsave(paste0(plotName,"I.png"),trajPlotInfections,width=9,height=5,units="in")
  ggsave(paste0(plotName,"cumI.png"),trajPlotCumI,width=9,height=5,units="in")
  ggsave(paste0(plotName,"newI.png"),trajPlotNewI,width=9,height=5,units="in")
  
  return(result)
}
