# covid19_SimInf_ageGroups06.16.2020.R

##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################
##################### NEED TO TEST THIS SCRIPT #######################

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
# Added in functionality to represent different age groups. 0-17; 18-29, 30-64, 65+
# This will allow for different mixing rates, different hospitalization rates; different mortality rates, etc.
# Removed natural birth rate and death rate to cut down on number of transitions
# Removed the random assignment of S0 to nodes; everything is uniform

# changes from covid19_SimInf_06.06.2020.R
# Added in functionality to represent variable transmission rates: 80% of new infections caused by 20% of infectious; 20% of new caused by 80%, with variability

# changes from covid19_SimInf_05.30.2020.R
# Added post-infectious compartment to represent monitoring symptoms. Post infectious captures all exposed and a proportion of infectious.
# Added unknown-infectious compartment to represent infections that remain unknown
# Changed pts_function to calculate prevalence based only on post-infectious (i.e. observed infections).
# Removed Rtee from local data since it is too hard to calculate exactly and can be reproduced from phi and global R0I, R0P, R0U.
# Added option to pass data frames for different events into the function instead of always using the built in data frame build
# Added option to plot different percentiles instead of just median
# updated kbDates to reflect reality (phase 1 on 5/15, phase 2 on 6/05, moved up phase 3 to anticipated 7/03)
# updated kbPhase = 1 to reflect starting 6/01/2020 in phase 1
# updated simDate = "2020-06-01"
# added kbPhase parameter to allow for updated phase of reopening

# changes from covid19_SimInf_5.16.2020.R
# Changed the parachute distribution function from a janky chi-square to beta distribution with user-defined shape parameters.
# Removed limit on number of trials
# Added an option for uniform distribution of initial susceptibles among all nodes in a trial

# changes from covid19_SimInf_5.14.2020.R
# Added option to choose which nodeGroups get initial infectious and to limit parachuter/super-spreader/massEntry events to certain nodeGroups

# changes from covid19_SimInf_5.12.2020.R
# Added node groups to represent structured sub-populations (e.g. cities, communities with higher mixing, etc.)
# Fixed error in initial distribution of infectious - old versions distributed uniformly instead of according to maxInodeProp

# changes from covid19_SimInf_5.10.2020.R
# wrapped into a function for more efficient scenario generation

# changes from covid19_SimInf_5.8.2020.R
# added hospitalized compartment
# changed death transition from I -> M to H --> M and I -> M to represent that most fatalities occur in individuals who were hospitalized
# improve external (i.e. node-to-node) transfer algorithm
### add option to cluster nodes (higher inGroupTransferRate within nodes compared to across nodes)

# changes from covid19_SimInf_5.6.2020.R
# added option for a "super spreader" through the events
# changed file name for massEntryEVentFunction script to massEventFunction, and added super spreader function as another function in that script
# added option to plot I+Is to represent all infected individuals

# changes from covid19_SimInf_5.5.200.R
# added isolation = Is compartment that exposed can enter to represent enhanced contact tracing
# added cumulative infections = cumI compartment
# added code to simInfPlottingFunction to show daily new infections
# changed plotting function to be able to plot continuous variables as well as compartments

# changes from covid19_SimInf_5.1.2020.R
# reorganized parameters
#   moved constant parameters from V0 to ldata if they varied across trial
#   moved constant parameters from v0/ldata to gdata/pts_function parameters if they didn't vary across trial
#   requires significant changes in pts_fun
# added in pdDecay - decay parameter to represent the slow decay of physical distancing as people return to normal interactions
# added confidence interval parameter for plotting different spreads

# changes from covid19_SimInf_4.30.2020.R
# changed R0 and infectiousPeriod to match IDM 4/22 data
# increased baseline phi to represent impact of physical distancing and contact tracing even in absence of interventions
# created parmsList to organize parameters

# changes from covid19_SimInf_4.24.2020.R
# Split phiMovement into an up movement and a down movement, which seems more realistic to me.

# changes from covid19_SimInf_4.24.2020.R
# Starting simulation with high intervention to match reality, using continuous variable KB in pts_function

# changes from covid19_SimInf_4.23.2020.R
# adding local data for beta to introduce more variability

# changes from covid19_SimInf_4.18.2020.R
# moved pts_fun to a separate script for space and annotation
# added more massEntry functionality
# added local data parameters for more variability

# changes from covid19_SimInf_4.9.2020.R
# removed all the efforts before Effort 15

covidWrapper <- function(
  ### Simulation parameters
  folderPath = NULL,                       # folder path for subroutines and output
  simID = "covid.mo.day.2020.Generic",     # Simulation ID
  simDate = "2020-06-01",                  # Start date of simulation
  maxT = 10,                              # Length of simulation in days
  numTrials = 2,                           # Number of trials in simulation
  
  ### population parameters
  trialPop = 1000,                         # Total population in each trial
  N = 3,                                  # Number of population nodes
  nodeGroupList = NULL,                    # List of group IDs for node groups. Default is 1 group of nodes.
  S0Prop = c(.25,.25,.25,.25),                 # Initial population: 0-17
  I0Prop = rep(.1,4),                   # Initial proportion of infectious: 0-17
  maxINodeProp = 1/10,                     # Maximum proportion of nodes that intially have one or more infectious
  I0nodeGroups = list(NULL,NULL,NULL,NULL), # Node groups where initial infectious are distributed by age group
  E0Pop = rep(0,4),                               # Initial number of exposed; all ages
  uI0Pop = rep(0,4),                              # Initial number of unknown-infectious
  R0Pop = rep(0,4),                               # Initial number of recovered
  Im0Pop = rep(0,4),                              # Initial number of immune
  pI0Pop = rep(0,4),                              # Initial number of post-infectious
  H0Pop = rep(0,4),                               # Initial number of hospitalized
  Is0Pop = rep(0,4),                              # Initial number of isolated
  M0Pop = rep(0,4),                               # Initial number of dead
  
  ### gdata parameters
  R0I = 1.9,                               # Basic reproduction number for first infectiousPeriod days
  R0P = .125,                              # Basic reproduction number for post infectious (after successful monitoring)
  R0U = 1,                                 # Basic reproduction number for remaining infectious period if monitoring fails
  RFactor = rep(1,10),                     # Factor to adjust the reproductive number within and between different age groups.
  initInfectiousPeriod = 1/4,              # Reciprocal of initial infectious period
  postInfectiousPeriod = 1/10,             # Reciprocal of post-infectious period
  unknownInfectiousPeriod = 1/6,           # Reciprocal of remaining, unknown-infectious period
  exposedPeriod = 1/4,                     # Reciprocal of exposed period
  isoRates = rep(.125,4),                  # Proportion of exposed who are identified and isolated before they become infectious; each age group
  monitoringSuccess = rep(.25,4),          # Proportion of infectious whose symptoms are identified and moved to post-infectious; each age group
  hospRateExp = .033,                      # Proportion of isolated that are hospitalized - should match overall population hospitalization rate
  hospRatePost = .16,                      # Proportion of successfully identified cases who are hospitalized - should match observed hospitalization rate
  hospRateFactor = rep(1,4),               # Factor to adjust hospitalization rate for each age group.
  hospPeriod = 1/10,                       # Reciprocal of length of hospitalization period
  nonHospDeathRate = 0,                    # Non-hospitalized fatality rate
  hospDeathRate = .125,                    # Hospitalized fatality rate NOT CASE FATALITY RATE
  deathRateFactor = rep(1,4),              # Factor to adjust mortality rate for each age group.
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = 1/100,                   # Reciprocal of temporary immunity period, after which R becomes Im or S
  
  ### ldata parameters
  R0Spread = .1,                           # Uniform variation in R0 across trials (measured as %change from R0)
  
  ### continuous initialized parameters (v0)
  phi0 = 2.9/.9,                           # Initial phi under the March 23rd stay-at-home orders
  highSpreadProp = .5,                     # Average proportion of infectious who are high spreaders
  highSpreadVariance = 0,                  # Uniform variance of high spreader proportion
  highSpreadRate = .5,                     # Proportion of new cases caused by high spreaders - called rate to differentiate from highSpreadProp; defaults to no difference
  
  ### pts_fun parameters
  cosAmp = 0.25,                           # Amplitude of seasonal variation in beta
  RPhysicalDistancing = 2,                 # Ongoing baseline Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  RNoAction = NULL,                        # Ongoing baseline Rt if no actions at are all taken
  RTarget1 = 1.5,                          # Target for the reduction in R under minor intervention
  RTarget2 = .9,                           # Target for the reduction in R under major intervention
  maxPrev1 = .0005,                         # Maximum prevalence before instituting minor intervention
  maxPrev2 = .001,                         # Maximum prevalence before instituting major intervention
  upDelay = 10,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 28,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = .25,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = .25,                       # Rate at which phi decreases when interventions are lifted
  pdDecay = 30,                            # Number of days until RPhysicalDistancing decays toward RNoAction in the absence of interventions.
  ########################################## Represents gradual relaxation of physical distancing as people return to normal.
  ########################################## pdDecay = -1 removes this decay factor
  
  ### other date parameters
  kbDay1 = "2020-05-15",                   # Date of first phase of lifting stay-at-home orders
  kbDay2 = "2020-06-05",                   # Date of second phase of lifting stay-at-home orders
  kbDay3 = "2020-07-03",                   # Date of third phase of lifting stay-at-home orders
  kbPhase = 1,                             # Phase of reopening, defaults to full stay-at-home orders
  switchOffPolicies = 0,                   # Indicator if intervention policies will cease after a certain day
  switchOffDay = "2020-07-03",             # Date intervention policies will cease if indicator == 1
  
  ### parachuter event parameters
  paraEventsDF = NULL,                     # pre-built parachuter events data frame if the events will be preset before the script is run
  paraMu = 1,                              # First shape parameter for beta function for timing of parachute events
  paraSig = 1,                             # Second shape parameter for beta function for timing of parachute events
  parachuteRate = 1/21,                    # Reciprocal of expected waiting time for a parachute event
  parachuteNum = 1,                        # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### transfer event parameters
  transferEventsDF = NULL,                 # pre-built transfer events data frame if the events will be preset before the script is run
  inGroupTransferRate = list(1/7,1/7,1/7,1/7),               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = list(.1,.1,.1,.1),             # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = list(.1,.1,.1,.1),             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = list(.3,.3,.3,.3),             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = list(1/15,1/15,1/15,1/15),             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = list(.05,.05,.05,.05),           # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = list(.05,.05,.05,.05),           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = list(.15,.15,.15,.15),           # Maximum proportion of node population that transfers in each event
  
  ### super spreader parameters - list to allow multiple events
  superEventsDF = NULL,                    # pre-built super spreader events data frame if the events will be preset before the script is run
  superInfections = c(),                   # Number of infections caused by the super spreader
  superNodes = c(),                        # Number of nodes that the super spreader contacts
  superNodeGroups = NULL,                  # Which node groups the super spreader contacts. Must use list() syntax for multiple events
  superDate = c(),                         # Date the super spreader lands. Date can also be numeric i.e. 200
  superSpread = c(),                       # Symmetric spread in days of super spreader infections
  superAgeGroup = c(),                     # Which age group is infected; 1=k,2=y,3=a,4=s
  
  ### massEntry parameters
  massEntryEventsDF = NULL,                # pre-built mass entry events data frame if the events will be preset before the script is run
  massEntryReturnDate = "2020-09-21",      # Date of mass Entry
  massEntryReturnSpread = 3,               # Symmetric spread of days to spread out mass entry
  massEntryPop = 100,                    # Mass entry population
  maxMassEntryNodes = 0,                   # Maximum number of nodes that individuals enter. Set to 0 to remove massEntry event
  massEntryNodeGroups = NULL,              # Which node groups the individuals enter
  mSProp = rep(.9,4)/4,                      # Proportion of individuals who are susceptible
  mEProp = rep(.0001,4)/4,                          # Proportion of individuals who are exposed
  mIProp = rep(.001,4)/4,                           # Proportion of individuals who are infectious
  muIProp = rep(0,4)/4,                             # Proportion of individuals who are unknown-infectious; other compartments are calculated
  
  ### plot parameters
  plotCompList = "cumI",                   # List of compartments that will be plotted
  rollM = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  allTraj = FALSE,                         # Logical if all simulation trajectories are plotted or just median and spread
  plotRandomTrajs = 0,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  percentile = .5,                         # Percentile of point-in-time simulations to plot, default to median
  confIntv = .95,                          # two-sided confidence interval for plotting spread
  plotGroups = NULL,                       # Which node groups to plot. NULL plots all
  sumGroups = TRUE,                          # Whether to sum the node groups.
  dateBreaks = "1 month",                  # Plot parameter, x-axis displays months
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
  
  library(ggplot2)
  library(reshape2)
  library(SimInf)
  library(data.table)
  library(dplyr)
  library(simFrame)
  library(Matrix)
  library(zoo)
  library(viridis)
  
  if(!is.null(folderPath)) {setwd(folderPath)}
  
  if(is.null(folderPath)) {setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")}
  source("eventFunctionsAgeGroups_06.19.2020.R")
  source("simInfPlottingFunction_05.30.2020.R")
  source("pts_funScriptAgeGroups06.19.2020.R")
  source("transitionsScriptAgeGroups06.19.2020.R")
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
  # Date current physical distancing is relaxed
  kbDay1 <- as.numeric(as.Date(kbDay1)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
  kbDay2 <- as.numeric(as.Date(kbDay2)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
  kbDay3 <- as.numeric(as.Date(kbDay3)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
  
  # Time span
  tspan <- 1:maxT
  
  # additional plotting parameters
  fileName <- paste0("data",simID,".txt")
  plotName <- paste0("plot",simID)
  if(is.null(cString)) {
    cString <- plotName
  }
  
  # Trial and node parameters
  numTrials <- min(100,numTrials) # capping number of trials at 100
  N <- min(1000,N) # capping number of nodes per trial at 1000
  if(is.null(nodeGroupList)) {nodeGroupList <- rep(1,N)} # default is one nodeGroup
  numGroups <- length(unique(nodeGroupList))
  NnumTrials <- N*numTrials
  nodeTrialMat <- data.table(node=c(1:(NnumTrials)),
                             nodeGroup = rep(nodeGroupList,numTrials),
                             trial=rep(1:numTrials,each=N))
  NList <- split(nodeTrialMat$node,nodeTrialMat$trial)
  nodeAndGroupList <- split(nodeTrialMat,nodeTrialMat$nodeGroup)
  nodeAndGroupList <- lapply(nodeAndGroupList,function(x) split(x$node,x$trial))
  I0nodeGroups <- lapply(I0nodeGroups, function(x) if(is.null(x)) {
                                                                  x <- unique(nodeGroupList)
                                                                  } else {x <- x} )
  maxNodePop <- floor(trialPop/N) # maximum node population
  
  # Model parameters
  # compartments
  compartments <- paste0(rep(c("S","E","hI","lI","uhI","ulI","R","Im","pI","H","cumI","M"),each=4),rep(c("k","y","a","s"),12))
  # Susceptible, Exposed, high Infectious, low Infectious, unknown high Infectious, unknown low Infectious,
  # Recovered, Immune, post Infectious, Hospitalized, cumulative Infections, Mortality; kid, youngadult, adult senior
  
  # compartments that don't participate in transfers
  xCompt <- 16
  
  # More model parameters
  # global parameters
  gdata = data.frame(
    betaI = R0I*initInfectiousPeriod,                      # Transmission rate among initial infectious = R0/infectiousPeriod
    betaP = R0P*postInfectiousPeriod,                      # Transmission rate among post infectious and isolated
    betaU = R0U*unknownInfectiousPeriod,                   # Transmission rate among unknown infectious
    spreadRatio = highSpreadRate/(1-highSpreadRate),       # Ratio of transmissibility of high spreaders to low spreaders
    initInfectiousPeriod = initInfectiousPeriod,           # Reciprocal of infectious period
    postInfectiousPeriod = postInfectiousPeriod,           # Reciprocal of post infectious period and isolation period
    unknownInfectiousPeriod = unknownInfectiousPeriod,     # Reciprocal of unknown infectious period
    exposedPeriod = exposedPeriod,                         # Reciprocal of exposed period
    mSk = monitoringSuccess[1],                            # How successful symptom monitoring is to transition from hI/lI to pI: 0-17
    mSy = monitoringSuccess[2],                            # How successful symptom monitoring is to transition from hI/lI to pI: 18-29
    mSa = monitoringSuccess[3],                            # How successful symptom monitoring is to transition from hI/lI to pI: 30-64
    mSs = monitoringSuccess[4],                            # How successful symptom monitoring is to transition from hI/lI to pI: 65+
    isoK = isoRates[1],                                    # Acknowledges that all hospitalized will be isolated: 0-17
    isoY = isoRates[2],                                    # Acknowledges that all hospitalized will be isolated: 18-29
    isoA = isoRates[3],                                    # Acknowledges that all hospitalized will be isolated: 30-64
    isoS = isoRates[4],                                    # Acknowledges that all hospitalized will be isolated: 65+
    hospExpK = hospRateExp*hospRateFactor[1],              # Hospitalization rate among exposed (lower): 0-17
    hospExpY = hospRateExp*hospRateFactor[2],              # Hospitalization rate among exposed (lower): 18-29
    hospExpA = hospRateExp*hospRateFactor[3],              # Hospitalization rate among exposed (lower): 30-64
    hospExpS = hospRateExp*hospRateFactor[4],              # Hospitalization rate among exposed (lower): 65+
    hospPostK = hospRatePost*hospRateFactor[1],            # Hospitalization rate among symptomatic who are identified (higher): 0-17
    hospPostY = hospRatePost*hospRateFactor[2],            # Hospitalization rate among symptomatic who are identified (higher): 18-29
    hospPostA = hospRatePost*hospRateFactor[3],            # Hospitalization rate among symptomatic who are identified (higher): 30-64
    hospPostS = hospRatePost*hospRateFactor[4],            # Hospitalization rate among symptomatic who are identified (higher): 65+
    hospPeriod = hospPeriod,                               # Length of hospitalization
    nonHospDeathK = nonHospDeathRate*deathRateFactor[1],   # non-hospitalized fatality rate: 0-17
    nonHospDeathY = nonHospDeathRate*deathRateFactor[2],   # non-hospitalized fatality rate: 18-29
    nonHospDeathA = nonHospDeathRate*deathRateFactor[3],   # non-hospitalized fatality rate: 30-64
    nonHospDeathS = nonHospDeathRate*deathRateFactor[4],   # non-hospitalized fatality rate: 65+
    hospDeathK = hospDeathRate*deathRateFactor[1],     # Hospitalized fatality rate NOT CASE FATALITY RATE: 0-17
    hospDeathY = hospDeathRate*deathRateFactor[2],     # Hospitalized fatality rate NOT CASE FATALITY RATE: 18-29
    hospDeathA = hospDeathRate*deathRateFactor[3],     # Hospitalized fatality rate NOT CASE FATALITY RATE: 30-64
    hospDeathS = hospDeathRate*deathRateFactor[4],     # Hospitalized fatality rate NOT CASE FATALITY RATE: 65+
    reSuscepRate = reSuscepRate,                           # Proportion of recovered who become re-susceptible
    tempImmPeriod = tempImmPeriod,                         # Reciprocal of recovered period before re-susceptibility
    kkf = RFactor[1],                                     # kid-kid factor
    kyf = RFactor[2],                                     # kid-youngadult factor
    kaf = RFactor[3],                                     # kid-adult factor
    ksf = RFactor[4],                                     # kid-senior factor
    yyf = RFactor[5],                                     # youngadult-youngadult factor
    yaf = RFactor[6],                                     # youngadult-adult factor
    ysf = RFactor[7],                                     # youngadult-senior factor
    aaf = RFactor[8],                                     # adult-adult factor
    asf = RFactor[9],                                     # adult-senior factor
    ssf = RFactor[10]                                    # senior-senior factor
  )
  
  # local parameters
  ldata <- data.frame(
    firstNode = sort(rep.int(0:(numTrials-1)*N,N),method="quick"), # first node in a trial, for computing trial prevalence
    betaRandomizer = rep(runif(numTrials,min=1-R0Spread,max=1+R0Spread),each=N)  # randomizer for beta or R0
  )
 
  # Transitions
  transitions <- transitionsScript()
  
  # Initial state
  # Uniformly distribute population between nodes
  S0Prop <- S0Prop/sum(S0Prop)
  S0 <- lapply(S0Prop, function(x) round(rep(x*trialPop/N,NnumTrials),0))
  eligibleInodes <- lapply(I0nodeGroups, function(x) nodeTrialMat[which(nodeTrialMat$nodeGroup %in% x & nodeTrialMat$trial==1),"node"]$node)
  
  InodeFunction <- function(ageGroupInodes,soprop,ioprop) {
    Inodes <- lapply(0:(numTrials-1),function(x) (x*N)+sample(ageGroupInodes,ceiling(maxINodeProp*length(ageGroupInodes))))
    I0indices <- lapply(Inodes,function(x) if(length(x)>1){sample(x,rbinom(1,2*round(trialPop*soprop*ioprop,0),.5),replace=TRUE)}
                        else{rep(x,round(trialPop*soprop*ioprop,0))})
    I0indices <- unlist(lapply(I0indices,function(x) table(x)))
    I0 <- matrix(0,nrow=NnumTrials)
    I0[as.numeric(names(I0indices))] <- I0indices
    return(I0)
  }
  
  I0 <- lapply(1:length(eligibleInodes), function(x) InodeFunction(eligibleInodes[[x]],S0Prop[x],I0Prop[x]))

  # Initial state data frame
  u0 <- data.frame(
    Sk = S0[[1]],
    Sy = S0[[2]],
    Sa = S0[[3]],
    Ss = S0[[4]],
    Ek = rep(E0Pop[1],NnumTrials),
    Ey = rep(E0Pop[2],NnumTrials),
    Ea = rep(E0Pop[3],NnumTrials),
    Es = rep(E0Pop[4],NnumTrials),
    hIk = floor(I0[[1]]*highSpreadProp),
    hIy = floor(I0[[2]]*highSpreadProp),
    hIa = floor(I0[[3]]*highSpreadProp),
    hIs = floor(I0[[4]]*highSpreadProp),
    lIk = ceiling(I0[[1]]*(1-highSpreadProp)),
    lIy = ceiling(I0[[2]]*(1-highSpreadProp)),
    lIa = ceiling(I0[[3]]*(1-highSpreadProp)),
    lIs = ceiling(I0[[4]]*(1-highSpreadProp)),
    uhIk = rep(round(uI0Pop[1]*highSpreadProp,0),NnumTrials),
    uhIy = rep(round(uI0Pop[2]*highSpreadProp,0),NnumTrials),
    uhIa = rep(round(uI0Pop[3]*highSpreadProp,0),NnumTrials),
    uhIs = rep(round(uI0Pop[4]*highSpreadProp,0),NnumTrials),
    ulIk = rep(round(uI0Pop[1]*(1-highSpreadProp),0),NnumTrials),
    ulIy = rep(round(uI0Pop[2]*(1-highSpreadProp),0),NnumTrials),
    ulIa = rep(round(uI0Pop[3]*(1-highSpreadProp),0),NnumTrials),
    ulIs = rep(round(uI0Pop[4]*(1-highSpreadProp),0),NnumTrials),
    Rk = rep(R0Pop[1],NnumTrials),
    Ry = rep(R0Pop[2],NnumTrials),
    Ra = rep(R0Pop[3],NnumTrials),
    Rs = rep(R0Pop[4],NnumTrials),
    Imk = rep(Im0Pop[1],NnumTrials),
    Imy = rep(Im0Pop[2],NnumTrials),
    Ima = rep(Im0Pop[3],NnumTrials),
    Ims = rep(Im0Pop[4],NnumTrials),
    pIk = rep(pI0Pop[1],NnumTrials),
    pIy = rep(pI0Pop[2],NnumTrials),
    pIa = rep(pI0Pop[3],NnumTrials),
    pIs = rep(pI0Pop[4],NnumTrials),
    Hk = rep(H0Pop[1],NnumTrials),
    Hy = rep(H0Pop[2],NnumTrials),
    Ha = rep(H0Pop[3],NnumTrials),
    Hs = rep(H0Pop[4],NnumTrials),
    cumIk = I0[[1]],
    cumIy = I0[[2]],
    cumIa = I0[[3]],
    cumIs = I0[[4]],
    Mk = rep(M0Pop[1],NnumTrials),
    My = rep(M0Pop[2],NnumTrials),
    Ma = rep(M0Pop[3],NnumTrials),
    Ms = rep(M0Pop[4],NnumTrials)
  )
  
  # Initialized continuous variables
  v0 = data.frame(
    phi = rep(phi0,NnumTrials),                              # initial beta reduction factor (larger phi = more reduction)
    prevUp01 = rep(upDelay+1,NnumTrials),                    # initialized variable for capturing delay in response when prevalence increases
    prevUp12 = rep(upDelay+1,NnumTrials),                    # variable for capturing delay in response when prevalence increases
    prevDown21 = rep(downDelay+1,NnumTrials),                # variable for capturing delay in response when prevalence decreases
    prevDown10 = rep(downDelay+1,NnumTrials),                # variable for capturing delay in response when prevalence decreases
    kbPhase = rep(kbPhase,NnumTrials),                       # 0 means current physical distancing, 1-3 is phase of lifting.
    policy = rep(1,NnumTrials),                              # logical: 1 means policies can take effect; 0 means they don't
    season = rep(1,NnumTrials),                              # seasonality factor for beta
    observedPrev = rep(0,NnumTrials),                        # prevalence tracker - observed from pI. Does not include hI, lI, uhI or ulI by definition.
    previousState = rep(2,NnumTrials),                       # previous state: 0 = baseline, 1 = low intervention, 2 = high intervention, used for low intervention logic
    pdCounter = rep(0,NnumTrials),                           # counter for pdDecay
    highSpreadProp = rep(highSpreadProp,NnumTrials)          # Average proportion of infectious who are high spreaders
  )
  
  if(is.null(RNoAction)) {RNoAction <- R0I+R0U}
  
  # building pts_fun
  pts_fun <- pts_funScript(
    phiPhysicalDistancing = (R0I+R0U)/RPhysicalDistancing,    # Phi reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
    phiNoAction = (R0I+R0U)/RNoAction,                        # Phi reflecting no actions
    maxPrev1 = maxPrev1,                                      # Maximum prevalence before instituting minor intervention
    maxPrev2 = maxPrev2,                                      # Maximum prevalence before instituting major intervention
    phiFactor1 = (R0I+R0U)/RTarget1,                          # Target for the reduction in R0 under minor intervention
    phiFactor2 = (R0I+R0U)/RTarget2,                          # Target for the reduction in R0 under major intervention
    cosAmp = cosAmp,                                          # Amplitude of seasonal variation in beta
    startDay = startofSimDay,                                 # start of simulation
    kbDay1 = kbDay1,                                          # date of first phase
    kbDay2 = kbDay2,                                          # date of second phase
    kbDay3 = kbDay3,                                          # date of third phase
    enn = N,                                                  # number of nodes in a trial
    numComp = length(compartments),                           # number of compartments
    pICompt = which(compartments=="pIk")-1,                   # which compartment is pI for prevalence monitoring
    prevType = ifelse(maxPrev1 < 1,1,0),                      # prevalence type. 0 = count, 1 = proportion
    upDelay = upDelay,                                        # Number of days after prevalence passes threshold until minor/major intervention
    downDelay = downDelay,                                    # Number of days after prevalence drops below threshold until intervention lifted
    phiMoveUp = phiMoveUp,                                    # Rate at which phi increases when interventions are imposed
    phiMoveDown = phiMoveDown,                                # Rate at which phi decreases when interventions are lifted
    pdDecay = pdDecay,                                        # Rate at which phi decreases toward 1 in the absence of interventions. Represents gradual relaxation of physical distancing
    highSpreadProp0 = highSpreadProp,                         # Baseline high spreader proportion that gets re-randomized each time step
    highSpreadVar = highSpreadVariance,                       # Uniform variance of high spreader proportion
    switchOffPolicies = switchOffPolicies,                    # Logical variable to switch off policies (used for counterfactuals)
    switchOffDay = as.numeric(as.Date(switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
  )
  
  
  
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  
  # create events: occasional infections parachute into the nodes; occasional transfers between the nodes. individuals arrive in a wave
  
  # Events matrix
  # leaves off Is, H, cumI, and M from last column for each age group
  E <- cbind(diag(length(compartments)),
             c(rep(c(1,0,0,0),(length(compartments)-xCompt)/4),rep(0,xCompt)),
             c(rep(c(0,1,0,0),(length(compartments)-xCompt)/4),rep(0,xCompt)),
             c(rep(c(0,0,1,0),(length(compartments)-xCompt)/4),rep(0,xCompt)),
             c(rep(c(0,0,0,1),(length(compartments)-xCompt)/4),rep(0,xCompt)),
             c(rep(1,length(compartments)-xCompt),rep(0,xCompt)))
  dimnames(E) <- list(compartments,c(1:ncol(E)))
  
  # Shift matrix; shifts each age group S to the corresponding age group E
  Nmat <- matrix(c(rep(4,4),rep(0,length(compartments)-4)),ncol=1)
  dimnames(Nmat) <- list(compartments,1)
  
  ###### PARACHUTE EVENTS ######
  # Follows a poisson process
  
  if(is.null(paraEventsDF)) {
  
    # Allowable node groups: all are allowed if none are selected
    if(is.null(parachuteNodeGroups)) {parachuteNodeGroups <- unique(nodeGroupList)}
    
    # number of events in each trial
    numParachuteList <- lapply(1:numTrials,function(x) rpois(1,parachuteRate*max(tspan)))
    parachuteNTM <- nodeTrialMat[nodeTrialMat$nodeGroup %in% parachuteNodeGroups,]
    pNList <- split(parachuteNTM$node,parachuteNTM$trial)
    parachuteSelect <- lapply(numParachuteList, function(x) ups(8,x,prob = c(highSpreadProp*S0Prop,(1-highSpreadProp)*S0Prop),replace=TRUE))
    parachuteSelect <- unlist(lapply(parachuteSelect, function(x) compartments[x+8]))
    parachuteSelect <- unlist(lapply(parachuteSelect,function(x) which(compartments==x)))
    
    if(max(unlist(numParachuteList)) > 0){
      # parachute events data frame
      parachuteEvents <- data.frame(
        event = "enter",
        time = unlist(lapply(numParachuteList, function(x) ceiling(rbeta(x,paraMu,paraSig)*max(tspan)))),
        node = unlist(lapply(c(1:numTrials),function(x) sample(pNList[[x]],numParachuteList[[x]],replace=TRUE))),
        dest = 0,
        n = parachuteNum,
        proportion = 0,
        select = parachuteSelect,
        shift = 0
      )
    }
  } else {parachuteEvents <- paraEventsDF}
  
  ###### TRANSFER EVENTS ######
  if(is.null(transferEventsDF)) {
    
    transferEvents <- data.frame(NULL)
    
    for(i in 1:4) {
      
      iGTR <- inGroupTransferRate[[i]]
      iGTNN <- inGroupTransferNodeNum[[i]]
      iGTMinP <- inGroupTransferMinProp[[i]]
      iGTMaxP <- inGroupTransferMaxProp[[i]]
      oGTR <- outGroupTransferRate[[i]]
      oGTNN <- outGroupTransferNodeNum[[i]]
      oGTMinP <- outGroupTransferMinProp[[i]]
      oGTMaxP <- outGroupTransferMaxProp[[i]]
      
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
      }
    
      ###### IN-GROUP TRANSFER EVENTS ######
      # Follows a poisson distribution
      # Partition transfers by trials and groups
      if(iGTNN[1]<1) {iGTNN <- ceiling(iGTNN*N)} # converts proportion to number
      
      numInTransferList <- lapply(1:numGroups,function(x) 
        lapply(1:numTrials, function(y) rpois(1,iGTR[x]*iGTNN[x]*max(tspan))))
      
      inTransferEventsList <- lapply(1:numGroups, transferFunction,
                                   numTList = numInTransferList,
                                   tS = tspan,
                                   nT = numTrials,
                                   nList = nodeAndGroupList,
                                   minProp = iGTMinP,
                                   maxProp = iGTMaxP,
                                   selectCol = length(compartments)+i,
                                   nTM = nodeTrialMat,
                                   u0S = u0$S,
                                   outLogic = FALSE
                                   )
      inTransferEvents <- bind_rows(inTransferEventsList[lengths(inTransferEventsList) != 0])
       
      ###### OUT-GROUP TRANSFER EVENTS ######
      if(numGroups > 1) {
        if(oGTNN[1]<1) {oGTNN<- ceiling(oGTNN*N)} # converts proportion to number
        
        numOutTransferList <- lapply(1:numGroups,function(x) 
          lapply(1:numTrials, function(y) rpois(1,oGTR[x]*oGTNN[x]*max(tspan))))
        
        outTransferEventsList <- lapply(1:numGroups, transferFunction,
                                        numTList = numOutTransferList,
                                        tS = tspan,
                                        nT = numTrials,
                                        nList = nodeAndGroupList,
                                        minProp = oGTMinP,
                                        maxProp = oGTMaxP,
                                        selectCol = length(compartments)+i,
                                        nTM = nodeTrialMat,
                                        u0S = u0$S,
                                        outLogic = TRUE
                                        )
        outTransferEvents <- bind_rows(outTransferEventsList[lengths(outTransferEventsList) != 0])
        transferEvents <- rbind(transferEvents,inTransferEvents,outTransferEvents)
        rm(outTransferEventsList)
        rm(outTransferEvents)
      } else {transferEvents <- rbind(transferEvents,inTransferEvents)}
      rm(inTransferEventsList)
      rm(inTransferEvents)
    } # end of i in 1:4 loop
  } else {transferEvents <- transferEventsDF}
  
  ###### SUPER SPREADER EVENTS #####
  if(is.null(superEventsDF)) {
  
    if(length(superInfections) > 0){
      
      if(is.null(superNodeGroups)) {superNodeGroups <- lapply(superInfections,function(x) unique(nodeGroupList))}
      
      eligibleSuperNodes <- lapply(superNodeGroups,function(x) nodeTrialMat[which(nodeTrialMat$nodeGroup %in% x & nodeTrialMat$trial==1),"node"]$node)
  
      superNodeList <- lapply(1:length(superNodes),
                              function(x) if(length(eligibleSuperNodes[[x]])>1)
                                {sort(sample(eligibleSuperNodes[[x]],min(length(eligibleSuperNodes[[x]]),superNodes[x])))}
                                else{eligibleSuperNodes[[x]]})
    
      superEventList <- lapply(1:length(superInfections), superFunction,
                               IList=superInfections,
                               ageGroup = superAgeGroup,
                               nodeList=superNodeList,
                               dayList=superDate,
                               spreadList=superSpread,
                               startDay=startofSimDay,
                               nT = numTrials,
                               enn = N)
      superEvents <- bind_rows(superEventList[lengths(superEventList) != 0])
    }
  } else {superEvents <- superEventsDF}
    
  ###### MASS ENTRY EVENTS ######
  
  if(is.null(massEntryEventsDF)) {
  
    # massEntry event parameters - relies on compartments
    if(is.null(massEntryNodeGroups)) {massEntryNodeGroups <- unique(nodeGroupList)}
    
    massEntryReturnDate <- as.numeric(as.Date(massEntryReturnDate)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
    maxMassEntryNodes <- min(length(which(nodeGroupList %in% massEntryNodeGroups)),maxMassEntryNodes) # maximum number of nodes individuals can enter
    mhIProp <- mIProp*highSpreadProp
    mlIProp <- mIProp*(1-highSpreadProp)
    muhIProp <- muIProp*highSpreadProp
    mulIProp <- muIProp*(1-highSpreadProp)
    mRProp <- (1-(mSProp+mEProp+mIProp))*(1-reSuscepRate) # proportion of individuals who are recovered
    mImProp <- (1-(mSProp+mEProp+mIProp))*reSuscepRate # proportion of individuals who are immune
    mpIProp <- rep(0,4)
    mHProp <- rep(0,4)
    mcumIProp <- mIProp
    mMProp <- rep(0,4)
    massEntryPropTable <- data.frame(compartment = compartments,
                                     frac = c(mSProp,mEProp,mhIProp,mlIProp,
                                              muhIProp,mulIProp,mRProp,mImProp,mpIProp,
                                              mHProp,mcumIProp,mMProp))
    
    if(maxMassEntryNodes > 0) {
      # Event for day(s) individuals return, with random returns around massEntryReturnDate
      eligibleMassEntryNodes <- nodeTrialMat[which(nodeTrialMat$nodeGroup %in% massEntryNodeGroups & nodeTrialMat$trial==1),"node"]$node
      massEntryNodeList <- sort(sample(eligibleMassEntryNodes,min(length(eligibleMassEntryNodes),maxMassEntryNodes)))
      
      # exclude isolated, hospitalized,  cumulative, and deceased by length(compartments)-xCompt
      massEntryEventList <- lapply(c(1:(length(compartments)-xCompt)),massEntryEventFunction,
                                 fracTable = massEntryPropTable,
                                 nodeList = massEntryNodeList,
                                 pop = massEntryPop,
                                 nT = numTrials,
                                 enn = N,
                                 sRD = massEntryReturnDate,
                                 sRS = massEntryReturnSpread)
      massEntryEvents <- bind_rows(massEntryEventList[lengths(massEntryEventList) != 0])
      
      massEntryCumIFunction <- function(x, mEE=massEntryEvents, compts=compartments) {
        mECI <- mEE[which(mEE$select %in% which(compts %in% paste0(c("hI","lI","uhI","ulI"),x))),]
        if(nrow(mECI) > 0){
          mECI$select <- which(compts==paste0("cumI",x))
          return(mECI)
        } else {
          return(NULL)
        }
      }
    
      massEntryCumIList <- lapply(c("k","y","a","s"), massEntryCumIFunction)
      massEntryCumI <- bind_rows(massEntryCumIList[lengths(massEntryEventList) != 0])
      massEntryEvents <- rbind(massEntryEvents,massEntryCumI)
    }
  } else {massEntryEvents <- massEntryEventsDF}
  
  ###### ALL EVENTS ######
  allEvents <- data.frame(NULL)
  if(nrow(parachuteEvents) > 0) { allEvents <- rbind(allEvents,parachuteEvents) }
  if(nrow(transferEvents) > 0) { allEvents <- rbind(allEvents,transferEvents) }
  if(nrow(superEvents) > 0 ) { allEvents <- rbind(allEvents,superEvents)}
  if(nrow(massEntryEvents) > 0) { allEvents <- rbind(allEvents,massEntryEvents) }
  
  
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
  
  # drawing plot of infectious +  isolated to get total infected
  trajPlotInfections <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("hI_lI_uhI_ulI_pI_H"),           # compartments that will be summed and plotted
    groups = plotGroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                     # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    confIntv = confIntv,                       # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Active infections in ",titleString),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )
  
  trajPlotPhi <- simInfPlottingFunction(
    result = result,                           # model result
    table = "V",                               # which table: U or V
    compts= "phi",                             # compartments that will be plotted
    groups = plotgroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                       # Automatically sums over all groups
    uNames = NULL,                             # list of compartments
    vNames = names(v0),                        # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = "1 month",                    # plot parameter: Date axis format
    titleString = paste0("Interventions in ",titleString),   # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Intensity of intervention",          # plot parameter: Title of y axis
    lString = "Intervention metric",                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )

  if(hospRateExp > 0 | hospRatePost > 0){
    trajPlotHosp <- simInfPlottingFunction(
      result = result,                           # model result
      table = "U",                               # which table: U or V
      compts= c("H"),                            # compartments that will be plotted
      groups = plotGroups,                       # List of groups to aggregate and plot
      sumGroups = sumGroups,                       # Automatically sums over all groups
      uNames = names(u0),                        # list of compartments
      vNames = NULL,                             # list of variables
      rollM = 7,                             # number of days for rolling mean
      allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
      plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
      percentile = percentile,                   # percentile for central (or other) tendency
      confIntv = confIntv,              # confidence interval for plotting spread
      nTM = nodeTrialMat,                        # node-Trial matrix
      tS = tspan,                                # length of simulation
      enn = N,                                   # number of nodes per trial
      nT = numTrials,                            # number of trials in simulation
      startDate = startofSimDay,                 # start date of simulation
      dateBreaks = dateBreaks,          # plot parameter: Date axis format
      titleString = paste0("Hospitalizations, 7-day average, ",titleString),          # plot parameter: Title of plot
      xString = "Date",                          # plot parameter: Title of x axis
      yString = "Number of hospital beds",          # plot parameter: Title of y axis
      lString = lString,                         # plot parameter: Title of legend
      cString = cString                          # plot parameter: Plot caption
    )
  }

  trajPlotnewI <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("newI"),                         # compartments that will be plotted
    groups = plotGroups,                       # List of groups to aggregate and plot
    sumGroups = sumGroups,                          # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
    plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = percentile,                   # percentile for central (or other) tendency
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Daily new infections in ",titleString),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = lString,                         # plot parameter: Title of legend
    cString = cString                          # plot parameter: Plot caption
  )

  if(hospDeathRate+nonHospDeathRate > 0){
    trajPlotDeaths <- simInfPlottingFunction(
      result = result,                           # model result
      table = "U",                               # which table: U or V
      compts= c("M"),                            # compartments that will be plotted
      groups = plotGroups,                       # List of groups to aggregate and plot
      sumGroups = sumGroups,                       # Automatically sums over all groups
      uNames = names(u0),                        # list of compartments
      vNames = NULL,                             # list of variables
      rollM = rollM,                             # number of days for rolling mean
      allTraj = allTraj,                         # Logical for plotting all simulation trajectories or just median and spread
      plotRandomTrajs = plotRandomTrajs,         # If allTraj = true, can plot random trajectories, select the number desired
      percentile = percentile,                   # percentile for central (or other) tendency
      confIntv = confIntv,              # confidence interval for plotting spread
      nTM = nodeTrialMat,                        # node-Trial matrix
      tS = tspan,                                # length of simulation
      enn = N,                                   # number of nodes per trial
      nT = numTrials,                            # number of trials in simulation
      startDate = startofSimDay,                 # start date of simulation
      dateBreaks = dateBreaks,          # plot parameter: Date axis format
      titleString = paste0("Cumulative deaths in ",titleString),        # plot parameter: Title of plot
      xString = "Date",                          # plot parameter: Title of x axis
      yString = "Cumulative deaths",             # plot parameter: Title of y axis
      lString = lString,                         # plot parameter: Title of legend
      cString = cString                          # plot parameter: Plot caption
    )
  }
  
  ggsave(paste0(plotName,"I.png"),trajPlotInfections,width=9,height=5,units="in")
  ggsave(paste0(plotName,"phi.png"),trajPlotPhi,width=9,height=5,units="in")
  ggsave(paste0(plotName,"newI.png"),trajPlotnewI,width=9,height=5,units="in")
  if(hospRateExp > 0 | hospRatePost > 0){ ggsave(paste0(plotName,"hosp.png"),trajPlotHosp,width=9,height=5,units="in") }
  if(hospDeathRate + nonHospDeathRate > 0){ ggsave(paste0(plotName,"deaths.png"),trajPlotDeaths,width=9,height=5,units="in") }
  
  return(result)
}
