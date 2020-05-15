# covid19_SimInf_5.14.2020.R

# copyright Peter Banwarth Benton County Health Department 2020

# Using package SimInf to model stochastic disease spread
# https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
# https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
# https://github.com/stewid/SimInf

# Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
# Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

# Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
# International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

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
# changed file name for studentEVentFunction script to massEventFunction, and added super spreader function as another function in that script
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
# added more student functionality
# added local data parameters for more variability

# changes from covid19_SimInf_4.9.2020.R
# removed all the efforts before Effort 15

covidWrapper <- function(
  ### Simulation parameters
  folderPath = NULL,                       # folder path for subroutines and output
  simID = "covid.mo.day.2020.Generic",     # Simulation ID
  simDate = "2020-05-01",                  # Start date of simulation
  maxT = 365,                              # Length of simulation in days
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  trialPop = 65000,                        # Total population in each trial
  N = 500,                                 # Number of population nodes
  nodeGroup = 1,                           # List of group IDs for node groups. Default is 1 group of nodes.
  I0Pop = round(.0023*65000,0),            # Initial number of infectious
  maxINodeProp = 1/10,                     # Maximum proportion of nodes that intially have one or more infectious
  E0Pop = 0,                               # Initial number of exposed
  R0Pop = 0,                               # Initial number of recovered
  Im0Pop = 0,                              # Initial number of immune
  H0Pop = 0,                               # Initial number of hospitalized
  Is0Pop = 0,                              # Initial number of isolated
  M0Pop = 0,                               # Initial number of dead
  
  ### gdata parameters
  R0 = 2.9,                                # Basic reproduction number
  infectiousPeriod = 1/8,                  # Reciprocal of infectious period
  exposedPeriod = 1/4,                     # Reciprocal of exposed period
  isoRate = .125,                          # Proportion of exposed who are identified and isolated before they become infectious
  isoPeriod = 1/10,                        # Reciprocal of length of isolation period
  RIsolated = .125,                        # Reproduction number of isolated infectious
  hospRate = .033,                         # Proportion of infectious/isolated that are hospitalized
  hospPeriod = 1/14,                       # Reciprocal of length of hospitalization period
  nonHospDeathRate = 0,                    # Non-hospitalized fatality rate
  hospDeathRate = .125,                    # Hospitalized fatality rate NOT CASE FATALITY RATE
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = 1/100,                   # Reciprocal of temporary immunity period, after which R becomes Im or S
  mu = 0,                                  # Natural birth/susceptible immigration rate
  nu = 0,                                  # Natural non-Covid death rate
  
  ### ldata parameters
  R0Spread = .1,                           # Uniform variation in R0 across trials (measured as %change from R0)
  
  ### continuous initialized parameters (v0)
  phi0 = 2.9/.9,                           # Initial phi under current (March 23rd) stay-at-home orders
  
  ### pts_fun parameters
  RPhysicalDistancing = 2,                 # Ongoing baseline Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  RNoAction = 2.9,                         # Ongoing baseline Rt if no actions at are all taken
  maxPrev1 = .001,                         # Maximum prevalence before instituting minor intervention
  maxPrev2 = .002,                         # Maximum prevalence before instituting major intervention
  RTarget1 = 1.5,                          # Target for the reduction in R under minor intervention
  RTarget2 = .9,                           # Target for the reduction in R under major intervention
  cosAmp = 0.25,                           # Amplitude of seasonal variation in beta
  upDelay = 10,                            # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = 28,                          # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = .25,                         # Rate at which phi increases when interventions are imposed
  phiMoveDown = .25,                       # Rate at which phi decreases when interventions are lifted
  pdDecay = 30,                            # Number of days until RPhysicalDistancing decays toward RNoAction in the absence of interventions.
  ########################################## Represents gradual relaxation of physical distancing as people return to normal.
  ########################################## pdDecay = -1 removes this decay factor
  
  ### other date parameters
  kbDay1 = "2020-05-25",                   # Date of first phase of lifting stay-at-home orders
  kbDay2 = "2020-06-20",                   # Date of second phase of lifting stay-at-home orders
  kbDay3 = "2020-07-15",                   # Date of third phase of lifting stay-at-home orders
  switchOffPolicies = 0,                   # Indicator if intervention policies will cease after a certain day
  switchOffDay = "2020-07-15",             # Date intervention policies will cease if indicator == 1
  
  ### event parameters
  paraChi_df = 4,                          # Distribution parameter for parachute events
  parachuteRate = 1/21,                    # Reciprocal of expected waiting time for a parachute event
  parachuteNum = 1,                        # Number of infectious in each parachute event
  inGroupTransferRate = 1/7,               # Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = .1,      # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = .1,             # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = .3,             # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = 1/15,             # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = .05,    # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = .05,           # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = .15,           # Maximum proportion of node population that transfers in each event
  
  ### super spreader parameters - list to allow multiple events
  superInfections = c(),                   # Number of infections caused by the super spreader
  superNodes = c(),                        # Number of nodes that the super spreader contacts
  superDate = c(),                         # Date the super spreader lands. Date can also be numeric i.e. 200
  superSpread = c(),                       # Symmetric spread in days of super spreader infections
  
  ### student parameters
  studentReturnDate = "2020-09-21",        # Date OSU student return
  studentReturnSpread = 3,                 # Symmetric spread of days to spread out student return
  studentPop = 25000,                      # OSU student population
  maxStudentNodes = 0,                     # Maximum number of nodes that students enter. Set to 0 to remove studentEvent
  sSProp = .9,                             # Proportion of students who are susceptible
  sEProp = .0001,                          # Proportion of students who are exposed
  sIProp = .0015,                          # Proportion of students who are infectious; R, Im, and M are calculated
  
  ### plot parameters
  plotCompList = "I",                      # List of compartments that will be plotted
  rollM = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  confIntv = .95,                          # Confidence interval for plotting spread
  dateBreaks = "1 month",                  # Plot parameter, x-axis displays months
  titleString = "Generic Title",           # Title of plot
  xString = "Date",                        # Title of x-axis
  yString = "Frequency",                   # Title of y-axis
  lString = "Compartment"                  # Title of legend
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
  library(Matrix)
  library(zoo)
  
  if(!is.null(folderPath)) {setwd(folderPath)}
  
  if(is.null(folderPath)) {setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")}
  source("eventFunctions5.14.2020.R")
  source("simInfPlottingFunction5.14.2020.R")
  source("pts_funScript5.5.2020.R")
  if(is.null(folderPath)) {setwd("../")}
  
  # saving parameters
  if(is.null(folderPath)) {setwd("./parameters")}
  sink(paste0("parms",simID,".txt"))
  print(paste("Simulation:",simID,sep=" "))
  print(paste("Population:",trialPop,sep=" "))
  print(paste("Initial infected:",I0Pop,sep=" "))
  print(paste("Number of nodes:",N,sep=" "))
  print(paste("Population per node:",ceiling(trialPop/N),sep=" "))
  print(paste("Proportion of exposed who are identified through contact tracing:",isoRate,sep=" "))
  if(switchOffPolicies==1) {print(paste("Policies switched off on",switchOffDay,sep=" "))}
  if(length(superNodes)>0) {
    print(paste("Super spreader introduced on",superDate,sep=" "))
    print(paste("Super spreader adds",superInfections,"infections to population",sep=" "))
    }
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
  
  # Trial and node parameters
  numTrials <- min(100,numTrials) # capping number of trials at 100
  N <- min(1000,N) # capping number of nodes per trial at 1000
  numGroups <- length(unique(nodeGroup))
  NnumTrials <- N*numTrials
  nodeTrialMat <- data.table(node=c(1:(NnumTrials)),
                             nodeGroup = rep(nodeGroup,numTrials),
                             trial=rep(1:numTrials,each=N))
  NList <- split(nodeTrialMat$node,nodeTrialMat$trial)
  nodeAndGroupList <- split(nodeTrialMat,nodeTrialMat$nodeGroup)
  nodeAndGroupList <- lapply(nodeAndGroupList,function(x) split(x$node,x$trial))
  maxNodePop <- floor(trialPop/N) # maximum node population
  
  # parachute and transfer parameters
  
  # random chisquare distribution for selecting parachute events; truncate at no more than tspan;
  # assumes parachuting starts slow (because of travel restrictions), grows, then drops off as global prevalence decreases
  # decrease paraChi_df to left-skew; increase df to right-skew
  parachuteDist <- ceiling(rchisq(10000,df=paraChi_df)*max(tspan)/10)
  parachuteDist <- parachuteDist[parachuteDist < max(tspan)]
  
  # Model parameters
  # compartments
  compartments <- c("S","E","I","R","Im","H","Is","cumI","M")
  
  # student event parameters - relies on compartments
  # September 20th (day students return) = day 264 in calendar year
  studentReturnDate <- as.numeric(as.Date(studentReturnDate)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
  studentReturnSpread <- studentReturnSpread # symmetric spread around return date in number of days
  maxStudentNodes <- min(ceiling(N/3),maxStudentNodes) # maximum number of nodes students can enter
  sSProp <- sSProp # proportion of students who are susceptible
  sEProp <- sEProp # proportion of students who are exposed
  sIProp <- sIProp # proportion of students who are infectious
  sRProp <- (1-(sSProp+sEProp+sIProp))*(1-reSuscepRate) # proportion of students who are recovered
  sImProp <- (1-(sSProp+sEProp+sIProp))*reSuscepRate # proportion of students who are immune
  sHProp <- 0
  sIsProp <- 0
  scumIProp <- sIProp
  sMProp <- 0
  studentPropTable <- data.frame(compartment = compartments,
                                 frac = c(sSProp,sEProp,sIProp,sRProp,sImProp,sHProp,sIsProp,scumIProp,sMProp))
  
  # More model parameters
  # global parameters
  gdata = data.frame(
    beta = R0*infectiousPeriod,    # Transmission rate = R0/infectiousPeriod
    infectiousPeriod = infectiousPeriod,               # Reciprocal of infectious period
    exposedPeriod = exposedPeriod,               # Reciprocal of exposed period
    isoRate = max(isoRate,hospRate),     # Acknowledges that all hospitalized will be isolated
    isoPeriod = isoPeriod,
    hospRate = hospRate,
    hospPeriod = hospPeriod,
    betaIsolated = RIsolated*isoPeriod,
    nonHospDeathRate = nonHospDeathRate, # non-hospitalized fatality rate
    hospDeathRate = hospDeathRate,                   # Hospitalized fatality rate NOT CASE FATALITY RATE
    reSuscepRate = reSuscepRate,               # Proportion of recovered who become re-susceptible
    tempImmPeriod = tempImmPeriod,               # Reciprocal of recovered period before re-susceptibility
    mu = mu,                     # Natural birth/immigration rate
    nu = nu                      # Non-Covid death rate
  )
  
  # Transitions
  transitions <- c(
    "@ -> mu*(S+E+I+R+Im) -> S",                                        # natural birth rate
    "R -> reSuscepRate*tempImmPeriod*R -> S",                           # Recovereds who become susceptible again
    "S -> nu*S -> @",                                                   # Non-Covid deaths among susceptibles
    "S -> ((1/phi)*beta*betaRandomizer*(season)*I+betaIsolated*Is)*S/(S+E+I+R+Im)-> E",      # Exposures, allows for effective or ineffective policies and seasonal variation in beta
    "E -> nu*E -> @",                                                   # Non-Covid deaths among exposed
    "E -> (1-isoRate)*exposedPeriod*E -> I + cumI",                     # Development of disease among non-isolated exposed
    "E -> (isoRate-hospRate)*exposedPeriod*E -> Is + cumI",             # Isolation of exposed (non-hospitalized)
    "E -> hospRate*exposedPeriod*E -> H",                               # New hospitalizations
    "I -> nu*I -> @",                                                   # Non-Covid deaths among infectious
    "I -> (1-nonHospDeathRate)*infectiousPeriod*I -> R",                # Recovery among infectious, assumes non-hosp infectious all recover
    "I -> nonHospDeathRate*infectiousPeriod*I -> M",                    # Deaths among non-hospitalized COVID patients
    "H -> (1-hospDeathRate)*hospPeriod*H -> R",                         # Recovery among hospitalized
    "H -> hospDeathRate*H -> M",                                        # Deaths among hospitalized COVID patients
    "Is -> nu*Is -> @",                                                 # Non-Covid deaths among isolated
    "Is -> isoPeriod*Is -> R",                                          # Recovery from COVID among isolated, assumes non-hosp isolated all recover
    "R -> nu*R -> @",                                                   # Non-Covid deaths among recovered
    "R -> (1-reSuscepRate)*tempImmPeriod*R -> Im",                      # Development of immunity
    "Im -> nu*Im -> @"                                                  # Non-Covid deaths among immune
  )
  
  # Initial state
  # Randomly distribute population between nodes with bounds
  S0 <- matrix(floor(maxNodePop/log(maxNodePop)*log(runif(NnumTrials,min=1,max=maxNodePop))),ncol=numTrials)
  S0 <- matrix(apply(S0,2,function(x) x+round(((trialPop-sum(x))/length(x)),0)))
  Inodes <- lapply(0:(numTrials-1),function(x) (x*N)+sample(N,ceiling(maxINodeProp*N)))
  I0indices <- lapply(Inodes,function(x) sample(x,rbinom(1,2*I0Pop,.5),replace=TRUE))
  I0indices <- unlist(lapply(I0indices,function(x) table(x)))
  I0 <- matrix(0,nrow=NnumTrials)
  I0[as.numeric(names(I0indices))] <- I0indices
  
  # Initial state matrix
  u0 <- data.frame(
    S = S0,
    E = rep(E0Pop,NnumTrials),
    I = I0, # random number of initial infectious, allows for 0, capped at maxIPop for each trial
    R = rep(R0Pop,NnumTrials),
    Im = rep(Im0Pop,NnumTrials),
    H = rep(H0Pop,NnumTrials),
    Is = rep(Is0Pop,NnumTrials),
    cumI = I0,
    M=rep(M0Pop,NnumTrials)
  )
  
  # local parameters
  ldata <- data.frame(
    firstNode = sort(rep.int(0:(numTrials-1)*N,N),method="quick"), # first node in a trial, for computing trial prevalence
    betaRandomizer = rep(runif(numTrials,min=1-R0Spread,max=1+R0Spread),each=N)  # randomizer for beta or R0
  )
  
  # Initialized continuous variables
  v0 = data.frame(
    phi = rep(phi0,NnumTrials),                # initial beta reduction factor (larger phi = more reduction)
    prevUp01 = rep(upDelay+1,NnumTrials),      # initializedvariable for capturing delay in response when prevalence increases
    prevUp12 = rep(upDelay+1,NnumTrials),      # variable for capturing delay in response when prevalence increases
    prevDown21 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
    prevDown10 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
    kbPhase = rep(0,NnumTrials),                        # logical: 1 means current physical distancing, 0 means current physical distancing completely lifted
    policy = rep(1,NnumTrials),                         # logical: 1 means policies can take effect; 0 means they don't
    season = rep(1,NnumTrials),                         # seasonality factor for beta
    prevalence = rep(0,NnumTrials),                     # prevalence tracker
    previousState = rep(2,NnumTrials),                  # previous state: 0 = baseline, 1 = low intervention, 2 = high intervention, used for low intervention logic
    RTee = rep(R0,NnumTrials),                 # tracking effective Rt
    pdCounter = rep(0,NnumTrials)                  # counter for pdDecay
  )
  
  # building pts_fun
  pts_fun <- pts_funScript(
    phiPhysicalDistancing = (gdata$beta/gdata$infectiousPeriod)/RPhysicalDistancing, # Ongoing Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
    phiNoAction = (gdata$beta/gdata$infectiousPeriod)/RNoAction, # Ongoing Rt, reflecting no actions
    maxPrev1 = maxPrev1,                              # Maximum prevalence before instituting minor intervention
    maxPrev2 = maxPrev2,                              # Maximum prevalence before instituting major intervention
    phiFactor1 = (gdata$beta/gdata$infectiousPeriod)/RTarget1,   # Target for the reduction in R under minor intervention
    phiFactor2 = (gdata$beta/gdata$infectiousPeriod)/RTarget2,   # Target for the reduction in R under major intervention
    cosAmp = cosAmp,                                  # Amplitude of seasonal variation in beta
    startDay = startofSimDay,                                  # start of simulation
    kbDay1 = kbDay1, # date of first phase
    kbDay2 = kbDay2, # date of second phase
    kbDay3 = kbDay3, # date of third phase
    enn = N,                                                   # number of nodes in a trial
    numComp = length(compartments),                            # number of compartments
    prevType = ifelse(maxPrev1 < 1,1,0),              # prevalence type. 0 = count, 1 = proportion
    upDelay = upDelay,                                # Number of days after prevalence passes threshold until minor/major intervention
    downDelay = downDelay,                            # Number of days after prevalence drops below threshold until intervention lifted
    phiMoveUp = phiMoveUp,                            # Rate at which phi increases when interventions are imposed
    phiMoveDown = phiMoveDown,                        # Rate at which phi decreases when interventions are lifted
    pdDecay = pdDecay,                         # Rate at which phi decreases toward 1 in the absence of interventions. Represents gradual relaxation of physical distancing
    switchOffPolicies = switchOffPolicies,            # Logical variable to switch off policies (used for counterfactuals)
    switchOffDay = as.numeric(as.Date(switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
  )
  
  
  
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  ############### EVENTS ##################
  
  # create events: occasional infections parachute into the nodes; occasional transfers between the nodes. Students arrive in a wave
  
  # Events matrix
  # leaves of M from column 7
  E <- cbind(diag(length(compartments)),c(rep(1,length(compartments)-4),rep(0,4)))
  dimnames(E) <- list(compartments,c(1:ncol(E)))
  
  ###### PARACHUTE EVENTS ######
  # Follows a poisson process
  # number of events in each trial
  numParachuteList <- lapply(1:numTrials,function(x) rpois(1,parachuteRate*max(tspan)))
  
  if(max(unlist(numParachuteList)) > 0){
    # events data frame
    parachuteEvents <- data.frame(
      event = "enter",
      time = unlist(lapply(numParachuteList, function(x) sample(parachuteDist,x,replace=TRUE))),
      node = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numParachuteList[[x]],replace=TRUE))),
      dest = 0,
      n = parachuteNum,
      proportion = 0,
      select = which(compartments == "I"),
      shift = 0
    )
  }
  
  ###### TRANSFER EVENTS ######
  # Making list of transfer parameters if they are the same across groups
  if(numGroups > 1 & length(inGroupTransferRate == 1))
  {
    inGroupTransferRate = rep(inGroupTransferRate,numGroups)
    inGroupTransferNodeNum = rep(inGroupNodeNum,numGroups)
    inGroupTransferMinProp = rep(inGroupMinProp,numGroups)
    inGroupTransferMaxProp = rep(inGroupMaxProp,numGroups)
    outGroupTransferRate = rep(outGroupTransferRate,numGroups)
    outGroupTransferNodeNum = rep(outGroupTransferNodeNum,numGroups)
    outGroupTransferMinProp = rep(outGroupTransferMinProp,numGroups)
    outGroupTransferMaxProp = rep(outGroupTransferMaxProp,numGroups)  
  }
  
  ###### IN-GROUP TRANSFER EVENTS ######
  # Follows a poisson process
  # Partition transfers by trials and groups
  if(inGroupTransferNodeNum[1]<1) {inGroupTransferNodeNum <- ceiling(inGroupTransferNodeNum*N)} # converts proportion to number
  
  numInTransferList <- lapply(1:numGroups,function(x) 
    lapply(1:numTrials, function(y) rpois(1,inGroupTransferRate[x]*inGroupTransferNodeNum[x]*max(tspan))))
  
  inTransferEventsList <- lapply(1:numGroups, transferFunction,
                               numTList = numInTransferList,
                               tS = tspan,
                               nT = numTrials,
                               nList = nodeAndGroupList,
                               minProp = inGroupTransferMinProp,
                               maxProp = inGroupTransferMaxProp,
                               Emat = E,
                               nTM = nodeTrialMat,
                               u0S = u0$S,
                               outLogic = FALSE
                               )
  inTransferEvents <- bind_rows(inTransferEventsList[lengths(inTransferEventsList) != 0])
   
  ###### OUT-GROUP TRANSFER EVENTS ######
  if(numGroups > 1) {
    if(outGroupTransferNodeNum[1]<1) {outGroupTransferNodeNum <- ceiling(outGroupTransferNodeNum*N)} # converts proportion to number
    
    numOutTransferList <- lapply(1:numGroups,function(x) 
      lapply(1:numTrials, function(y) rpois(1,outGroupTransferRate[x]*outGroupTransferNodeNum[x]*max(tspan))))
    
    outTransferEventsList <- lapply(1:numGroups, transferFunction,
                                    numTList = numOutTransferList,
                                    tS = tspan,
                                    nT = numTrials,
                                    nList = nodeAndGroupList,
                                    minProp = outGroupTransferMinProp,
                                    maxProp = outGroupTransferMaxProp,
                                    Emat = E,
                                    nTM = nodeTrialMat,
                                    u0S = u0$S,
                                    outLogic = TRUE
                                    )
    outTransferEvents <- bind_rows(outTransferEventsList[lengths(outTransferEventsList) != 0])
    transferEvents <- rbind(inTransferEvents,outTransferEvents)
    rm(outTransferEventsList)
    rm(outTransferEvents)
  } else {transferEvents <- inTransferEvents}
  rm(inTransferEventsList)
  rm(inTransferEvents)
  
  # Capping transfer events proportion so large populations don't overwhelm small populations
  # ceiling on transfer populations so large populations don't overwhelm small populations
  if(nrow(transferEvents)>0){
    Sgrid <- data.frame(node=1:length(u0$S),nodePop = u0$S)
    transferEvents <- merge(transferEvents,Sgrid,by="node",all.x=TRUE)
    names(Sgrid) <- c("dest","destPop")
    transferEvents <- merge(transferEvents,Sgrid,by="dest",all.x=TRUE)
    transferEvents$maxProp <- transferEvents$destPop/transferEvents$nodePop
    transferEvents[transferEvents$maxProp<1,"proportion"] <- transferEvents[transferEvents$maxProp<1,"maxProp"]*transferEvents[transferEvents$maxProp<1,"proportion"]
    transferEvents <- transferEvents[,c("event","time","node","dest","n","proportion","select","shift")]
  }
  
  ###### SUPER SPREADER EVENTS #####
  # shift matrix
  Nmat <- matrix(c(1,rep(0,length(compartments)-1)),ncol=1)
  dimnames(Nmat) <- list(compartments,1)
  
  if(length(superInfections) > 0){
    superNodes <- lapply(superNodes, function(x) sort(sample(1:N,min(N,x))))
    superEventList <- lapply(1:length(superInfections), superFunction,
                             IList=superInfections,
                             nodeList=superNodes,
                             dayList=superDate,
                             spreadList=superSpread,
                             startDay=startofSimDay,
                             nT = numTrials,
                             enn = N)
    superEvents <- bind_rows(superEventList[lengths(superEventList) != 0])
  }
  ###### STUDENT ARRIVAL EVENTS ######
  
  if(maxStudentNodes > 0) {
    # Event for day(s) students return, with random returns around studentReturnDate
    studentNodes <- sort(sample(1:N,min(N,maxStudentNodes)))
    
    # exclude isolated, cumulative, and deceased by length(compartments)-3
    studentEventList <- lapply(c(1:(length(compartments)-4)),studentEventFunction,
                               fracTable = studentPropTable,
                               nodeList = studentNodes,
                               pop = studentPop,
                               nT = numTrials,
                               enn = N,
                               sRD = studentReturnDate,
                               sRS = studentReturnSpread)
    studentEvents <- bind_rows(studentEventList[lengths(studentEventList) != 0])
    studentCumI <- studentEvents[which(studentEvents$select == which(compartments=="I")),]
    studentCumI$select <- which(compartments=="cumI")
    studentEvents <- rbind(studentEvents,studentCumI)
  }
  
  ###### ALL EVENTS ######
  allEvents <- data.frame(NULL)
  if(max(unlist(numParachuteList)) > 0) { allEvents <- rbind(allEvents,parachuteEvents) }
  if(nrow(transferEvents)>0) { allEvents <- rbind(allEvents,transferEvents) }
  if(maxStudentNodes > 0) { allEvents <- rbind(allEvents,studentEvents) }
  if(length(superInfections)>0) { allEvents <- rbind(allEvents,superEvents)}
  
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
  
  # trajTable <- data.table(trajectory(result))
  # fwrite(trajTable,file=fileName,sep=",")
  
  # drawing plot of infectious +  isolated to get total infected
  trajPlotInfections <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("I_Is_H"),                        # compartments that will be plotted
    groups = NULL,                             # List of groups to aggregate and plot
    sumGroups = TRUE,                          # Automatically sums over all groups
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Active infections in ",titleString," County"),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = "Compartment"                    # plot parameter: Title of legend
  )
  
  trajPlotPhi <- simInfPlottingFunction(
    result = result,                           # model result
    table = "V",                               # which table: U or V
    compts= "phi",                             # compartments that will be plotted
    groups = NULL,                           # List of groups to aggregate and plot
    sumGroups = TRUE,                          # Automatically sums over all groupsuNames = NULL,                             # list of compartments
    vNames = names(v0),                        # list of variables
    rollM = rollM,                             # number of days for rolling mean
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = "1 month",                    # plot parameter: Date axis format
    titleString = paste0("Interventions in ",titleString," County"),   # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Intensity of intervention",          # plot parameter: Title of y axis
    lString = "Intervention metric"                    # plot parameter: Title of legend
  )

  trajPlotHosp <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("H"),                            # compartments that will be plotted
    groups = NULL,                             # List of groups to aggregate and plot
    sumGroups = TRUE,                          # Automatically sums over all groupsuNames = names(u0),                        # list of compartments
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = 7,                             # number of days for rolling mean
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Hospitalizations, 7-day average, ",titleString," County"),          # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of hospital beds",          # plot parameter: Title of y axis
    lString = "Compartment"                    # plot parameter: Title of legend
  )

  trajPlotnewI <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("newI"),                            # compartments that will be plotted
    groups = NULL,                             # List of groups to aggregate and plot
    sumGroups = TRUE,                          # Automatically sums over all groupsuNames = names(u0),                        # list of compartments
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Daily new infections in ",titleString," County"),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Number of infections",          # plot parameter: Title of y axis
    lString = "Compartment"                    # plot parameter: Title of legend
  )

  trajPlotDeaths <- simInfPlottingFunction(
    result = result,                           # model result
    table = "U",                               # which table: U or V
    compts= c("M"),                            # compartments that will be plotted
    groups = NULL,                             # List of groups to aggregate and plot
    sumGroups = TRUE,                          # Automatically sums over all groupsuNames = names(u0),                        # list of compartments
    uNames = names(u0),                        # list of compartments
    vNames = NULL,                             # list of variables
    rollM = rollM,                             # number of days for rolling mean
    confIntv = confIntv,              # confidence interval for plotting spread
    nTM = nodeTrialMat,                        # node-Trial matrix
    tS = tspan,                                # length of simulation
    enn = N,                                   # number of nodes per trial
    nT = numTrials,                            # number of trials in simulation
    startDate = startofSimDay,                 # start date of simulation
    dateBreaks = dateBreaks,          # plot parameter: Date axis format
    titleString = paste0("Cumulative deaths in ",titleString," County"),        # plot parameter: Title of plot
    xString = "Date",                          # plot parameter: Title of x axis
    yString = "Cumulative deaths",             # plot parameter: Title of y axis
    lString = "Compartment"                    # plot parameter: Title of legend
  )

  ggsave(paste0(plotName,"I.png"),trajPlotInfections,width=9,height=5,units="in")
  ggsave(paste0(plotName,"phi.png"),trajPlotPhi,width=9,height=5,units="in")
  ggsave(paste0(plotName,"hosp.png"),trajPlotHosp,width=9,height=5,units="in")
  ggsave(paste0(plotName,"newI.png"),trajPlotnewI,width=9,height=5,units="in")
  ggsave(paste0(plotName,"deaths.png"),trajPlotDeaths,width=9,height=5,units="in")
  
  return(result)
}