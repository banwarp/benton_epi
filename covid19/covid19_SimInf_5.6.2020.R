# covid19_SimInf_5.6.2020.R

# copyright Peter Banwarth Benton County Health Department 2020

# changes from covid19_SimInf_5.5.200.R
# added isolation = Is compartment that exposed can enter to represent enhanced contact tracing

# changes from covid19_SimInf_5.1.2020.R
# reorganized parameters
#   moved constant parameters from V0 to ldata if they varied across trial
#   moved constant parameters from v0/ldata to gdata/pts_function parameters if they didn't vary across trial
#   requires significant changes in pts_fun
# added in pdDecay - decay parameter to represent the slow decay of physical distancing as people return to normal interactions
# added confidence interval parameter for plotting different spreads

# changes from covid19_SimInf_4.30.2020.R
# changed R0 and gamma to match IDM 4/22 data
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

# Using package SimInf to model stochastic disease spread
# https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
# https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
# https://github.com/stewid/SimInf

# Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
# Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

# Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
# International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723


library(ggplot2)
library(reshape2)
library(SimInf)
library(data.table)
library(dplyr)
library(Matrix)
library(zoo)

setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("studentEventFunction.R")
source("simInfPlottingFunction5.5.2020.R")
source("pts_funScript5.5.2020.R")

setwd("../")

#############################################
#############################################
################ parameters #################
#############################################
#############################################

parmList = list(
  ### Simulation parameters
  "simID" = "covid5.6.2020.1",            # Simulation ID
  "simDate" = "2020-04-01",                  # Start date of simulation
  "maxT" = 365,                              # Length of simulation in days
  "numTrials" = 100,                         # Number of trials in simulation
  
  ### population parameters
  "trialPop" = 65000,                        # Total population in each trial
  "N" = 500,                                 # Number of population nodes
  "I0Pop" = 40,                              # Initial number of infectious
  "maxINodeProp" = 1/10,                     # Maximum proportion of nodes that intially have one or more infectious
  "E0Pop" = 0,                               # Initial number of exposed
  "Is0Pop" = 0,                              # Initial number of isolated
  "R0Pop" = 0,                               # Initial number of recovered
  "Im0Pop" = 0,                              # Initial number of immune
  "M0Pop" = 0,                               # Initial number of dead
  
  ### gdata parameters
  # disease dynamics
  "R0" = 2.9,                                # Basic reproduction number
  "gamma" = 1/8,                             # Reciprocal of infectious period
  "sigma" = 1/4,                             # Reciprocal of exposed period
  "rho" = .25,                               # Proportion of exposed who are identified and isolated before they become infectious
  "lambda" = 1/10,                           # Reciprocal of length of isolation period
  "RIsolated" = .125,                        # Reproduction number of isolated infectious
  "eta" = .024,                              # Case fatality rate
  "delta" = .1,                              # Proportion of recovereds who eventually become susceptible again
  "kappa" = 1/100,                           # Reciprocal of temporary immunity period, after which R becomes Im or S
  "mu" = 0,                                  # Natural birth/susceptible immigration rate
  "nu" = 0,                                  # Natural non-Covid death rate
  ### ldata parameters
  "R0Spread" = .1,                           # Uniform variation in R0 across trials (measured as %change from R0)
  ### continuous initialized parameters (v0)
  "phi0" = 2.9/.9,                           # Initial phi under current (March 23rd) stay-at-home orders
  
  ### pts_fun parameters
  "RPhysicalDistancing" = 2,                 # Ongoing baseline Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  "RNoAction" = 2.9,                         # Ongoing baseline Rt if no actions at are all taken
  "maxPrev1" = 15,                           # Maximum prevalence before instituting minor intervention
  "maxPrev2" = 30,                           # Maximum prevalence before instituting major intervention
  "RTarget1" = 1.5,                          # Target for the reduction in R under minor intervention
  "RTarget2" = .9,                           # Target for the reduction in R under major intervention
  "cosAmp" = .25,                            # Amplitude of seasonal variation in beta
  "upDelay" = 10,                            # Number of days after prevalence passes threshold until minor/major intervention
  "downDelay" = 28,                          # Number of days after prevalence drops below threshold until intervention lifted
  "phiMoveUp" = .25,                         # Rate at which phi increases when interventions are imposed
  "phiMoveDown" = .25,                       # Rate at which phi decreases when interventions are lifted
  "pdDecay" = 30,                            # Number of days until phi decays toward 1 in the absence of interventions.
                                             # Represents gradual relaxation of physical distancing as people return to normal.
                                             # pdDecay = -1 removes this decay factor
 
   ### other date parameters
  "kbDay1" = "2020-05-25",                   # Date of first phase of lifting stay-at-home orders
  "kbDay2" = "2020-06-14",                   # Date of second phase of lifting stay-at-home orders
  "kbDay3" = "2020-07-01",                   # Date of third phase of lifting stay-at-home orders
  "switchOffPolicies" = 0,                   # Indicator if intervention policies will cease after a certain day
  "switchOffDay" = "2021-03-01",             # Date intervention policies will cease if indicator == 1
  
  ### event parameters
  "paraChi_df" = 4,                          # Distribution parameter for parachute events
  "parachuteRate" = 1/30,                    # Reciprocal of expected waiting time for a parachute event
  "parachuteNum" = 1,                        # Number of infectious in each parachute event
  "transferRate" = 1/7,                      # Reciprocal of expected waiting time for transfer event
  "transferMinProp" = .01,                   # Minimum proportion of node population that transfers in each event
  "transferMaxProp" = .1,                    # Maximum proportion of node population that transfers in each event
  
  ### student parameters
  "studentReturnDate" = "2020-09-21",        # Date OSU student return
  "studentReturnSpread" = 3,                 # Symmetric spread of days to spread out student return
  "studentPop" = 25000,                      # OSU student population
  "maxStudentNodes" = 100,                   # Maximum number of nodes that students enter
  "sSProp" = .9,                             # Proportion of students who are susceptible
  "sEProp" = .0001,                          # Proportion of students who are exposed
  "sIProp" = .001,                           # Proportion of students who are infectious; R, Im, and M are calculated
  
  ### plot parameters
  "plotCompList" = "I",                      # List of compartments that will be plotted
  "rollM" = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  "confIntv" = .95,                          # Confidence interval for plotting spread
  "dateBreaks" = "1 month",                  # Plot parameter, x-axis displays months
  "titleString" = "Daily active infections in Benton County", # Title of plot
  "xString" = "Date",                        # Title of x-axis
  "yString" = "Frequency",                   # Title of y-axis
  "lString" = "Compartment"                  # Title of legend
)

# saving parameter list
setwd("./parameters")
sink(paste0("parms",parmList$simID,".txt"))
parmList
sink()
setwd("../")

# Start of simulation
startofSimDay <- as.numeric(as.Date(parmList$simDate)) - as.numeric(as.Date("2020-01-01"))
# Date current physical distancing is relaxed
kbDay1 <- as.numeric(as.Date(parmList$kbDay1)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
kbDay2 <- as.numeric(as.Date(parmList$kbDay2)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
kbDay3 <- as.numeric(as.Date(parmList$kbDay3)) - as.numeric(as.Date("2020-01-01")) - startofSimDay

# Time span
tspan <- 1:parmList$maxT

# plotting parameters
plotCompList <- parmList$plotCompList   # plots these compartments
rollM <- parmList$rollM                   # number of days for rolling mean
dateBreaks <- parmList$dateBreaks      # date axis formatting
titleString <- parmList$titleString  # plot title
xString <-  parmList$xString           # x axis title
yString <- parmList$yString       # y axis title
fileName <- paste0("data",parmList$simID,".txt")
plotName <- paste0("plot",parmList$simID)

# Trial and node parameters
numTrials <- min(100,parmList$numTrials) # capping number of trials at 100
N <- min(1000,parmList$N) # capping number of nodes per trial at 1000
NnumTrials <- N*numTrials
nodeTrialMat <- data.table(node=c(1:(NnumTrials)),trial=rep(1:numTrials,each=N))
NList <- split(nodeTrialMat$node,nodeTrialMat$trial)
trialPop <- parmList$trialPop # total population of a trial
maxNodePop <- floor(trialPop/N) # maximum node population
I0Pop <- parmList$I0Pop  # population that is infectious in intial state
maxINodeProp <- parmList$maxINodeProp # proportion of nodes with one or more infectious in initial state

# parachute and transfer parameters

# random chisquare distribution for selecting parachute events; truncate at no more than tspan;
# assumes parachuting starts slow (because of travel restrictions), grows, then drops off as global prevalence decreases
# decrease paraChi_df to left-skew; increase df to right-skew
parachuteDist <- ceiling(rchisq(10000,df=parmList$paraChi_df)*max(tspan)/10)
parachuteDist <- parachuteDist[parachuteDist < max(tspan)]

# Model parameters
# compartments
compartments <- c("S","E","I","Is","R","Im","M")

# student event parameters - relies on compartments
# September 20th (day students return) = day 264 in calendar year
studentReturnDate <- as.numeric(as.Date(parmList$studentReturnDate)) - as.numeric(as.Date("2020-01-01")) - startofSimDay
studentReturnSpread <- parmList$studentReturnSpread # symmetric spread around return date in number of days
studentPop <- parmList$studentPop
maxStudentNodes <- min(ceiling(N/3),parmList$maxStudentNodes) # maximum number of nodes students can enter
sSProp <- parmList$sSProp # proportion of students who are susceptible
sEProp <- parmList$sEProp # proportion of students who are exposed
sIProp <- parmList$sIProp # proportion of students who are infectious
sRProp <- (1-(sSProp+sEProp+sIProp))*(1-parmList$delta) # proportion of students who are recovered
sImProp <- (1-(sSProp+sEProp+sIProp))*parmList$delta # proportion of students who are immune
sIsProp <- 0
sMProp <- 0
studentPropTable <- data.frame(compartment = compartments, frac = c(sSProp,sEProp,sIProp,sIsProp,sRProp,sImProp,sMProp))

# More model parameters
# global parameters
gdata = data.frame(
  beta = parmList$R0*parmList$gamma,    # Transmission rate = R0/gamma
  gamma = parmList$gamma,               # Reciprocal of infectious period
  sigma = parmList$sigma,               # Reciprocal of exposed period
  rho = parmList$rho,
  lambda = parmList$lambda,
  betaIsolated = parmList$RIsolated*parmList$lambda,
  eta = parmList$eta,                   # Case fatality rate
  delta = parmList$delta,               # Proportion of recovered who become re-susceptible
  kappa = parmList$kappa,               # Reciprocal of recovered period before re-susceptibility
  mu = parmList$mu,                     # Natural birth/immigration rate
  nu = parmList$nu                      # Non-Covid death rate
)

# Transitions
transitions <- c(
  "@ -> mu*(S+E+I+R+Im) -> S",                                        # natural birth rate
  "R -> delta*kappa*R -> S",                                          # Recovereds who become susceptible again
  "S -> nu*S -> @",                                                   # Non-Covid deaths among susceptibles
  "S -> ((1/phi)*beta*betaRandomizer*(season)*I+betaIsolated*Is)*S/(S+E+I+R+Im)-> E",      # Exposures, allows for effective or ineffective policies and seasonal variation in beta
  "E -> nu*E -> @",                                                   # Non-Covid deaths among exposed
  "E -> (1-rho)*sigma*E -> I",                                        # Development of disease among non-isolated exposed
  "E -> rho*sigma*E -> Is",                                           # Isolation of exposed
  "I -> nu*I -> @",                                                   # Non-Covid deaths among infectious
  "I -> (1-eta)*gamma*I -> R",                                        # Recovery among infectious
  "I -> eta*gamma*I -> M",                                            # Death from COVID
  "Is -> nu*Is -> @",                                                 # Non-Covid deaths among isolated
  "Is -> (1-eta)*lambda*Is -> R",                                     # Recovery from COVID among isolated
  "Is -> eta*lambda*Is -> M",                                         # Death from COVID among isolated
  "R -> nu*R -> @",                                                   # Non-Covid deaths among recovered
  "R -> (1-delta)*kappa*R -> Im",                                     # Development of immunity
  "Im -> nu*Im -> @"                                                  # Non-Covid deaths among immune
)

# Initial state
# Randomly distribute population between nodes with bounds
S0 <- matrix(floor(maxNodePop/log(maxNodePop)*log(runif(NnumTrials,min=1,max=maxNodePop))),ncol=numTrials)
S0 <- matrix(apply(S0,2,function(x) x+round(((trialPop-sum(x))/length(x)),0)))
I0indices <- lapply(0:(numTrials-1),function(x) (x*N)+sample(N,rbinom(1,2*I0Pop,.5),replace=TRUE))
I0indices <- unlist(lapply(I0indices,table))
I0 <- matrix(0,nrow=NnumTrials)
I0[as.numeric(names(I0indices))] <- I0indices

# Initial state matrix
u0 <- data.frame(
  S = S0,
  E = rep(0,NnumTrials),
  I = I0, # random number of initial infectious, allows for 0, capped at maxIPop for each trial
  Is = rep(0,NnumTrials),
  R = rep(0,NnumTrials),
  Im = rep(0,NnumTrials),
  M=rep(0,NnumTrials)
)

# local parameters
ldata <- data.frame(
  firstNode = sort(rep.int(0:(numTrials-1)*N,N),method="quick"), # first node in a trial, for computing trial prevalence
  betaRandomizer = rep(runif(numTrials,min=1-parmList$R0Spread,max=1+parmList$R0Spread),each=N)  # randomizer for beta or R0
)

# Initialized continuous variables
v0 = data.frame(
  phi = rep(parmList$phi0,NnumTrials),                # initial beta reduction factor (larger phi = more reduction)
  prevUp01 = rep(parmList$upDelay+1,NnumTrials),      # initializedvariable for capturing delay in response when prevalence increases
  prevUp12 = rep(parmList$upDelay+1,NnumTrials),      # variable for capturing delay in response when prevalence increases
  prevDown21 = rep(parmList$downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  prevDown10 = rep(parmList$downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  kbPhase = rep(0,NnumTrials),                        # logical: 1 means current physical distancing, 0 means current physical distancing completely lifted
  policy = rep(1,NnumTrials),                         # logical: 1 means policies can take effect; 0 means they don't
  season = rep(1,NnumTrials),                         # seasonality factor for beta
  prevalence = rep(0,NnumTrials),                     # prevalence tracker
  previousState = rep(2,NnumTrials),                  # previous state: 0 = baseline, 1 = low intervention, 2 = high intervention, used for low intervention logic
  RTee = rep(parmList$R0,NnumTrials),                 # tracking effective Rt
  pdCounter = rep(0,NnumTrials)                  # counter for pdDecay
)

# building pts_fun
pts_fun <- pts_funScript(
  phiPhysicalDistancing = (gdata$beta/gdata$gamma)/parmList$RPhysicalDistancing, # Ongoing Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  phiNoAction = (gdata$beta/gdata$gamma)/parmList$RNoAction, # Ongoing Rt, reflecting no actions
  maxPrev1 = parmList$maxPrev1,                              # Maximum prevalence before instituting minor intervention
  maxPrev2 = parmList$maxPrev2,                              # Maximum prevalence before instituting major intervention
  phiFactor1 = (gdata$beta/gdata$gamma)/parmList$RTarget1,   # Target for the reduction in R under minor intervention
  phiFactor2 = (gdata$beta/gdata$gamma)/parmList$RTarget2,   # Target for the reduction in R under major intervention
  cosAmp = parmList$cosAmp,                                  # Amplitude of seasonal variation in beta
  startDay = startofSimDay,                                  # start of simulation
  kbDay1 = as.numeric(as.Date(parmList$kbDay1)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of first phase
  kbDay2 = as.numeric(as.Date(parmList$kbDay2)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of second phase
  kbDay3 = as.numeric(as.Date(parmList$kbDay3)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of third phase
  enn = N,                                                   # number of nodes in a trial
  numComp = length(compartments),                            # number of compartments
  prevType = ifelse(parmList$maxPrev1 < 1,1,0),              # prevalence type. 0 = count, 1 = proportion
  upDelay = parmList$upDelay,                                # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = parmList$downDelay,                            # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = parmList$phiMoveUp,                            # Rate at which phi increases when interventions are imposed
  phiMoveDown = parmList$phiMoveDown,                        # Rate at which phi decreases when interventions are lifted
  pdDecay = parmList$pdDecay,                         # Rate at which phi decreases toward 1 in the absence of interventions. Represents gradual relaxation of physical distancing
  switchOffPolicies = parmList$switchOffPolicies,            # Logical variable to switch off policies (used for counterfactuals)
  switchOffDay = as.numeric(as.Date(parmList$switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
)



############### EVENTS ##################
############### EVENTS ##################
############### EVENTS ##################
############### EVENTS ##################
############### EVENTS ##################

# create events: occasional infections parachute into the nodes; occasional transfers between the nodes. Students arrive in a wave

# Events matrix
# leaves of M from column 7
E <- cbind(diag(length(compartments)),c(rep(1,length(compartments)-2),0,0))
dimnames(E) <- list(compartments,c(1:ncol(E)))

###### PARACHUTE EVENTS ######
# Follows a poisson process
# number of events in each trial
numParachuteList <- lapply(1:numTrials,function(x) rpois(1,parmList$parachuteRate*max(tspan)))

if(max(unlist(numParachuteList)) > 0){
  # events data frame
  parachuteEvents <- data.frame(
    event = "enter",
    time = unlist(lapply(numParachuteList, function(x) sample(parachuteDist,x,replace=TRUE))),
    node = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numParachuteList[[x]],replace=TRUE))),
    dest = 0,
    n = parmList$parachuteNum,
    proportion = 0,
    select = 2,
    shift = 0
  )
}
###### TRANSFER EVENTS ######
# Follows a poisson process
# Partition transfers by trials
numTransferList <- lapply(1:numTrials,function(x) rpois(1,parmList$transferRate*max(tspan)))

if(max(unlist(numTransferList))>0) {
  # Transfer data frame
  # first non-students
  transferEvents <- data.frame(
    event = "extTrans",
    time = unlist(lapply(numTransferList, function(x) sample(tspan,x,replace=TRUE))),
    node = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numTransferList[[x]],replace=TRUE))),
    dest = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numTransferList[[x]],replace=TRUE))),
    n = 0,
    proportion = unlist(lapply(numTransferList,function(x) round(runif(x,min=parmList$transferMinProp,max=parmList$transferMaxProp),3))),
    select = ncol(E),
    shift = 0
  )
  
  # Capping transfer events proportion so large populations don't overwhelm small populations
  # ceiling on transfer populations so large populations don't overwhelm small populations
  Sgrid <- data.frame(node=1:length(u0$S),nodePop = u0$S)
  transferEvents <- merge(transferEvents,Sgrid,by="node",all.x=TRUE)
  names(Sgrid) <- c("dest","destPop")
  transferEvents <- merge(transferEvents,Sgrid,by="dest",all.x=TRUE)
  transferEvents$maxProp <- transferEvents$destPop/transferEvents$nodePop
  transferEvents[transferEvents$maxProp<1,"proportion"] <- transferEvents[transferEvents$maxProp<1,"maxProp"]*transferEvents[transferEvents$maxProp<1,"proportion"]
  transferEvents <- transferEvents[,c("event","time","node","dest","n","proportion","select","shift")]
}

###### STUDENT ARRIVAL EVENTS ######

if(maxStudentNodes > 0) {
  # Event for day(s) students return, with random returns around studentReturnDate
  studentNodes <- sort(sample(1:N,min(N,maxStudentNodes)))
  
  # exclude deceased by length(compartments)-1
  studentEventList <- lapply(c(1:(length(compartments)-1)),studentEventFunction)
  studentEvents <- bind_rows(studentEventList[lengths(studentEventList) != 0])
}

###### ALL EVENTS ######
allEvents <- data.frame(NULL)
if(max(unlist(numParachuteList)) > 0) { allEvents <- rbind(allEvents,parachuteEvents) }
if(max(unlist(numTransferList))>0) { allEvents <- rbind(allEvents,transferEvents) }
if(maxStudentNodes > 0) { allEvents <- rbind(allEvents,studentEvents) }

###### BUILDING MODEL ######
###### BUILDING MODEL ######
###### BUILDING MODEL ######
###### BUILDING MODEL ######
###### BUILDING MODEL ######

if(nrow(allEvents)>0) {
  model <- mparse(transitions = transitions, compartments = compartments, events = allEvents, E=E,
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

set.seed(123)
result <- run(model)

###### PLOTTING RESULT ######
###### PLOTTING RESULT ######
###### PLOTTING RESULT ######
###### PLOTTING RESULT ######
###### PLOTTING RESULT ######

setwd("./trajectories")

# trajTable <- data.table(trajectory(result))
# fwrite(trajTable,file=fileName,sep=",")

# drawing plot of trajectory
trajPlotI <- simInfPlottingFunction(
  result = result,                           # model result
  table = "U",                               # which table: U or V
  compts= c("I"),                            # compartments that will be plotted
  uNames = names(u0),                        # list of compartments
  vNames = NULL,                             # list of variables
  rollM = rollM,                             # number of days for rolling mean
  confIntv = parmList$confIntv,              # confidence interval for plotting spread
  nTM = nodeTrialMat,                        # node-Trial matrix
  tS = tspan,                                # length of simulation
  enn = N,                                   # number of nodes per trial
  nT = numTrials,                            # number of trials in simulation
  startDate = startofSimDay,                 # start date of simulation
  dateBreaks = parmList$dateBreaks,          # plot parameter: Date axis format
  titleString = parmList$titleString,        # plot parameter: Title of plot
  xString = "Date",                          # plot parameter: Title of x axis
  yString = "Number of infections",          # plot parameter: Title of y axis
  lString = "Compartment"                    # plot parameter: Title of legend
)

trajPlotPhi <- simInfPlottingFunction(
  result = result,                           # model result
  table = "V",                               # which table: U or V
  compts= "phi",                             # compartments that will be plotted
  uNames = NULL,                             # list of compartments
  vNames = names(v0),                        # list of variables
  rollM = rollM,                             # number of days for rolling mean
  nTM = nodeTrialMat,                        # node-Trial matrix
  tS = tspan,                                # length of simulation
  enn = N,                                   # number of nodes per trial
  nT = numTrials,                            # number of trials in simulation
  startDate = startofSimDay,                 # start date of simulation
  dateBreaks = "1 month",                    # plot parameter: Date axis format
  titleString = "Physical distancing interventions",   # plot parameter: Title of plot
  xString = "Date",                          # plot parameter: Title of x axis
  yString = "Intensity of intervention",          # plot parameter: Title of y axis
  lString = "Intervention metric"                    # plot parameter: Title of legend
)

# trajPlotI <- trajPlotI+
#     geom_vline(xintercept=as.numeric(as.Date("2020-07-07")),color="darkseagreen",size=1)+
#     geom_text(aes(x=as.Date("2020-07-07"), label="Lift all intv.", y=60), angle=90, vjust = -1, color = "black")
# 
# trajPlotPhi <- trajPlotPhi+
# geom_vline(xintercept=as.numeric(as.Date("2020-07-07")),color="darkseagreen",size=1.5)
# 
# ggsave(paste0(plotName,"I.png"),trajPlotI,width=9,height=5,units="in")
# ggsave(paste0(plotName,"phi.png"),trajPlotPhi,width=9,height=5,units="in")
