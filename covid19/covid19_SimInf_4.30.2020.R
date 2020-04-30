# covid19_SimInf_4.30.2020.R

# copyright Peter Banwarth Benton County Health Department 2020

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
source("ldataFunction.R")
source("studentEventFunction.R")
source("simInfPlottingFunction.R")
source("pts_funScript.R")

setwd("../")

#############################################
#############################################
################ parameters #################
#############################################
#############################################

# Start day of simulation: April 1st
simDate <- as.Date("2020-04-01")
startofSimDay <- as.numeric(simDate) - as.numeric(as.Date("2020-01-01"))
# Date current physical distancing is relaxed
kbDay1 <- as.numeric(as.Date("2020-05-25")) - as.numeric(as.Date("2020-01-01")) - startofSimDay
kbDay2 <- as.numeric(as.Date("2020-06-14")) - as.numeric(as.Date("2020-01-01")) - startofSimDay
kbDay3 <- as.numeric(as.Date("2020-07-01")) - as.numeric(as.Date("2020-01-01")) - startofSimDay

# Time span
tspan <- 1:365

# plotting parameters
plotCompList <- c("I","R")   # plots these compartments
rollM <- 1                   # number of days for rolling mean
dateBreaks <- "1 month"      # date axis formatting
titleString <- "Test title"  # plot title
xString <- "Date"            # x axis title
yString <- "Frequency"       # y axis title
fileName <- "covid04.28.2020.3.txt"
plotName <- "trajPlot04.28.2020.3."

# Trial and node parameters
numTrials <- 50 # number of trials
numTrials <- min(100,numTrials) # capping number of trials at 100
N <- 500 # number of main nodes - there must be at least two per trial for ldata
N <- min(N,500) # capping number of nodes per trial
NnumTrials <- N*numTrials
nodeTrialMat <- data.table(node=c(1:(NnumTrials)),trial=rep(1:numTrials,each=N))
NList <- split(nodeTrialMat$node,nodeTrialMat$trial)
trialPop <- 65000 # total population of a trial
maxNodePop <- floor(trialPop/N) # maximum node population
I0Pop <- 40    # population that is infected in intial state
maxINodeProp <- 1/10 # proportion of nodes with one or more infected in initial state

# parachute and transfer parameters

# random chisquare distribution for selecting parachute events; truncate at no more than tspan;
paraChi_df <- 4 # assumes parachuting starts slow (because of travel restrictions), grows, then drops off as global prevalence decreases
                # decrease paraChi_df to left-skew; increase df to right-skew
parachuteDist <- ceiling(rchisq(10000,df=paraChi_df)*max(tspan)/10)
parachuteDist <- parachuteDist[parachuteDist < max(tspan)]

parachuteRate <- 1/30 # expected rate of parachuters
parachuteNum <- 1 # assume parachuting one infected at a time
transferRate <- 1/7 # expected rate of transfers
transferMaxProp <- .1 # max proportion of people to transfer nodes
transferMinProp <- .01 # min proportion of people to transfer nodes

# Model parameters
# compartments
compartments <- c("S","E","I","R","Im","M")

# student event parameters - relies on compartments
# September 20th (day students return) = day 264 in calendar year
studentReturnDate <- as.numeric(as.Date("2020-09-21")) - as.numeric(as.Date("2020-01-01")) - startofSimDay
studentReturnSpread <- 3 # symmetric spread around return date in number of days
studentPop <- 25000
maxStudentNodes <- min(ceiling(N/3),100) # maximum number of nodes students can enter
sSProp <- .9 # proportion of students who are susceptible
sEProp <- .0001 # proportion of students who are exposed
sIProp <- .001 # proportion of students who are infected
sRProp <- (1-(sSProp+sEProp+sIProp))*.9 # proportion of students who are recovered
sImProp <- (1-(sSProp+sEProp+sIProp))*.1 # proportion of students who are immune
sMProp <- 0
studentPropTable <- data.frame(compartment = compartments, frac = c(sSProp,sEProp,sIProp,sRProp,sImProp,sMProp))

# More model parameters
# global parameters
gdata = data.frame(
  beta = 2.25/12,    # Transmission rate = R0/gamma
  gamma = 1/12,      # Reciprocal of infectious period
  sigma = 1/4,       # Reciprocal of exposed period
  eta = .024,        # Case fatality rate
  delta = .1,        # Proportion of recovered who become re-susceptible
  kappa = .01,       # Reciprocal of recovered period before re-susceptibility
  mu = 0,            # Natural birth/immigration rate
  nu = 0             # Non-Covid death rate
)

# Transitions
transitions <- c(
  "@ -> mu*(S+E+I+R+Im) -> S",                                        # natural birth rate
  "R -> delta*kappa*R -> S",                                          # Recovereds who become susceptible again
  "S -> nu*S -> @",                                                   # Non-Covid deaths among susceptibles
  "S -> (1/phi)*beta*betaFactor*(season)*S*I/(S+E+I+R+Im) -> E",      # Exposures, allows for effective or ineffective policies and seasonal variation in beta
  "E -> nu*E -> @",                                                   # Non-Covid deaths among exposed
  "E -> sigma*E -> I",                                                # Development of disease
  "I -> nu*I -> @",                                                   # Non-Covid deaths among infected
  "I -> (1-eta)*gamma*I -> R",                                        # Recovery
  "I -> eta*gamma*I -> M",                                            # Death from COVID
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
  I = I0, # random number of initial infected, allows for 0, capped at maxIPop for each trial
  R = rep(0,NnumTrials),
  Im = rep(0,NnumTrials),
  M=rep(0,NnumTrials)
)

# local parameters. Used only in pts_fun (not in transitions)
# Constructing ldata
# diffusion matrix

destC <- rep.int(0:(NnumTrials-1),N)
trialC <- rep.int(1:numTrials,N)
trialC <- sort(trialC,method="quick")
trialC <- rep.int(trialC,N)
rankC <- c(1:(N*NnumTrials))
nodeC <- rep.int(0:(NnumTrials-1),N)
nodeC <- sort(nodeC,method="quick")

dMat <- data.table(dest = destC)
dMat[,trial:=trialC]
# dMat <- data.table(dest=rep.int(0:(NnumTrials-1),N),trial=rep.int(trialC,N))
dMat[,rk:=rankC]
setorder(dMat,trial,rk)
dMat[,node:=nodeC]
dMat <- dMat[,.(node,dest)]
dMat <- dMat[node!=dest]

rm(destC)
rm(trialC)
rm(rankC)
rm(nodeC)

# destC <- rep.int(0:(NnumTrials-1),N)
# trialC <- sort(rep.int(1:numTrials,N),method="quick")
# trialC <- rep.int(trialC,N)
# rankC <- sort(rep.int(1:NnumTrials,N),method="quick")
# nodeC <- sort(rep.int(0:(NnumTrials-1),N),method="quick")
d <- rbind(matrix(rep(dMat$dest,each=2),nrow=(N-1)*2),matrix(rep(c(-1,0),NnumTrials),nrow=2))
d[c(1:(N-1))*2,] <- 1
rm(dMat)

# ldata # If changing ldata variables need to update pts_funScript
parms <- data.frame(
  maxPrev1 = rep(10,NnumTrials),                 # maximum prevalence allowed before minor intervention
  maxPrev2 = rep(20,NnumTrials),                 # maximum prevalence before major intervention
  phiFactor1 = (gdata$beta/gdata$gamma)/1.15,     # intermediate reduction in beta from policy response
  phiFactor2 = (gdata$beta/gdata$gamma)/.75       # large reduction in beta from policy response
)

names(parms) <- NULL

ldata <- rbind(t(parms),d)
rm(d)

# post time step function parameters
cosAmp <- 0            # amplitude of cosine wave for seasonality of beta
# startofSimDay
upDelay <- 7             # number of days that beta-reduction measures are delayed after prevalence passes threshold
downDelay <- 28           # number of days that beta-reduction measures are relaxed after prevalence decreases past threshold
switchOffPolicies <- 1    # logical variable to switch off policies (used for counterfactuals)
switchOffDay <- as.numeric(as.Date("2020-07-02")) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
phiMoveUp <- .25        # rate at which phi converges up to phiFactor; in [0,1], bigger is faster convergence
phiMoveDown <- .1       # rate at which phi converge down to phiFactor; in [0,1], bigger is faster convergence

# Initialized continuous variables
v0 = data.frame(
  phi = rep(3,NnumTrials),                   # initial beta reduction factor (larger phi = more reduction)
  prevUp01 = rep(upDelay+1,NnumTrials),      # initializedvariable for capturing delay in response when prevalence increases
  prevUp12 = rep(upDelay+1,NnumTrials),      # variable for capturing delay in response when prevalence increases
  prevDown21 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  prevDown10 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  kbPhase = rep(0,NnumTrials),               # logical: 1 means current physical distancing, 0 means current physical distancing completely lifted
  policy = rep(1,NnumTrials),                # logical: 1 means policies can take effect; 0 means they don't
  season = rep(1,NnumTrials),                # seasonality factor for beta
  prev = rep(0,NnumTrials),                  # prevalence tracker
  prevState = rep(2,NnumTrials),             # previous state: 0 = baseline, 1 = low intervention, 2 = high intervention, used for low intervention logic
  betaFactor = rep(runif(numTrials,min=.9,max=1.1),each=N), # randomizer for beta,
  baytat = rep(2.25/12,NnumTrials)
)

# building pts_fun
pts_fun <- pts_funScript(
                         cosAmp = cosAmp,
                         startDay = startofSimDay,
                         kbDay1 = kbDay1,
                         kbDay2 = kbDay2,
                         kbDay3 = kbDay3,
                         parmCol = ncol(parms),
                         numComp = length(compartments),
                         prevType = 0,
                         upDelay = upDelay,
                         downDelay = downDelay,
                         switchOffPolicies = switchOffPolicies,
                         switchOffDay = switchOffDay,
                         phiMoveUp = phiMoveUp,
                         phiMoveDown = phiMoveDown
                         )



############### END OF PARAMETERS ##################
############### END OF PARAMETERS ##################
############### END OF PARAMETERS ##################
############### END OF PARAMETERS ##################

# Events
# create events: occasional infections parachute into the nodes; occasional transfers between the nodes. Students arrive in a wave

# Events matrix
# leaves of M from column 7
E <- cbind(diag(length(compartments)),c(rep(1,length(compartments)-2),0,0))
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
    select = 2,
    shift = 0
  )
}
###### TRANSFER EVENTS ######
# Follows a poisson process
# Partition transfers by trials
numTransferList <- lapply(1:numTrials,function(x) rpois(1,transferRate*max(tspan)))

if(max(unlist(numTransferList))>0) {
  # Transfer data frame
  # first non-students
  transferEvents <- data.frame(
    event = "extTrans",
    time = unlist(lapply(numTransferList, function(x) sample(tspan,x,replace=TRUE))),
    node = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numTransferList[[x]],replace=TRUE))),
    dest = unlist(lapply(c(1:numTrials),function(x) sample(NList[[x]],numTransferList[[x]],replace=TRUE))),
    n = 0,
    proportion = unlist(lapply(numTransferList,function(x) round(runif(x,min=transferMinProp,max=transferMaxProp),3))),
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
if(nrow(allEvents)>0) {
  model <- mparse(transitions = transitions, compartments = compartments, events = allEvents, E=E,
                gdata = gdata,ldata=ldata, u0 = u0, v0=v0, tspan = tspan, pts_fun = pts_fun)
} else {
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = gdata,ldata=ldata, u0 = u0, v0=v0, tspan = tspan, pts_fun = pts_fun)
}

set.seed(123)
set_num_threads(1)
result <- run(model)

trajTable <- data.table(trajectory(result))

setwd("./trajectories")
fwrite(trajTable,file=fileName,sep=",")

# drawing plot of trajectory
trajPlotI <- simInfPlottingFunction(
  dt= trajTable,          # trajectory of the result
  compSelect= c("I"),                        # compartments that will be plotted
  rollM = rollM,                             # number of days for rolling mean
  nTM = nodeTrialMat,                        # node-Trial matrix
  tS = tspan,                                # length of simulation
  enn = N,                                   # number of nodes per trial
  nT = numTrials,                            # number of trials in simulation
  startDate = startofSimDay,                 # start date of simulation
  dateBreaks = "1 month",                    # plot parameter: Date axis format
  titleString = "Daily active infections in Benton County",   # plot parameter: Title of plot
  xString = "Date",                          # plot parameter: Title of x axis
  yString = "Number of infections",          # plot parameter: Title of y axis
  lString = "Compartment"                    # plot parameter: Title of legend
)

ggsave(paste0(plotName,"I.png"),trajPlotI,width=9,height=5,units="in")

trajPlotPhi <- simInfPlottingFunction(
  dt= trajTable,          # trajectory of the result
  compSelect= c("phi"),                        # compartments that will be plotted
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

trajPlotPhi <- trajPlotPhi+
              geom_vline(xintercept=as.numeric(as.Date("2020-07-07")),color="darkseagreen",size=1.5)

ggsave(paste0(plotName,"phi.png"),trajPlotPhi,width=9,height=5,units="in")
