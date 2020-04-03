# covid19_4.2.2020.R
library(deSolve)
library(R0)
library(ggplot2)
library(reshape2)
setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19")

# data from https://usafacts.org/visualizations/coronavirus-covid-19-spread-map/
covid <- read.csv("https://usafactsstatic.blob.core.windows.net/public/data/covid-19/covid_confirmed_usafacts.csv",header=T)

# list of states
stateList <- sort(unique(covid$State))
# filter data to just Oregon
covid <- covid[covid$State == "OR",]


# total population in Oregon
totalPop <- 4.19e6
# population of geography of interest: Benton County
subPop <- 95000

# cumulative number of cases in Oregon
stateCasesRaw <- colSums(covid[,c(5:ncol(covid))])
# most recent total number of cases
mostRecentTotal <- stateCasesRaw[length(stateCasesRaw)]
# keeping only non-zero cumulative case totals
stateCasesRaw <- as.vector(stateCasesRaw[stateCasesRaw>0])

# number of days with >0 cumulative cases
numDays <- length(stateCasesRaw)
# daySeq for use in modeling 4day, 6day, and 10day doubling rate case counts
daySeq <- c((-numDays+1):0)

# modeled number of cumulative cases based on different doubling rate assumptions.
stateCases4day <- round(mostRecentTotal*2^(daySeq/4),0)
stateCases6day <- round(mostRecentTotal*2^(daySeq/6),0)
stateCases10day <- round(mostRecentTotal*2^(daySeq/10),0)
underCountDoublingRate <- numDays*log(2)/log(mostRecentTotal/5)
stateCases6dayHalfCount <- round(mostRecentTotal*2*2^(daySeq/underCountDoublingRate),0)

# combining 4 different vectors into a list
stateCasesList <- list(stateCasesRaw,stateCases4day,stateCases6day,stateCases10day,stateCases6dayHalfCount)

# setting up parameters for SEIR model
# building a grid of parameters for machine testing
times <- c(0:80)
N <- 4.19e+06 # Total population
# Initial conditions
I_0 <- c(1:15) # Initial number of infectious
R_0 <- 0 # Initial number of recovered/immune
# S_0 <- N-I_0-R_0 # Initial number of susceptible
E_0 <- 0 # Initial number of exposed
D_0 <- 0 # Initial number of post-infectious who will die
M_0 <- 0 # Initial number of dead

R0 <- seq(from=1.5,to=3.5,by=.2) # basic reproduction number
mu <- 0 # birth rate
nu <- 0 # death rate other than COVID
sigma <- 1/c(2:6) # rate at which exposed become infected = reciprocal of exposed period
gamma <- 1/c(6:12) # reciprocal of infectious period

#############
#############
#############
# ###### FOR TESTING OTHERWISE COMMENT OUT
# # setting up parameters
# times <- c(0:80)
# N <- 4.19e+06 # Total population
# # Initial conditions
# I_0 <- c(10:12) # Initial number of infectious
# R_0 <- 0 # Initial number of recovered/immune
# # S_0 <- N-I_0-R_0 # Initial number of susceptible
# E_0 <- 0 # Initial number of exposed
# D_0 <- 0 # Initial number of post-infectious who will die
# M_0 <- 0 # Initial number of dead
# 
# R0 <- seq(from=2,to=2.4,by=.2) # basic reproduction number
# mu <- 0 # birth rate
# nu <- 0 # death rate other than COVID
# sigma <- 1/c(4:5) # rate at which exposed become infected = reciprocal of exposed period
# gamma <- 1/c(8:9) # reciprocal of infectious period
#############
#############
#############

parameterGrid <- expand.grid(I_0,R0,mu,nu,sigma,gamma)
names(parameterGrid) <- c("I_0","R0","mu","nu","sigma","gamma")
parameterGrid$bayta <- parameterGrid$R0*parameterGrid$gamma # rate at which susceptible become exposed

# number of models to keep; the n best models
numModels <- 5

# A function that lets Rhat decrease over time
RhatFunction <- function(t,rNaught) {
  if(t <numDays) {rhat <- rNaught}
  else if(t<numDays+10) {rhat <- rNaught*.9}
  else if(t<numDays+20) {rhat <- rNaught*.75}
  else {rhat <- exp((60-t)*log(1.04))*(rNaught*.75-.9)+.9}
}

# Function that builds dim(parameterGrid) number of models and chooses the 5 models with the smallest deviation from the data
# Then the model uses those parameters as the initial state and predicts the Oregon dynamics based on letting R0 decrease.
# Output is the set of 5 best parameters, 5 SEIR models with constant R0, and 5 SEIR models with stepped down R0.
# seirFunction runs on each set of stateCases: Raw, 4day, 6day, 10day
seirFunction <- function(stateCasesVector,numModels,totalPop) {
  errorGrid <- data.frame(index = c(1:nrow(parameterGrid)),error=0)
  
  for(i in 1:nrow(parameterGrid)) {
    I_0 <- parameterGrid[i,"I_0"]
    S_0 <- N-I_0-R_0
    bayta <- parameterGrid[i,"bayta"]
    mu <- parameterGrid[i,"mu"]
    nu <- parameterGrid[i,"nu"]
    sigma <- parameterGrid[i,"sigma"]
    gamma <- parameterGrid[i,"gamma"]
    
    init <- c(S=S_0,E=E_0,I=I_0,R=R_0) # initial conditions
    parameters <- c(N,bayta,mu,nu,sigma,gamma) # parameters
    
    seirModel <- function(time,state,parameters) {
      with(as.list(c(state,parameters)), {
        dS = mu*N-nu*S-(bayta*S*I)/N
        dE = (bayta*S*I)/N - nu*E - sigma*E
        dI = sigma*E - gamma*I - nu*I
        dR = gamma*I
        
        return(list(c(dS,dE,dI,dR)))
      })
    }
    
    seirModel1 <- data.frame(ode(y = init, times = times, func = seirModel, parms = parameters))
    
    infectionTime <- 1/sigma+1/gamma
    
    infectedTotal <- seirModel1$R[c(1:(numDays+infectionTime))]
    
    errorGrid[i,2] <- norm(matrix(infectedTotal-c(rep(1,infectionTime),stateCasesVector)),type="f")/infectionTime
  }
  
  errorGrid <- errorGrid[order(errorGrid$error),]
  
  bestGrid <- parameterGrid[errorGrid[c(1:numModels),"index"],]
  
  times <- c(0:185) # six months
  N <- totalPop # Total population
  # Initial conditions
  R_0 <- 0 # Initial number of recovered/immune
  E_0 <- 0 # Initial number of exposed
  D_0 <- 0 # Initial number of post-infectious who will die
  M_0 <- 0 # Initial number of dead
  t_0 <- 0 # Initial time
  
  constantR0Grid <- matrix(data=0,nrow=186,ncol=numModels)
  variableR0Grid <- matrix(data=0,nrow=186,ncol=numModels)
  
  for(i in 1:nrow(bestGrid)) {
    I_0 <- bestGrid[i,"I_0"]
    S_0 <- N-I_0-R_0
    R0 <- bestGrid[i,"R0"]
    bayta <- bestGrid[i,"bayta"]
    mu <- bestGrid[i,"mu"]
    nu <- bestGrid[i,"nu"]
    sigma <- bestGrid[i,"sigma"]
    gamma <- bestGrid[i,"gamma"]
    
    initConstant <- c(S=S_0,E=E_0,I=I_0,R=R_0) # initial conditions
    initVariable <- c(S=S_0,E=E_0,I=I_0,R=R_0,t=t_0) # initial conditions
    
    parameters <- c(N,bayta,mu,nu,sigma,gamma) # parameters
    
    seirConstantModel <- function(time,state,parameters) {
      with(as.list(c(state,parameters)), {
        dS = mu*N-nu*S-(bayta*S*I)/N
        dE = (bayta*S*I)/N - nu*E - sigma*E
        dI = sigma*E - gamma*I - nu*I
        dR = gamma*I
        
        return(list(c(dS,dE,dI,dR)))
      })
    }
    
    seirVariableModel <- function(time,state,parameters) {
      with(as.list(c(state,parameters)), {
        bayta = RhatFunction(t,R0)*gamma
        dS = mu*N-nu*S-(bayta*S*I)/N
        dE = (bayta*S*I)/N - nu*E - sigma*E
        dI = sigma*E - gamma*I - nu*I
        dR = gamma*I
        dt = 1
        
        return(list(c(dS,dE,dI,dR,dt)))
      })
    }
    
    seirConstantModel1 <- data.frame(ode(y = initConstant, times = times, func = seirConstantModel, parms = parameters))
    seirVariableModel1 <- data.frame(ode(y = initVariable, times = times, func = seirVariableModel, parms = parameters))
    
    
    constantR0Grid[,i] <- seirConstantModel1$I
    variableR0Grid[,i] <- seirVariableModel1$I
  }
  
  return(list(bestGrid,constantR0Grid,variableR0Grid))
    
  
}

# applying seirFunction to each set of stateCases
# This is when seirFunction is actually run for each of the 4 stateCases vectors
seirResultList <- lapply(stateCasesList,seirFunction,numModels=numModels,totalPop=totalPop)

# Captions for plotting
captionList <- c("Using raw daily totals, with R0 dropping from initial value to .9",
                 "Assuming 4 day doubling rate with initial missed cases, with R0 dropping from initial value to .9 beginning  ",
                 "Assuming 6 day doubling rate with initial missed cases, with R0 dropping from initial value to .9",
                 "Assuming 10 day doubling rate with initial missed cases, with R0 dropping from initial value to .9",
                 "Assuming 12 initial cases with exponential growth rate to twice the diagnosed cases, with R0 dropping from initial value to .9")

# plotting the 5 best models for each of the 5 stateCases vectors

setwd("./images")

plotFunction <- function(index,captionList) {
  indexList <- seirResultList[[index]]
  parameterSet <- data.frame(indexList[[1]])
  predictions <- data.frame(indexList[[3]])
  colnames(predictions) <- c("First","Second","Third","Fourth","Fifth")
  predictions$Actual <- NA
  infectionTime <- 1/parameterSet[1,"sigma"] + 1/parameterSet[1,"gamma"]
  predictions[c((infectionTime+1):(infectionTime+length(stateCasesList[[index]]))),"Actual"] <- stateCasesList[[index]]
  
  minDate <- as.Date("2020/2/29")-min(which(!is.na(predictions$Actual)))+1
  predictions$Date <- seq(as.Date(minDate),by="day",length.out=186)
  
  predictions <- melt(predictions,id.vars="Date")
  
  predictions$value <- predictions$value*subPop/totalPop
  predictions$hospitalizations <- predictions$value*.15
    
  captionText <- captionList[index]
  
  p1 <- ggplot(predictions, aes(x=Date,y=value))+
    geom_point(data=predictions)+
    geom_point(data=predictions[predictions$variable=="Actual",],aes(color="red"))+
    labs(title="COVID-19 daily infections",x="Date",y="Number of infections",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  p2 <- ggplot(predictions, aes(x=Date,y=hospitalizations))+
    geom_point(data=predictions)+
    geom_point(data=predictions[predictions$variable=="Actual",],aes(color="red"))+
    labs(title="COVID-19 daily hospital bed use, assuming 15% of infections need hospital",x="Date",y="Number of hospital beds",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  ggsave(paste0("covidPlot",index,".png"),p1,width=12,height=8,units="in")
  ggsave(paste0("covidHosps",index,".png"),p2,width=12,height=8,units="in")
}

lapply(c(1:length(stateCasesList)),plotFunction,captionList=captionList)
