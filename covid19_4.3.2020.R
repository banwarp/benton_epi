# covid19_4.3.2020.R
# key changes from covid19_4.2.2020.R:
# Updated data (one more day of data)
# Introduced case fatality modeling
# Modeling Benton County using N = 95000 instead of N=4.19e6 then multiplying resulting data by 95k/4.19M
# Changed hospitalization rates: 20.8% hospitalization rate and 6.3% ICU rate
# Changed Rhat function to match IDM model
# Changed destination folder for images to images/4.3.2020


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

totalDays <- 185 # six months
totalDaysPlusOne <- totalDays+1 # for correct date lengths
eta <- 5883/242004 # case fatality rate = 2.4%
hospitalizationRate <- .25 # proportion of cases in need of hospitalization (includes ICU) from: https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
icuRate <- .07 # proportion of cases in need of ICU from: https://www.cdc.gov/mmwr/volumes/69/wr/mm6912e2.htm
hospitalizationPeriodProportion <- .75 # estimated proportion of symptomatic period spent in hospital for hospitalized; to scale the numbers of beds needed
icuPeriodProportion <- .5 # estimated proportion of symptomatic period spent in ICU for those who need ICU

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
  else if(t<numDays+7) {rhat <- rNaught*.9}
  else if(t<numDays+14) {rhat <- rNaught*.8}
  else if(t<numDays+21) {rhat <- rNaught*.65}
  else if(t<numDays+28) {rhat <- rNaught*.45}
  else {rhat <- exp((60-t)*log(1.04))*(rNaught*.45-.9)+.9}
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
  
  times <- c(0:totalDays)
  N <- subPop # Population of specific geography
  # Initial conditions
  R_0 <- 0 # Initial number of recovered/immune
  E_0 <- 0 # Initial number of exposed
  M_0 <- 0 # Initial number of dead
  t_0 <- 0 # Initial time
  
  # constantR0Infections <- matrix(data=0,nrow=totalDaysPlusOne,ncol=numModels)
  variableR0Infections <- matrix(data=0,nrow=totalDaysPlusOne,ncol=numModels)
  # constantR0Deaths <- constantR0Infections
  variableR0Deaths <- variableR0Infections
  
  for(i in 1:numModels) {
    I_0 <- ceiling(bestGrid[i,"I_0"]*N/totalPop) # scaling to geography population, minimum I_0 = 1
    S_0 <- N-I_0-R_0
    R0 <- bestGrid[i,"R0"]
    bayta <- bestGrid[i,"bayta"]
    mu <- bestGrid[i,"mu"]
    nu <- bestGrid[i,"nu"]
    sigma <- bestGrid[i,"sigma"]
    gamma <- bestGrid[i,"gamma"]
    
    # initConstant <- c(S=S_0,E=E_0,I=I_0,R=R_0,M_0) # initial conditions
    initVariable <- c(S=S_0,E=E_0,I=I_0,R=R_0,M=M_0,t=t_0) # initial conditions
    
    parameters <- c(N,bayta,mu,nu,sigma,gamma,eta) # parameters
    
    # seirConstantModel <- function(time,state,parameters) {
    #   with(as.list(c(state,parameters)), {
    #     dS = mu*N-nu*S-(bayta*S*I)/N
    #     dE = (bayta*S*I)/N - nu*E - sigma*E
    #     dI = sigma*E - gamma*I - nu*I
    #     dR = gamma*I*(1-eta)
    #     dM = gamma*I*eta
    #     
    #     return(list(c(dS,dE,dI,dR,dM)))
    #   })
    # }
    
    seirVariableModel <- function(time,state,parameters) {
      with(as.list(c(state,parameters)), {
        bayta = RhatFunction(t,R0)*gamma
        dS = mu*N-nu*S-(bayta*S*I)/N
        dE = (bayta*S*I)/N - nu*E - sigma*E
        dI = sigma*E - gamma*I - nu*I
        dR = gamma*I*(1-eta)
        dM = gamma*I*eta
        dt = 1
        
        return(list(c(dS,dE,dI,dR,dM,dt)))
      })
    }
    
    # seirConstantModel1 <- data.frame(ode(y = initConstant, times = times, func = seirConstantModel, parms = parameters))
    seirVariableModel1 <- data.frame(ode(y = initVariable, times = times, func = seirVariableModel, parms = parameters))
    
    
    # constantR0Infections[,i] <- seirConstantModel1$I
    variableR0Infections[,i] <- seirVariableModel1$I
    variableR0Deaths[,i] <- seirVariableModel1$M
  }
  
  return(list(bestGrid,variableR0Infections,variableR0Deaths))
    
  
}

# applying seirFunction to each set of stateCases
# This is when seirFunction is actually run for each of the 5 stateCases vectors
seirResultList <- lapply(stateCasesList,seirFunction,numModels=numModels,totalPop=totalPop)


######################
######################
######################
# estimating hospitalizations and plotting curves
# Captions for plotting
captionList <- c("Using raw daily totals, with R0 dropping from initial value to .9",
                 "Assuming 4 day doubling rate with initial missed cases, with R0 dropping from initial value to .9 beginning  ",
                 "Assuming 6 day doubling rate with initial missed cases, with R0 dropping from initial value to .9",
                 "Assuming 10 day doubling rate with initial missed cases, with R0 dropping from initial value to .9",
                 "Assuming 1 local initial cases with exponential growth rate to twice the diagnosed cases, with R0 dropping from initial value to .9")

GSBeds <- 158 # regular number of Good Sam beds
GSSurge <- 225 # surged number of Good Sam beds
GSMax <- 250 # Estimated max number of clinical beds

# plotting the 5 best models for each of the 5 stateCases vectors

setwd("./images4.3.2020")

plotFunction <- function(index,captionList) {
  indexList <- seirResultList[[index]]
  parameterSet <- data.frame(indexList[[1]])
  
  #Infections and Hospitalizations
  infectionPredictions <- data.frame(indexList[[2]])
  colnames(infectionPredictions) <- c("Raw","4day","6day","10day","6dayHalfCount")
  
  # cumulative infections
  cumulativeInfectionPredictions <- cumsum(data.frame(t(apply(infectionPredictions,1,function(x) x*parameterSet$gamma))))
  
  infectionPredictions$Actual <- NA
  infectionTime <- 1/parameterSet[1,"sigma"] + 1/parameterSet[1,"gamma"]
  infectionPredictions[c((infectionTime+1):(infectionTime+length(stateCasesList[[index]]))),"Actual"] <- stateCasesList[[index]]*subPop/totalPop
  
  minDate <- as.Date("2020/2/29")-min(which(!is.na(infectionPredictions$Actual)))+1
  infectionPredictions$Date <- seq(as.Date(minDate),by="day",length.out=totalDaysPlusOne)
  
  # cumulative infections date range
  cumulativeInfectionPredictions$Date <- seq(as.Date(minDate),by="day",length.out=totalDaysPlusOne)
  
  infectionPredictions <- melt(infectionPredictions,id.vars="Date")
  cumulativeInfectionPredictions <- melt(cumulativeInfectionPredictions,id.vars="Date")
  
  # infectionPredictions$value <- infectionPredictions$value*subPop/totalPop
  infectionPredictions$hospitalizations <- infectionPredictions$value*hospitalizationRate*hospitalizationPeriodProportion
  infectionPredictions$ICU <- infectionPredictions$value*icuRate*icuPeriodProportion
  infectionPredictions <- infectionPredictions[-which(infectionPredictions$variable=="Actual"),] # dropping the actual data since it isn't scaled
    
  
  #########################
  # Deaths
  deathPredictions <- data.frame(indexList[[3]])
  colnames(deathPredictions) <- c("Raw","4day","6day","10day","6dayHalfCount")
  dailyDeathPredictions <- data.frame(apply(deathPredictions,2,diff))
  
  deathPredictions$Date <- seq(as.Date(minDate),by="day",length.out=totalDaysPlusOne)
  dailyDeathPredictions$Date <- seq(as.Date(minDate),by="day",length.out=totalDays)
  
  deathPredictions <- melt(deathPredictions,id.vars="Date")
  dailyDeathPredictions <- melt(dailyDeathPredictions,id.vars="Date")
  
  
  ##### Plotting
  
  captionText <- captionList[index]
  
  infectionsPlot <- ggplot(infectionPredictions, aes(x=Date,y=value))+
    geom_point(data=infectionPredictions)+
    labs(title="COVID-19 daily infections",x="Date",y="Number of infections",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  cumulativeInfectionsPlot <- ggplot(cumulativeInfectionPredictions, aes(x=Date,y=value))+
    geom_point(data=cumulativeInfectionPredictions)+
    labs(title="COVID-19 cumulative infections",x="Date",y="Number of infections",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  hospitalizationsPlot <- ggplot(infectionPredictions, aes(x=Date,y=hospitalizations))+
    geom_point(data=infectionPredictions)+
    geom_hline(yintercept=GSBeds,color="blue",size=2)+
    geom_hline(yintercept=GSSurge,color="purple",size=2)+
    geom_hline(yintercept=GSMax,color="red",size=2)+
    geom_label(aes(x=minDate,y=GSBeds,label = paste0("Good Sam bed capacity: ",GSBeds),vjust=1.2, hjust=0))+
    geom_label(aes(x=minDate,y=GSSurge,label = paste0("Good Sam surge capacity: ",GSSurge),vjust=1.2, hjust=0))+
    geom_label(aes(x=minDate,y=GSMax,label = paste0("Good Sam max capacity: ",GSMax),vjust=-.1, hjust=0))+
    labs(title="COVID-19 daily hospital bed use, assuming 25% of infections need hospital",x="Date",y="Number of hospital beds",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  icuPlot <- ggplot(infectionPredictions, aes(x=Date,y=ICU))+
    geom_point(data=infectionPredictions)+
    labs(title="COVID-19 daily ICU beds, assuming 7% of infections need hospital",x="Date",y="Number of ICU beds",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  deathPlot <- ggplot(deathPredictions, aes(x=Date,y=value))+
    geom_point(data=deathPredictions)+
    labs(title="COVID-19 cumulative deaths, assuming 2.4% of infected die",x="Date",y="Number of deaths",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  dailyDeathPlot <- ggplot(dailyDeathPredictions, aes(x=Date,y=value))+
    geom_point(data=dailyDeathPredictions)+
    labs(title="COVID-19 daily deaths, assuming 2.4% of infected die",x="Date",y="Number of deaths",caption = captionText)+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(angle=45))+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          plot.caption=element_text(hjust=0))
  
  ggsave(paste0("covidInfections",index,".png"),infectionsPlot,width=12,height=8,units="in")
  ggsave(paste0("covidCumulativeInfections",index,".png"),cumulativeInfectionsPlot,width=12,height=8,units="in")
  ggsave(paste0("covidHosps",index,".png"),hospitalizationsPlot,width=12,height=8,units="in")
  ggsave(paste0("covidICU",index,".png"),icuPlot,width=12,height=8,units="in")
  ggsave(paste0("covidDeaths",index,".png"),deathPlot,width=12,height=8,units="in")
  ggsave(paste0("covidDailyDeaths",index,".png"),dailyDeathPlot,width=12,height=8,units="in")
  
}

lapply(c(1:length(stateCasesList)),plotFunction,captionList=captionList)
