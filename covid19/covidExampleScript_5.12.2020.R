# covidPassThroughScript.R

# Pass through script for more efficient scenario generation

library(data.table)

setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("covid19_SimInf_5.12.2020.R")

countyDT <- fread("https://raw.githubusercontent.com/banwarp/benton_epi/master/covid19/countyMetaData_2020.05.04.csv")

countyDT[,simID:=paste0("covid5.12.2020",County)]
countyDT[,I0Pop:= round(.003*Population,0)]
countyDT[,N:=ceiling(Population/1000)]

#### All County Scenarios - long run time

baseLineList <- vector(mode = "list", length = nrow(countyDT))
noPoliciesList <- baseLineList
superSpreaderList <- baseLineList

for(i in 1:nrow(countyDT)){
  baseLineList[[i]] <- covidWrapper(simID = paste0(countyDT[i,(simID)],".BaseLine"),
                         trialPop = countyDT[i,(Population)],
                         I0Pop = countyDT[i,(I0Pop)],
                         N = countyDT[i,(N)],
                         titleString = countyDT[i,(County)]
                       )
}

for(i in 1:nrow(countyDT)){
  noPoliciesList[[i]] <- covidWrapper(simID = paste0(countyDT[i,(simID)],".NoPolicies"),
                          trialPop = countyDT[i,(Population)],
                          I0Pop = countyDT[i,(I0Pop)],
                          N = countyDT[i,(N)],
                          titleString = countyDT[i,(County)],
                          switchOffPolicies = 1
                         )
}

for(i in 1:nrow(countyDT)){
  superSpreaderList[[i]] <- covidWrapper(simID = paste0(countyDT[i,(simID)],".SuperSpreader"),
                             trialPop = countyDT[i,(Population)],
                             I0Pop = countyDT[i,(I0Pop)],
                             N = countyDT[i,(N)],
                             titleString = countyDT[i,(County)],
                             superNodes = min(floor(countyDT[i,(N)]/2),5),
                             superDate = "2020-11-01",
                             superInfections = 50,
                             superSpread=2
                            )
}


#### Benton County Scenarios

i <- 2

bentonBaseLine <- covidWrapper(simID=paste0(countyDT[i,(simID)],".baseLine"),
                   trialPop = countyDT[i,(Population)],
                   I0Pop = countyDT[i,(I0Pop)],
                   N = 500,
                   titleString = countyDT[i,(County)]
)

bentonBestGuess <- covidWrapper(simID=paste0(countyDT[i,(simID)],".bestGuess"),
                     trialPop = countyDT[i,(Population)],
                     I0Pop = countyDT[i,(I0Pop)],
                     N = 500,
                     maxStudentNodes = 100,
                     titleString = countyDT[i,(County)]
)

bentonCounterFactual <- covidWrapper(simID=paste0(countyDT[i,(simID)],".CounterFactual"),
                         trialPop = 65000,
                         I0Pop = round(.0023*65000,0),
                         N = 500,
                         maxStudentNodes = 100,
                         isoRate = .125,
                         switchOffPolicies = 1,
                         titleString = countyDT[i,(County)]
)

bentonEnhancedContactTracing <- covidWrapper(simID=paste0(countyDT[i,(simID)],".enhancedContactTracing"),
                                 trialPop = 65000,
                                 I0Pop = round(.0023*65000,0),
                                 N = 500,
                                 maxStudentNodes = 100,
                                 isoRate = .125,
                                 titleString = countyDT[i,(County)]
)