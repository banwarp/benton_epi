
setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("covid19_SimInf_06.06.2020.R")

bentonStatusQuo <- covidWrapper(simID=paste0("covid",Sys.Date(),"BentonCounty.StatusQuo"),
                                             simDate = "2020-06-01",
                                             trialPop = 65000,
                                             I0Pop = round(.0015*65000,0),
                                             N=500,
                                             paraMu = 1.2,
                                             paraSig = 1.25,
                                             monitoringSuccess = .125,
                                             maxPrev1 = 30,
                                             maxPrev2 = 60,
                                             massEntryReturnDate = "2020-09-25",      # Date of mass Entry
                                             massEntryReturnSpread = 6,               # Symmetric spread of days to spread out mass entry
                                             massEntryPop = 18000,                    # Mass entry population
                                             maxMassEntryNodes=100,
                                             massEntryNodeGroups = NULL,              # Which node groups the individuals enter
                                             mSProp = .9,                             # Proportion of individuals who are susceptible
                                             mEProp = .0001,                          # Proportion of individuals who are exposed
                                             mIProp = .001,                           # Proportion of individuals who are infectious; R, Im, and M are calculated
                                             titleString = "Benton County: Status Quo"
)

bentonnoOSU <- covidWrapper(simID=paste0("covid",Sys.Date(),"BentonCounty.noOSU"),
                            simDate = "2020-06-01",
                            trialPop = 65000,
                            I0Pop = round(.0015*65000,0),
                            N = 500,
                            maxPrev1 = 30,
                            maxPrev2 = 60,
                            paraMu = 1.2,
                            paraSig = 1.25,
                            monitoringSuccess = .125,
                            maxMassEntryNodes = 0,
                            titleString = "Benton County: No OSU"
)

rm(bentonnoOSU)

bentonCounterFactual <- covidWrapper(simID=paste0("covid",Sys.Date(),"BentonCounty.CounterFactual"),
                                     simDate = "2020-06-01",
                                     trialPop = 65000,
                                     I0Pop = round(.0015*65000,0),
                                     N = 500,
                                     paraMu = 1.2,
                                     paraSig = 1.25,
                                     massEntryReturnDate = "2020-09-25",      # Date of mass Entry
                                     massEntryReturnSpread = 6,               # Symmetric spread of days to spread out mass entry
                                     massEntryPop = 18000,                    # Mass entry population
                                     maxMassEntryNodes=100,
                                     massEntryNodeGroups = NULL,              # Which node groups the individuals enter
                                     mSProp = .9,                             # Proportion of individuals who are susceptible
                                     mEProp = .0001,                          # Proportion of individuals who are exposed
                                     mIProp = .001,                          # Proportion of individuals who are infectious; R, Im, and M are calculated
                                     isoRate = .125,
                                     monitoringSuccess = .125,
                                     switchOffPolicies = 1,
                                     titleString = "Benton County: Counterfactual"
)

rm(bentonCounterFactual)

bentonPublicHealthChampions <- covidWrapper(simID=paste0("covid",Sys.Date(),"BentonCounty.publicHealthChampions"),
                                           simDate = "2020-06-01",
                                           trialPop = 65000,
                                           I0Pop = round(.0015*65000,0),
                                           N = 500,
                                           paraMu = 1.2,
                                           paraSig = 1.25,
                                           massEntryReturnDate = "2020-09-25",      # Date of mass Entry
                                           massEntryReturnSpread = 6,               # Symmetric spread of days to spread out mass entry
                                           massEntryPop = 20000,                    # Mass entry population
                                           maxMassEntryNodes=100,
                                           massEntryNodeGroups = NULL,              # Which node groups the individuals enter
                                           mSProp = .9,                             # Proportion of individuals who are susceptible
                                           mEProp = .0001,                          # Proportion of individuals who are exposed
                                           mIProp = .0015,                          # Proportion of individuals who are infectious; R, Im, and M are calculated
                                           isoRate = .5,
                                           monitoringSuccess = .5,
                                           pdDecay = -1,
                                           RPhysicalDistancing = 2.5,
                                           titleString = "Benton County: PHC"
)

rm(bentonPublicHealthChampions)
