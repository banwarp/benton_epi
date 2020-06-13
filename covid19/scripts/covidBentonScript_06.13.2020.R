# covidBentonScript_06.13.2020.R

# Pass through script for more efficient scenario generation

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
                            paraMu = 1.2,
                            paraSig = 1.25,
                            monitoringSuccess = .125,
                            maxMassEntryNodes = 0,
                            titleString = "Benton County: No OSU"
)

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



r6 <- covidWrapper(simID="NodeGroups",
                   nodeGroupList = c(rep(1,200),rep(2,150),rep(3,100),rep(4,50)), # List of group IDs for node groups. Default is 1 group of nodes.
                   parachuteNodeGroups = c(1,2),             # Which nodes groups the parachuters can land in
                   superInfections = c(100,50),              # Number of infections caused by the super spreader
                   superNodes = c(20,20),                    # Number of nodes that the super spreader contacts
                   superNodeGroups = c(3,4),                 # Which node groups the super spreader contacts. Must use list() syntax for multiple events
                   superDate = c("2020-10-01","2021-02-01"), # Date the super spreader lands. Date can also be numeric i.e. 200
                   superSpread = c(3,5),                     # Symmetric spread in days of super spreader infections
                   plotGroups = c(1,3),                      # Which node groups to plot
                   plotSum = FALSE,                          # Whether to sum across all node groups or sum them separately
                   lString = "Node Group"                    # Title of legend
)




b1red <- covidWrapper(simID=paste0("covid",Sys.Date(),"benton_Red"),
                               trialPop = 65000,
                               I0Pop = round(.0023*65000,0),
                               N = 500,
                               maxMassEntryNodes = 100,
                               superNodes = c(10,5),
                               superDate = c("2020-10-15","2020-11-07"),
                               superInfections = c(100,100),
                               superSpread = c(0,0),
                               isoRate = .25,
                               titleString = "Benton County"
)

b1orange <- covidWrapper(simID=paste0("covid",Sys.Date(),"benton_Orange"),
                      trialPop = 65000,
                      I0Pop = round(.0023*65000,0),
                      N = 500,
                      maxMassEntryNodes = 100,
                      superNodes = 5,
                      superDate = "2020-10-15",
                      superInfections = 100,
                      superSpread = 0,
                      isoRate = .25,
                      titleString = "Benton County"
)

b1yellow <- covidWrapper(simID=paste0("covid",Sys.Date(),"benton_Yellow"),
                         trialPop = 65000,
                         I0Pop = round(.0023*65000,0),
                         N = 500,
                         maxMassEntryNodes = 100,
                         superNodes = 5,
                         superDate = "2020-10-15",
                         superInfections = 50,
                         superSpread = 0,
                         isoRate = .25,
                         titleString = "Benton County"
)

b1green <- covidWrapper(simID=paste0("covid",Sys.Date(),"benton_Green"),
                         trialPop = 65000,
                         I0Pop = round(.0023*65000,0),
                         N = 500,
                         maxMassEntryNodes = 100,
                         superNodes = 5,
                         superDate = "2020-10-15",
                         superInfections = 25,
                         superSpread = 0,
                        isoRate = .25,
                        titleString = "Benton County"
)

b1noGames <- covidWrapper(simID=paste0("covid",Sys.Date(),"benton_noGames"),
                        trialPop = 65000,
                        I0Pop = round(.0023*65000,0),
                        N = 500,
                        maxMassEntryNodes = 100,
                        isoRate = .25,
                        titleString = "Benton County"
)