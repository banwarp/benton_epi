# schoolPassThrough_07302020.R

setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("covid19_SimInf_School_10.07.2020.R")

# parachute rate at 10 cases/100,000k per week for a community of 25,000
# assume 1/10 of cases are known, so multiply by 10
casesPerWeekBase = 5 * 10 / 100000
parachuteRateBase = casesPerWeekBase / 7

# model 1: fullmixing lowprev noquarantine nativescale
# model 2: fullmixing lowprev withquarantine nativescale
# model 2.1: fullmixing lowprev withquarantine model1scale
# model 3: fullcohorting lowprev withquarantine nativescale
# model 4: twodayaweek lowprev withquarantine nativescale
# model 4.1: twodayaweek lowprev withquarantine model3scale

set.seed(123)

# model 1
elem250FullMixingLowPrev <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullMixingLowPrevNoQNoTestNativeScale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  studentTeacherDiff = c(1,1,1,1), 
  classroomThreshold = NULL,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = NULL,                      # Number of active detected cases before grade quarantined
  schoolThreshold = NULL,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  infectionCount = c(1), 
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.25,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.25,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4), #
  
  miny = 0,
  maxy = 225,
  
  plotRandomTrajs = 3, 
  
  titleString = "Elementary school; 250 students; 50 staff; 5 days a week",
  subtitleString = "Full mixing; low prevalence; no quarantines; no testing after identified case"
)

################################################
################################################
################################################
################################################
################################################

# model 2
set.seed(111)
elem250FullMixingLowPrevWQ <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullMixingLowPrevWithQNoTestNativeScale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  infectionCount = c(1),
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  studentTeacherDiff = c(1,1,1,1),
  classroomThreshold = 1,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = NULL,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 2,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.25,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.25,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4), #
  
  plotRandomTrajs = 3,
  
  titleString = "Elementary school; 250 students; 50 staff; 5 days a week",
  subtitleString = "Full mixing; low prevalence; with quarantines; no testing after identified case"
)

################################################
################################################
################################################
################################################
################################################

# model 2.1
elem250FullMixingLowPrevWQ <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullMixingLowPrevWithQNoTestModel1Scale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  infectionCount = c(1),
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  studentTeacherDiff = c(1,1,1,1),
  classroomThreshold = 1,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = NULL,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 2,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.25,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.25,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4), #
  
  miny = 0,
  maxy = 225,
  
  plotRandomTrajs = 3,
  
  titleString = "Elementary school; 250 students; 50 staff; 5 days a week",
  subtitleString = "Full mixing; low prevalence; with quarantines; no testing after identified case"
)

################################################
################################################
################################################
################################################
################################################

# model 3
elem250FullCohortingLowPrevWQ <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullCohortingLowPrevWithQNoTestNativeScale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  infectionCount = c(1),
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  studentTeacherDiff = c(1,1,1,1),
  classroomThreshold = 1,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = NULL,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 2,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = rep(0,4),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = rep(0,4), # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = rep(0,4),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = rep(0,4),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4), #
  
  miny = 0,
  maxy = 12,
  
  plotRandomTrajs = 3,
  
  titleString = "Elementary school; 250 students; 50 staff; 5 days a week",
  subtitleString = "Full cohorting; low prevalence; with quarantines; no testing after identified case"
)

################################################
################################################
################################################
################################################
################################################

# model 4
elem250TwiceAWeekLowPrevWQ <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem2502xWeekLowPrevWithQNoTestNativeScale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  infectionCount = c(1),
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  studentTeacherDiff = c(1,1,1,1),
  classroomThreshold = 1,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = 2,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 4,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  weekDays = list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),
                  c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,0,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,0,0), # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.1,0,0,0),   # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.1,0,0,0),   # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = c(0,0,1/8,0),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = c(0,0,.25,0),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = c(0,0,.3,0),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = c(0,0,.3,0), #
  
  plotRandomTrajs = 3,
  
  titleString = "Elementary school; 250 students; 50 staff; 2 groups x 2 days a week",
  subtitleString = "Limited mixing; low prevalence; with quarantines; no testing after identified case"
)

################################################
################################################
################################################
################################################
################################################

# model 4.1
elem250TwiceAWeekLowPrevWQ <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem2xWeekLowPrevWithQNoTestModel3Scale"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 200,                         # Number of trials in simulation
  
  infectionCount = c(1),
  infectionNodes = c(1),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c("2020-09-01"),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(0),                   # Symmetric spread in days of infection spreader infections
  infectionST = c(1),                       # Infection group type: 1=stationary or 2=transitory
  
  R0Base = 4.1, #
  R0Spread = 0.5, #
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  studentTeacherDiff = c(1,1,1,1),
  classroomThreshold = 1,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = 2,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 4,                    # Number of active detected cases before school quarantined
  cohortTesting = 0,                         # Do not test cohort after identified case
  
  weekDays = list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),
                  c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,0,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,0,0), # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.1,0,0,0),   # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.1,0,0,0),   # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = c(0,0,1/8,0),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = c(0,0,.25,0),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = c(0,0,.3,0),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = c(0,0,.3,0), #
  
  miny = 0, 
  maxy = 12, 
  
  plotRandomTrajs = 3,
  
  titleString = "Elementary school; 250 students; 50 staff; 2 groups x 2 days a week",
  subtitleString = "Limited mixing; low prevalence; with quarantines; no testing after identified case"
)

