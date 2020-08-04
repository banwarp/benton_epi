# schoolPassThrough_07302020.R

setwd("L:/Health/Epidemiology/Banwarth_Epi/covid19/scripts")
source("covid19_SimInf_School_07.21.2020.R")

# parachute rate at 10 cases/100,000k per week for a community of 25,000
# assume 1/5 of cases are known, so multiply by 5
casesPerWeekBase = 5 * 10 / 100000
parachuteRateBase = casesPerWeekBase / 7

set.seed(123)
test <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".test"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 10,                         # Number of trials in simulation
  
  ### population parameters
  N = 24,                                   # Number of classrooms
  S0Pops = 600,                             # Initial number of susceptible - stationary
  S0Popt = 150,                             # Initial number of susceptible - transitory
  
  gradeThreshold = 8,                   # Not treating grades separately
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*600,parachuteRateBase*150),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.2,0,.2,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.05,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.15,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elementary600",
  
  preDetectSuccess = .1,                  # Probability of detecting a pre-symptomatic individual
  sympDetectSuccess = .8,                 # Probability of detecting a symptomatic individual
  postDetectSuccess = .1,                 # Probability of detecting a post-symptomatic individual
  asympDetectSuccess = .1,                # Probability of detecting an asymptomatic individual
)

preK <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".preK"),        # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 1,                                   # Number of classrooms
  S0Pops = 20,                             # Initial number of susceptible - stationary
  S0Popt = 5,                              # Initial number of susceptible - transitory
  
  gradeThreshold = NULL,                   # Not treating grades separately
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*20,parachuteRateBase*5),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = rep(0,4),          # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = rep(0,4),       # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = rep(0,4),       # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = rep(0,4),       # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "preK"
)

elem80FullCohorting <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem80FullCohorting"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 4,                                   # Number of classrooms
  S0Pops = 80,                             # Initial number of susceptible - stationary
  S0Popt = 20,                             # Initial number of susceptible - transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*80,parachuteRateBase*20),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem80FullCohorting",
)

elem80FullMixing <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem80FullMixing"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 4,                                   # Number of classrooms
  S0Pops = 80,                             # Initial number of susceptible - stationary
  S0Popt = 20,                             # Initial number of susceptible - transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*80,parachuteRateBase*20),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.2,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.2,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem80FullMixing",
)

elem80TwiceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem80TwiceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=4),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 8,                                   # Number of classrooms
  S0Pops = 80,                             # Initial number of susceptible - stationary
  S0Popt = 20,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,1),c(0,1),c(0,1),c(0,1),c(3,4),c(3,4),c(3,4),c(3,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*80,parachuteRateBase*20),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem80TwiceAWeek",
)

elem80OnceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem80OnceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=4),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 8,                                   # Number of classrooms
  S0Pops = 80,                             # Initial number of susceptible - stationary
  S0Popt = 20,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,0),c(0,0),c(0,0),c(0,0),c(4,4),c(4,4),c(4,4),c(4,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*80,parachuteRateBase*20),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem80OnceAWeek",
)

elem250FullCohorting <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullCohorting"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
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
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem250FullCohorting",
)

elem250FullMixing <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250FullMixing"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 10,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*250,parachuteRateBase*50),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.2,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.2,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem250FullMixing",
)

elem250TwiceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250TwiceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=10),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 20,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem250TwiceAWeek",
)

elem250OnceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem250OnceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=10),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 20,                                   # Number of classrooms
  S0Pops = 250,                             # Initial number of susceptible - stationary
  S0Popt = 50,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),
                  c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4)),                 # List of days for the weekend
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem250OnceAWeek",
)



elem600FullCohorting <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem600FullCohorting"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 24,                                   # Number of classrooms
  S0Pops = 600,                             # Initial number of susceptible - stationary
  S0Popt = 150,                             # Initial number of susceptible - transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*600,parachuteRateBase*150),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem600FullCohorting",
)

elem600FullMixing <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem600FullMixing"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  N = 24,                                   # Number of classrooms
  S0Pops = 600,                             # Initial number of susceptible - stationary
  S0Popt = 150,                             # Initial number of susceptible - transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*600,parachuteRateBase*150),  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  inGroupTransferRate = c(1/8,0,1/8,0),    # inSchoolstudent, outSchoolstudent, inSchoolteacher, outSchoolteacher; Reciprocal of expected waiting time for in-transfer event
  inGroupTransferNodeNum = c(.25,0,.25,0),   # Average number/proportion of nodes that transfer at each transfer event
  inGroupTransferMinProp = c(.2,0,.3,0),  # Minimum proportion of node population that transfers in each event
  inGroupTransferMaxProp = c(.2,0,.3,0),  # Maximum proportion of node population that transfers in each event
  outGroupTransferRate = rep(0,4),         # Reciprocal of expected waiting time for out-transfer event
  outGroupTransferNodeNum = rep(0,4),      # Average number/proportion of nodes that transfer at each transfer event
  outGroupTransferMinProp = rep(0,4),      # Minimum proportion of node population that transfers in each event
  outGroupTransferMaxProp = rep(0,4),      # Maximum proportion of node population that transfers in each event
  crossGroupTransferRate = rep(0,4),       # Reciprocal of expected waiting time for out-transfer event
  crossGroupTransferNodeNum = rep(0,4),    # Average number/proportion of nodes that transfer at each transfer event
  crossGroupTransferMinProp = rep(0,4),    # Minimum proportion of node population that transfers in each event
  crossGroupTransferMaxProp = rep(0,4),
  
  titleString = "elem600FullMixing",
)

elem600TwiceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem600TwiceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=24),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 48,                                   # Number of classrooms
  S0Pops = 600,                             # Initial number of susceptible - stationary
  S0Popt = 150,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),
                  c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4),c(3,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*600,parachuteRateBase*150),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem600TwiceAWeek",
)

elem600OnceAWeek <- covidSimulator(
  simID = paste0("covid.",Sys.Date(),".elem600OnceAWeek"),      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 100,                         # Number of trials in simulation
  
  ### population parameters
  nodeGroupList = rep(1,2,each=24),                    # List of group IDs for node groups. Default is 1 group of nodes.
  N = 48,                                   # Number of classrooms
  S0Pops = 600,                             # Initial number of susceptible - stationary
  S0Popt = 150,                             # Initial number of susceptible - transitory
  
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),
                  c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4),c(4,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(parachuteRateBase*600,parachuteRateBase*150),  # Reciprocal of expected waiting time in days for a parachute event
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
  crossGroupTransferMaxProp = c(0,0,.3,0),
  
  titleString = "elem600OnceAWeek",
)