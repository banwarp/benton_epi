### List of parameters with defaults for BentonCountyCOVID19Simulator_Education_08042020.R

```
  ### Simulation parameters
  folderPath = NULL,                       # folder path for subroutines and output
  simID = "covid.mo.day.2020.School",      # Simulation ID
  simDate = "2020-09-01",                  # Start date of simulation
  maxT = 100,                              # number of days in simulation
  numTrials = 50,                          # Number of trials in simulation
  
  ### population parameters
  N = 18,                                  # Number of classrooms
  nodeGroupList = NULL,                    # List of group IDs for node groups. Default is 1 group of nodes.
  Symp0nodeGroups = NULL,                  # Node groups where initial infectious are distributed
  Trans0nodes = NULL,                      # Nodes group where transitory individuals are distributed
  S0Pops = 25*18,                          # Initial number of susceptible - stationary; total per trial
  preSymp0Pops = 0,                        # Initial number of presymptomatic - stationary; count per node
  Symp0Pops = 0,                           # Initial number of symptomatic and asymptomatic - stationary; count per node
  postSymp0Pops = 0,                       # Initial number of post-symptomatic - stationary; count per node
  R0Pops = 0,                              # Initial number of recovered - stationary; count per node
  Im0Pops = 0,                             # Initial number of immune - stationary; count per node
  Iso0Pops = 0,                            # Initial number of isolated - stationary; count per nodey
  S0Popt = 10,                             # Initial number of susceptible - transitory; total per trial
  preSymp0Popt = 0,                        # Initial number of pre-symptomatic - transitory; count per node
  Symp0Popt = 0,                           # Initial number of symptomatic and asymptomatic - transitory; count per node
  postSymp0Popt = 0,                       # Initial number of post-symptomatic - transitory; count per node
  R0Popt = 0,                              # Initial number of recovered - transitory; count per node
  Im0Popt = 0,                             # Initial number of immune - transitory; count per node
  Iso0Popt = 0,                            # Initial number of isolated - transitory; count per node
  
  ### gdata parameters
  R0Base = 3,                              # Baseline R0; compromize betwen 2.5 from CDC and 3.28 from OHA
  R0preFrac = .4,                          # Presymptomatic R0 fraction (what proportion of new cases come from presymptomatic)
  R0sympFrac = .55,                        # Symptomatic R0 fraction
  R0postFrac = .05,                        # Postsymptomatic R0 fraction
  R0asympFrac = .4,                        # Asymptomatic R0 fraction
  R0ComplianceFrac = 1,                    # Proportionate reduction in R0 with high compliance to masks/distancing
  R0VentilationFrac = 1,                   # Proportionate reduction in R0 with good ventilation
  studentTeacherDiff = c(.33,1,.5,1),      # Differential R0 for stationary-stationary, transitory-transitory, and stationary-transitory interactions
  preSympPeriod = 1/2,                     # Reciprocal of presymptomatic period
  symptomaticPeriod = 1/4,                 # Reciprocal of symptomatic period
  postSympPeriod = 1/6,                    # Reciprocal of post-symptomatic period
  aSympPeriod = 1/6,                       # Reciprocal of asymptomatic period
  isoPeriod = 1/10,                        # Reciprocal of isolation period
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = 1/100,                   # Reciprocal of temporary immunity period, after which R becomes Im or S
  R0Spread = 0,                            # Uniform variation in R0 across trials (measured as %change from base R0)
  
  ### continuous initialized parameters (v0)
  symptomaticProp = c(.43,.65),            # Average proportion of infectious who are symptomatic; (stationary, transitory)
  symptomaticVariance = 0,                 # Uniform variance of symptomatic proportion
  schoolDay = 1,                           # Indicator for school day: 1 = school day, 0 = evening, list for each node or single entry
  schoolWeek = 1,                          # Indicator for school week: 1 = week, 0 = weekend, list for each node or single entry
  schoolTerm = 1,                          # Indicator for school term: 1 = school term, 0 = break, list for each node or single entry
  noQuarantine = 1,                        # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
  dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
  weekDays = list(c(0,4)),                 # List of days for the weekend
  breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
  
  ### pts_fun parameters
  cosAmp = 0.125,                          # Amplitude of seasonal variation in beta
  quarantineDaysClassroom = 14,            # Length in days of quarantine
  quarantineDaysGrade = 14,                # Length in days of quarantine
  quarantineDaysSchool = 28,               # Length in days of quarantine
  preDetectSuccess = 0.1,                  # Probability of detecting a pre-symptomatic individual
  sympDetectSuccess = 0.8,                 # Probability of detecting a symptomatic individual
  postDetectSuccess = 0.1,                 # Probability of detecting a post-symptomatic individual
  asympDetectSuccess = 0.1,                # Probability of detecting an asymptomatic individual
  detectionProbabilityVar = 0,             # Proportionate variable decrease in detection probability - only decreases probability, never increases
  classroomThreshold = 2,                  # Number of active detected cases before classroom quarantined
  gradeThreshold = 8,                      # Number of active detected cases before grade quarantined
  schoolThreshold = 16,                    # Number of active detected cases before school quarantined
  countClassrooms = TRUE,                  # Whether to count classrooms on quarantine when counting active detected cases
  noSchoolFactor = .2,                     # Factor to reduce beta when outside of school
  nightFactor = .05,                       # Transmission rate at night; incorporates afternoon also
  quarantineFactor = .1,                   # Factor to reduce beta when quarantine
  
  ### *day* specified infection events parameters - list to allow multiple events, use list structure to accommodate stationary/transitory
  ### take care with infection of transitory - the populations are so small that it risks trying to shift 1+ from a pop of 0
  ### might need to revisit this, perhaps have students, teachers, and transitory as three different compartments
  infectionEventsDF = NULL,                # pre-built infection spreader events data frame if the events will be preset before the script is run
  infectionCount = c(),                    # Number of infections caused by the infection spreader
  infectionNodes = c(),                    # Number of nodes that the infection spreader contacts
  infectionNodeGroups = NULL,              # Which node groups the infection spreader contacts. Must use list() syntax for multiple events
  infectionDate = c(),                     # Date the infection spreader lands. Date can also be numeric i.e. 200
  infectionSpread = c(),                   # Symmetric spread in days of infection spreader infections
  infectionCohort = c(),                   # Infection cohort: 1=stationary or 2=transitory
  
  ### *day* random infection events "parachuters" parameters; list of two for (stationary,transitory)
  paraEventsDF = NULL,                     # pre-built parachuter events data frame if the events will be preset before the script is run
  paraMu = c(1,1),                         # First shape parameter for beta function for timing of parachute events
  paraSig = c(1,1),                        # Second shape parameter for beta function for timing of parachute events
  parachuteRate = c(0,0),                  # Reciprocal of expected waiting time in days for a parachute event
  parachuteNum = c(1,1),                   # Number of infectious in each parachute event
  parachuteNodeGroups = NULL,              # Which node groups the parachuters can land in
  
  ### *hourly* transfer event parameters
  transferEventsDF = NULL,                 # pre-built *hourly* transfer events data frame if the events will be preset before the script is run
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
  crossGroupTransferMaxProp = rep(0,4),    # Maximum proportion of node population that transfers in each event
    
  ### plot parameters
  plotCompList = "cumI",                   # List of compartments that will be plotted
  rollM = 1,                               # Number of days for plotting rolling means, rollM = 1 means no rolling mean
  allTraj = TRUE,                          # Logical if all simulation trajectories are plotted or just median and spread
  plotRandomTrajs = 5,                     # Number of randomly selected trajectories to plot, used when allTraj == TRUE
  percentile = .5,                         # Percentile of point-in-time simulations to plot, default to median
  confIntv = .95,                          # two-sided confidence interval for plotting spread
  plotGroups = NULL,                       # Which node groups to plot. NULL plots all
  sumGroups = TRUE,                        # Whether to sum the node groups.
  dateBreaks = "1 month",                  # Plot parameter, x-axis displays months
  dataTitle = "Infections",                # Title of data that is plotted
  titleString = "Generic Title",           # Title of plot
  xString = "Date",                        # Title of x-axis
  yString = "Frequency",                   # Title of y-axis
  lString = "Median",                      # Title of legend
  cString = NULL                           # Plot caption
  ```
