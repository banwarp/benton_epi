### Readme file for school scripts
The school scripts are different enough from the other three sets of scripts that I decided to write this supplemental readme. It combines a description of parameters and walkthrough of the main code, pts_function, and event functions.  

The other readmes should be reviewed before reading through this readme, and the scripts themselves are commented in detail.

### Description of model
Instead of a large, generally homogenous population where different nodes may be poorly defined, a school has a relatively small population (<400 for an elementary school, up to 2000 for a high school), in which the nodes can be very defined (e.g. by classroom or by grade), and in which there are at least two distinct types of individual, students and teachers that share a classroom, and teachers that rotate between classrooms. The school scripts can model different types of school structures, for example:
- An elementary school where music, art, and PE teachers rotate between classrooms that otherwise don't mix.
- A middle school with little mixing across grades, but where students in the same grade that don't share a homeroom may encounter each other in electives.
- A high school with nearly homogeneous mixing within a grade and significant mixing across grades, where teachers encounter multiple different groups of students throughout the day.  

In addition to different school structures, I wanted to update the model to reflect a better understanding of symptomatic versus asymptomatic cases. In the school model, all individuals who are infected start out in the preSymptomatic compartment, then some proportion move to the Symptomatic compartment, and the rest move to the asymptomatic compartment. The symptomatic individuals move to the postsymptomatic compartment after a period of time. Each of these four compartments can have different R0s.  

I also wanted to represent different approaches to controlling spread of the coronavirus. For example, an elementary school that is practicing strict cohorting would put each classroom in its own node and each grade in its own node group. The mixing between nodes in a node group could be low, and the mixing across nodegroups could be zero.  

There are also different approaches if a case of the coronavirus is identified in the school. The school could choose to isolate just the positive individual, quarantine the classroom (which is to say, send everyone in the classroom home and do remote learning for that classroom for a number of days), quarantine the grade, or close the school. The user selects which approach to model.

Finally, I wanted to recognize that school days, weekends, and breaks all affect disease dynamics, and transitory individuals may change nodes multiple times a day. Therefore the time step for this model is every hour, and the effective R0 changes after school, on weekends, on breaks, and during quarantine.  

To simplify modeling, I removed parachuting infections and mass-entry. I assume that there are not mass-entry of new individuals into the relatively closed school population, and parachuted infections are replaced by user-defined infection events. The user could define a random set of infection events if desired to model parachuting.

#### Compartments
The compartments in the model are as follows:
```
compartments <- c("Ss",        # Susceptible stationary
                    "St",        # Susceptible transitory
                    "preSymps",  # Pre-symptomatic stationary stationary
                    "preSympt",  # Pre-symptomatic transitory
                    "Symps",     # symptomatic Infectious stationary
                    "Sympt",     # symptomatic Infectious transitory
                    "postSymps", # post-symptomatic stationary
                    "postSympt", # post-symptomatic transitory
                    "aSymps",    # asymptomatic Infectious stationary
                    "aSympt",    # asymptomatic Infectious transitory
                    "Rs",        # Recovered stationary
                    "Rt",        # Recovered transitory
                    "Ims",       # Immune stationary
                    "Imt",       # Immune transitory
                    "Isos",      # Isolated stationary
                    "Isot",      # Isolated transitory
                    "cumIs",     # cumulative stationary Infections (used for final tabulation)
                    "cumIt"      # cumulative transitory Infections (used for final tabulation)
                    )
```
The subscripts 's' and 't' refer to stationary and transitory individuals. Students are generally stationary, as are classroom teachers in elementary school. Elective teachers and other rotating staff are transitory.

#### Transitions
In order to save space with complex transitions, I moved the transitions string to a separate R script that is called by the main script. There are two options for the transitions string: Group quarantine (classroom, grade, or school) or Individual isolation. The transitions are slightly different for each.
##### Group quarantine
Under Group quarantine, an identified case will lower the effective R0 to represent that the relevant group has been sent home and is not spending as much time with each other as they would in school. This can be achieved without actually moving the affected group to a different compartment by multiplying the infection transition by the term `outOfSchool`. `outOfSchool` is a continuous variable stored in `v` and updated as necessary (to be explained later):
```
paste0("Ss ->",
           "outOfSchool*",
           "(ssD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "stD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*Ss/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSymps + cumIs"
    )
```

##### Individual isolation
Under individual isolation, any individual who is identified (through symptom monitoring or testing, for example) is shifted to an Isolation compartment. For example, if the proportion of presymptomatic people who are identified is `preDetectSuccess`, then at each time step, the transition would be: `"preSymps -> preDetectSuccess*preSympPeriod*preSymps -> Isos"`  
The infection transtion does not include `Isos` or `Isot` terms to represent that these individuals are isolated.

#### Continuous variables
As in the other scripts, the continuous variables are used to affect transitions and store different states.
```
# Initialized continuous variables
  v0 = data.frame(
    outOfSchool = rep(1,NnumTrials),                         # out-of-school factor for beta
    season = rep(1,NnumTrials),                              # seasonality factor for beta
    symptomaticProp = rep(symptomaticProp,NnumTrials),       # proportion of infections that are symptomatic
    schoolDay = rep(schoolDay,NnumTrials),                   # Indicator for school day: 1 = school day, 0 = evening 
    schoolWeek = rep(schoolWeek,NnumTrials),                 # Indicator for school week: 1 = week, 0 = weekend
    schoolTerm = rep(schoolTerm,NnumTrials),                 # Indicator for school term: 1 = school term, 0 = break
    noQuarantine = rep(noQuarantine,NnumTrials),             # Indicator for classroom level quarantine: 1 = no quarantine, 0 = quarantine
    infectionTimer = rep(0,NnumTrials),                      # Timer for length that at least one infection is present in node
    quarTimer = rep(0,NnumTrials),                           # Timer for length of quarantine
    quarCounter = rep(0,NnumTrials),                         # Counter for number of quarantine events
    preDetectSuccess = rep(preDetectSuccess,NnumTrials),     # Probability of detecting a pre-symptomatic individual
    sympDetectSuccess = rep(sympDetectSuccess,NnumTrials),   # Probability of detecting a symptomatic individual
    postDetectSuccess = rep(postDetectSuccess,NnumTrials),   # Probability of detecting a post-symptomatic individual
    asympDetectSuccess = rep(asympDetectSuccess,NnumTrials)  # Probability of detecting an asymptomatic individual
  )
```
There are two main differences between the school scripts and the other scripts, both of which are reflected in `v0`.  

First, since the time step is hourly and school is only in session about eight hours a day, five days a week, minus breaks like Winter Break, `schoolDay, schoolWeek, schoolTerm` track these different conditions and change the effective R0 accordingly.  

Second, the school script does not use policy interventions to reduce R0 when prevalence grows. Instead, schools monitor for symptoms and then react by isolating the individual or quarantining the classroom, grade, or school.  `preDetectSuccess, sympDetectSuccess, postDetectSuccess, asympDetectSuccess` determine the success of symptom monitoring, and `noQuarantine, infectionTimer, quarTimer, quarCounter` are used for the quarantine process.

#### Post time step function
There are four post time step functions to choose from, Individual isolation, and classroom, grade, or school quarantine.

##### Individual isolation post time step function.
Here are the arguments for `pts_funScriptIndividualQ`:
```
pts_fun <- pts_funScriptIndividualQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      detVar = detectionProbabilityVar,       # Proportionate variable decrease in probability
      bDays = breakDays,                      # List of school breaks, each element has start and end date
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor                     # proportionate reduction in beta at night
    )
```
The post time step function for individual isolation mostly just tracks if the time step is in school or if it is after school, on the weekend, or on a break. It also randomizes the probability of detecting symptoms at each time step, within bounds.

##### Quarantine post time step function (classroom quarantine).
Here are the argumnets for `pts_funScriptClassQ`:
```
pts_fun <- pts_funScriptClassQ(
      startDay = startofSimDay,               # start of simulation
      enn = N,                                # number of nodes in a trial
      numComp = length(compartments),         # number of compartments
      ncolV0 = ncol(v0),                      # number of continuous variables
      cosA = cosAmp,                          # Amplitude of seasonal variation in beta
      sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
      sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
      preDet = preDetectSuccess,              # Probability of detecting a pre-symptomatic individual
      sympDet = sympDetectSuccess,            # Probability of detecting a symptomatic individual
      postDet = postDetectSuccess,            # Probability of detecting a post-symptomatic individual
      asympDet = asympDetectSuccess,          # Probability of detecting an asymptomatic individual
      dsTimeMu = detectSuccessTimeMu,         # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeSig = detectSuccessTimeSig,       # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
      dsTimeLB = detectSuccessLowerBound,     # Lower bound on probability of detecting as time increases
      dsTimeUB = detectSuccessUpperBound,     # Upper bound on probability of detecting as time increases
      qDays = quarantineDays,                 # Length in days of quarantine
      bDays = breakDays,                      # List of school breaks, each element has start and end date
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor,                    # proportionate reduction in beta at night
      quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
    )
```
In addition to tracking the time step (in school, after school, weekend, break), `pts_funScriptClassQ` also attempts to detect infectious individuals once each morning in each classroom (including weekends, although I could change this).  

First the function counts the number of preSymptomatic, symptomatic, postSymptomatic, and asymptomatic individuals (combining stationary and transitory). The code generates a uniform random variable between 0 and 1. The detection will be successful if the random variable is above that day's threshold.  

The threshold is computed by the complement of the product of the probability of not detecting any infections in each of the four compartments. I.e. 1-(P(no detection of preSymptomatic)*P(no detection of symptomatic)*P(no detection of postsymptomatic)*P(no detection of asymptomatic)*TimeFactor.  

The Time Factor is included to allow for the probability of detection to depend not just on the number of infections, but also on how long there have been any infections in the classroom. This could reflect, for example, the gradual onset of symptoms or more careful monitoring if students seem fatigued. Time factor defaults to having no effect.

If the random variable is below the probability of detection, then the classroom is put into quarantine (remote learning) for the user-defined number of days. This changes `outOfSchool` to reflect that students/teachers are not mixing in a classroom.

##### Grade quarantine and school closure post time step functions
These functions are similar to the classroom quarantine, but additional code is included so that a detected infection in one classroom sends not just the node into quarantine, but (for grades) all nodes in the node group or (for school) all nodes in the trial.

#### Transfer events
The SimInf package is very good at modeling population movements/contacts between classrooms in a school. To represent contacts between stationary populations, I randomly swap stationary individuals from one node to another. In real life, these individuals would not actuall swap classroom, but it is a good approximation for lower levels of mixing. To represent transitory individuals, I rotate transitory individuals in an orbit of nodes. I defined three types of transfers (in-group, out-group, and cross-group) for both types of populations (stationary and transitory). The transfers events data frame is built in a subroutine `transferFunction` in eventFunctions_school_06.2020.R.

##### Transfers for stationary populations
Stationary populations can swap nodes within a group, to other groups, or across all groups. The process is largely the same. After defining the transfer rate, the code randomly selects nodes. The code selects a proportion of the population in each node, then randomly selects the destination within a node group (in-group), in any node group other than the node's group (out-group), or across all groups (cross-group). To balance populations, the same proportion is transfered from the destination node to the origination node.

##### Transfers for transitory populations
Transitory populations can orbit nodes within a group, to other groups (which is not likely, but included for completeness), or across all groups. For example, an elementary school music teacher would orbit across all node groups (grades), while a high school calculus teacher may orbit among just the seniors node group. After defining the transfer rate, the code generates the complete list of relevant nodes. The destination is always the next node in the group (in-group), in the set of other groups (out-group), or in the trial (cross-group). The proportion of individuals who transfer is user-defined. There is no balancing of populations for transitory transfers.

#### Running the model and plotting the results
The model is built and run the same way. I made some minor updates to the plotting script, specifically making the median trend an option instead of a default and allowing the user to define the name of the data, since with eight infectious compartments the name "preSymps_preSympt_Symps_Sympt_postSymps_postSympt_aSymps_aSympt" is both unwieldy and confusing.
