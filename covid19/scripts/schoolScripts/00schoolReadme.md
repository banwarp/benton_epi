### Readme file for school scripts
The school scripts are different enough from the other three sets of scripts that I decided to write this supplemental readme. It combines a description of parameters and walkthrough of the main code, pts_function, and event functions.  

The other readmes should be reviewed before reading through this readme, and the scripts themselves are commented in detail.

### Description of model
Instead of a large, generally homogenous population where different nodes may be poorly defined, a school has a relatively small population (<400 for an elementary school, up to 2000 for a high school), in which the nodes can be very defined (e.g. by classroom or by grade), and in which there are at least two distinct types of individual: students and teachers that share a classroom, and teachers that rotate between classrooms. The school scripts can model different types of school structures, for example:
- An elementary school where music, art, and PE teachers rotate between classrooms that otherwise don't mix.
- A middle school with little mixing across grades, but where students in the same grade that don't share a homeroom may encounter each other in electives.
- A high school with nearly homogeneous mixing within a grade and significant mixing across grades, where teachers encounter multiple different groups of students throughout the day.  

In addition to different school structures, I wanted to update the model to reflect a better understanding of symptomatic versus asymptomatic cases. In the school model, all individuals who are infected start out in the preSymptomatic compartment, then some proportion move to the Symptomatic compartment, and the rest move to the asymptomatic compartment. The symptomatic individuals move to the postsymptomatic compartment after a period of time. Each of these four compartments can have different R0s.  

I also wanted to represent different approaches to controlling spread of the coronavirus. For example, an elementary school that is practicing strict cohorting would put each classroom in its own node and each grade in its own node group. The mixing between nodes in a node group could be low, and the mixing across nodegroups could be zero. Furthermore, schools could run half-days, alternating days, or other staggered schedules to reduce onsite populations.

There are also different approaches if a case of the coronavirus is identified in the school. The school could choose to isolate just the positive individual, quarantine the classroom (which is to say, send everyone in the classroom home and do remote learning for that classroom for a number of days), quarantine the grade, or close the school. The user can select which approach to model.

Finally, I wanted to recognize that school days, weekends, and breaks all affect disease dynamics, and transitory individuals may change nodes multiple times a day. Therefore the time step for this model is every hour, and the effective R0 changes after school, on weekends, on breaks, and during quarantine. Also, there can be different start/stop times or different weekdays/weekends for half-day or alternating-day schedules.  

To simplify modeling, I removed parachuting infections and mass-entry. I assume that there are not mass-entry of new individuals into the relatively closed school population, and parachuted infections are replaced by user-defined infection events. The user could define a random set of infection events if desired to model parachuting.

#### Staggered schedules
The model incorporates the daily and weekly school schedule, and allows for different kinds of schedules. Here are the key parameters for setting up the school schedule:
```
dayTimes = list(c(8,16)),                # List of times between which school is in session for the given node
weekDays = list(c(0,4)),                 # List of days for the weekend
breakDays = list(),                      # List of school break days, each element has a start date and end date, can be Date or Numeric
```
`dayTimes` defines the school day in the 24 hour clock. The default is 8am to 4pm. `weekDays` defines the school week in days, starting with Monday = 0. The default is Monday-Friday. `breakDays` defines the dates of school breaks in Dates or Numeric. The default is no breaks.  

The code allows the user to simulate different schedules, like half-day school with half of the student body attending in the morning and half in the afternoon, or alternating days, with half of the students coming Monday, Tuesday, and half coming Wednesday, Thursday. These staggered schedules each have the own start/stop times or weekday/weekend days. The default is to assume the the whole school is on 8am-4pm; Monday-Friday, with no breaks. If there are no staggered schedules, the user can define just one set of parameters for `dayTimes`, `weekDays`, and `breakDays`. If there are staggered schedules, the code needs the schedule for each node. Here is are two alternate examples:  

Suppose there are 18 classrooms in the school:
- If half the students come in the morning each day, and half come in the afternoon, and everyone gets Thanksgiving Week off, there would be 36 nodes, and the parameters would be:
```
dayTimes =  unlist(list(rep(list(c(8,12)),18),rep(list(c(12,16)),18)),recursive=FALSE) # List of times between which school is in session for the given node
weekDays = list(c(0,4)),                          # List of days for the weekend
breakDays = list(c("2020-11-23","2020-11-27"),    # List of school break days, each element has a start date and end date, can be Date or Numeric
```
- If instead half the students come on Monday and Tuesday, and the other half come on Wednesday and Thursday, and everyone gets Thanskgiving Week off, the parameter would be:
```
dayTimes = list(c(8,16)),                         # List of times between which school is in session for the given node
weekDays = unlist(list(rep(list(c(0,1)),18),rep(list(c(2,3)),18)),recursive=FALSE)     # List of days for the weekend
breakDays = list(c("2020-11-23","2020-11-27"),    # List of school break days, each element has a start date and end date, can be Date or Numeric
```

Note: It is possible to have alternating weeks by using `breakDays`. For example, you could define every other week to be a break week, or every other Friday in combination with Monday and Tuesday or Wednesday and Thursday.  

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
Under individual isolation, any individual who is identified (through symptom monitoring or testing, for example) is shifted to an Isolation compartment. For example, if the proportion of symptomatic people who are identified is `sympDetectSuccess`, then at each time step, the transitions would include: `"Symps -> sympDetectSuccess*symptomaticPeriod*Symps -> Isos"` and `"Symps -> sympDetectFailure*symptomaticPeriod*Symps -> postSymps"`.  
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
    preDetectFailure = rep(1-preDetectSuccess,NnumTrials),   # Probability of failing to detect a pre-symptomatic individual
    sympDetectSuccess = rep(sympDetectSuccess,NnumTrials),   # Probability of detecting a symptomatic individual
    sympDetectFailure = rep(1-sympDetectSuccess,NnumTrials), # Probability of failing to detect a symptomatic individual
    postDetectSuccess = rep(postDetectSuccess,NnumTrials),   # Probability of detecting a post-symptomatic individual
    postDetectFailure = rep(1-postDetectSuccess,NnumTrials), # Probability of failing to detect a post-symptomatic individual
    asympDetectSuccess = rep(asympDetectSuccess,NnumTrials),  # Probability of detecting an asymptomatic individual
    asympDetectFailure = rep(1-asympDetectSuccess,NnumTrials) # Probability of failing to detect an asymptomatic individual
  )
```
There are two main differences between the school scripts and the other scripts, both of which are reflected in `v0`.  

First, since the time step is hourly and school is only in session about eight hours a day, five days a week, minus breaks like Winter Break, `schoolDay, schoolWeek, schoolTerm` track these different conditions and change the effective R0 accordingly.  

Second, the school script does not use policy interventions to reduce R0 when prevalence grows. Instead, schools monitor for symptoms and then react by isolating the individual or quarantining the classroom, grade, or school.  `noQuarantine, infectionTimer, quarTimer, quarCounter` are used for the quarantine process.  

`preDetectSuccess, sympDetectSuccess, postDetectSuccess, asympDetectSuccess` determine the success of symptom monitoring and `preDetectFailure, sympDetectFailure, postDetectFailure, asympDetectFailure` determine the failure of symptom monitoring. Theoretically, the success and failure probabilities sum to 1; however, the effect on the compartments of detection versus failure to detect are different: under the individual isolation model, successful detection of an infection results in an immediate removal to the isolation compartment (mathematically, it cancels out the reciprocal of the infectious period), but failure to detect does not immediately move the individual to the next compartment (mathematically, it retains the reciprocal of the infectious period). Therefore it is necessary to have two different parameters for the two different transitions. Note that `preDetectFailure`, etc., are only used in the individual isolation model.  

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
      dTimes = maxDayTimes,                   # number of school day time start/stops, for staggered schedules
      wDays = maxWeekDays,                    # number of school week start/stops, for staggered schedules
      bDays = maxBreakDays,                   # number of school break start/stops, for staggered schedules
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor                     # proportionate reduction in beta at night
    )
```
Every morning at 8am, the post time step function randomly identifies some proportion of infections, within bounds dictated by `preDet (etc.), detVar`. The detected infections will be moved to the isolation compartment at the next time step, while the unidentified infections will undergo normal transitions until the next day at 8am. The post time step function for individual isolation als tracks if the time step is in school or if it is after school, on the weekend, or on a break, for each of the staggered schedules.  

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
      dTimes = maxDayTimes,                   # number of school day time start/stops, for staggered schedules
      wDays = maxWeekDays,                    # number of school week start/stops, for staggered schedules
      bDays = maxBreakDays,                   # number of school break start/stops, for staggered schedules
      noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
      night = nightFactor,                    # proportionate reduction in beta at night
      quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
    )
```
In addition to tracking the time step (in school, after school, weekend, break) for each staggered schedule, `pts_funScriptClassQ` also attempts to detect infectious individuals once each morning in each classroom (including weekends, although I could change this).  

First the function counts the number of preSymptomatic, symptomatic, postSymptomatic, and asymptomatic individuals (combining stationary and transitory). The code generates a uniform random variable between 0 and 1. The detection will be successful if the random variable is above that day's threshold.  

The threshold is computed by the complement of the product of the probability of not detecting any infections in each of the four compartments. I.e. 1-(P(no detection of preSymptomatic)*P(no detection of symptomatic)*P(no detection of postsymptomatic)*P(no detection of asymptomatic))*TimeFactor.  

The Time Factor is included to allow for the probability of detection to depend not just on the number of infections, but also on how long there have been any infections in the classroom. This could reflect, for example, the gradual onset of symptoms or more careful monitoring if students seem fatigued. Time factor defaults to having no effect.

If the random variable is below the probability of detection, then the classroom is put into quarantine (remote learning) for the user-defined number of days. This changes `outOfSchool` to reflect that students/teachers are not mixing in a classroom.  

##### Grade quarantine and school closure post time step functions
These functions are similar to the classroom quarantine, but additional code is included so that a detected infection in one classroom sends not just the node into quarantine, but (for grades) all nodes in the node group or (for school) all nodes in the trial. Due to my lack of understanding of the complexities of the post-time step function, it takes an additional time step to propogate the quarantine from the classroom where the infection was detected to the grade/school, and again when quarantine is lifted. This is not ideal, but since the time step is only 1 hour, the error is negligible.

#### Transfer events
The SimInf package is very good at modeling population movements/contacts between classrooms in a school. To represent contacts between stationary populations, I randomly swap stationary individuals from one node to another. In real life, these individuals would not actually swap classrooms, but it is a good approximation for lower levels of mixing. To represent transitory individuals, I rotate transitory individuals in an orbit of nodes. I defined three types of transfers (in-group, out-group, and cross-group) for both types of populations (stationary and transitory). The transfers events data frame is built in a subroutine `transferFunction` in eventFunctions_school_06.2020.R.

##### Transfers for stationary populations
Stationary populations can swap nodes within a group, to other groups, or across all groups. The process is largely the same. After defining the transfer rate, the code randomly selects nodes. The code selects a proportion of the population in each node, then randomly selects the destination within a node group (in-group), in any node group other than the node's group (out-group), or across all groups (cross-group). To balance populations, the same proportion is transfered from the destination node to the origination node.

##### Transfers for transitory populations
Transitory populations can orbit nodes within a group, to other groups (which is not likely, but included for completeness), or across all groups. For example, an elementary school music teacher would orbit across all node groups (grades), while a high school calculus teacher may orbit among just the senior's node group. After defining the transfer rate, the code generates the complete list of relevant nodes. The destination is always the next node in the group (in-group), in the set of other groups (out-group), or in the trial (cross-group). The proportion of individuals who transfer is user-defined. There is no balancing of populations for transitory transfers.

#### Running the model and plotting the results
The model is built and run the same way. I made some minor updates to the plotting script, specifically making the median trend an option instead of a default and allowing the user to define the name of the data, since with eight infectious compartments the name "preSymps_preSympt_Symps_Sympt_postSymps_postSympt_aSymps_aSympt" is both unwieldy and confusing.
