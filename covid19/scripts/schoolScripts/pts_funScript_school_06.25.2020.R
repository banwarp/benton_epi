# pts_funScript_schoolIndividualQ_06.30.2020.R

# changes from pts_funScript06.25.2020.R
# added functionality for multiple school times, days, breaks

# changes from pts_funScript06.16.2020.R
# Removed phi and associated code
# Added evening and weekend continuous variables
# Added functionality for quarantining a classroom

# changes from pts_funScript06.09.2020.R
# Added in highSpreadProp randomizer for each time step to represent variability in transmissibility between individuals
# Fixed the defaults in pts_funScript to reflect that parmsList is deprecated.
# added H for calculating community prevalence
# Fixed error in summing compartments for community prevalence - forgot to exclude cumI from denominator

# changes from pts_funScript5.5.2020.R
# Use pI for prevalence, since I and uI are by definition not detectd.

# changes from pts_funScript5.2.2020.R
# revised the neighbor algorithm to be more efficient on the front(R)-end


# changes from pts_funScript5.1.2020.R
# reorganized parameters
#   moved constant parameters from V0 to ldata if they varied across trial
#   moved constant parameters from v0/ldata to gdata/pts_function parameters if they didn't vary across trial
# added Rt tracker
# added in pdDecay - decay parameter to represent the slow decay of physical distancing as people return to normal interactions


# changes from pts_funScript.R
# removed baytat
# increased baseline phi to represent impact of physical distancing and contact tracing
# changed defaults for cosAmp, prevType, upDelay, downDelay, switchOffPolicies, switchOffDay, phiMoveUp, phiMoveDown


########### INDIVIDUAL QUARANTINE ############
########### INDIVIDUAL QUARANTINE ############
########### INDIVIDUAL QUARANTINE ############
########### INDIVIDUAL QUARANTINE ############
########### INDIVIDUAL QUARANTINE ############


pts_funScriptIndividualQ <- function(
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
  dTimes = maxDayTimes,                   # number of school day time start/stops
  wDays = maxWeekDays,                    # number of school week start/stops
  bDays = maxBreakDays,                   # number of school break start/stops
  noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
  night = nightFactor                     # proportionate reduction in beta at night
) {
  pts_fun <- paste0(
    "const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const int Nn = ",
    enn,
    "; /* number of nodes in trial */
    const int Nc = ",
    numComp,
    "; /* number of compartments */
    const int Nv = ",
    ncolV0,
    "; /* number of continuous variables */
    const double cosAmp = ",
    cosA,
    "; /* seasonality factor */
    const double sympProp0 = ",
    sympProp0,
    "; /* Proportion of infections that are symptomatic */
    const double sympVar = ",
    sympVar,
    "; /* Variance of proportion of infections that are symptomatic */
    const double preDet = ",
    preDet,
    "; /* Probability of detecting a pre-symptomatic infection */
    const double sympDet = ",
    sympDet,
    "; /* Probability of detecting a symptomatic infection */
    const double postDet = ",
    postDet,
    "; /* Probability of detecting a post-symptomatic infection */
    const double asympDet = ",
    asympDet,
    "; /* Probability of detecting an asymptomatic infection */
    const double detVar = ",
    detVar,
    "; /* Proportionate variable decrease in probability */
    const int dTimes = ",
    dTimes,
    "; /* number of school day time start/stops */
    const int wDays = ",
    wDays,
    "; /* number of school week start/stops */
    const int bDays = ",
    bDays,
    "; /* number of school break start/stops */
    const double noSchool = ",
    noSchool,
    "; /* Factor to reduce beta when outside of school */
    const double night = ",
    night,
    "; /* Factor to reduce beta at night */
    
    const double pi = 3.141593;        
    
    const int tint = t;
    
    const double day = ceil(t / 24);
    
    const int dayint = day - 1;
    
    // continuous variables
    double outOfSchool = v[0];    // out-of-school factor
    double season = v[1];         // seasonality factor
    double sympProp = v[2];       // proportion of infections that are symptomatic
    double schoolDay = v[3];      // Is it during school hours
    double schoolWeek = v[4];     // Is it during week
    double schoolTerm = v[5];     // Is school in session (or on break)
    double rvDet;
    
    // Starting to update variables
    
    // Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1
    season = cosAmp * cos(pi * (day - 32 + startDay) / 182.5) + 1;
    
    // Generate random variable for proportion of new infectious who are symptomatic
    double rvSymp;
    rvSymp = (double)rand() / RAND_MAX;
    sympProp = sympProp0 + (rvSymp - .5) * 2 * sympVar;
    
    // Attempting to detect infections
    // Always test at 8am, for simplicity
    // Positive detection moves individual to isolation right away, hence dividing by _Period
    // Failure to detect does not move individual at higher rate, hence not dividing by _Period
    if(tint % 24 == 8) {
      rvDet = (double)rand() / RAND_MAX;
      v_new[10] = preDet * (1 + (rvDet - 1) * detVar) / gdata[4];
      v_new[11] = preDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[12] = sympDet * (1 + (rvDet - 1) * detVar) / gdata[5];
      v_new[13] = sympDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[14] = postDet * (1 + (rvDet - 1) * detVar) / gdata[6];
      v_new[15] = postDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[16] = asympDet * (1 + (rvDet - 1) * detVar) / gdata[7];
      v_new[17] = asympDet * (1 + (rvDet - 1) * detVar);
    } else {
      v_new[10] = 0;
      v_new[11] = 1;
      v_new[12] = 0;
      v_new[13] = 1;
      v_new[14] = 0;
      v_new[15] = 1;
      v_new[16] = 0;
      v_new[17] = 1;
    }
    
    // School in session or not
    
    int j; // counter for dayTimes, weekDays, and breakDays
    
    // school day
    // Safety so the first school time isn't overwritten by the test for the second time, etc.
    double outSchool = 0;
    for(j = 0; j < dTimes / 2; j++) {
      if(tint % 24 < ldata[4 + 2 * j] || tint % 24 >= ldata[4 + 2 * j + 1]) {
        schoolDay = 0;
        outSchool = 1;
      } else {
        if(outSchool == 0) {
          schoolDay = 1;
        }
      }
    }
    
    // week day
    // Safety so the first school day isn't overwritten by the test for the second day, etc.
    double weekend = 0;
    for(j = 0; j < wDays / 2; j++){
      if(dayint % 7 < ldata[4 + dTimes + 2 * j] || dayint % 7  > ldata[4 + dTimes + 2 * j + 1]) {
        schoolWeek = 0;
        weekend = 1;
      } else {
      if(weekend == 0) {
          schoolWeek = 1;
        }
      }
    }
    
    // school break
    // Safety so the first school break isn't overwritten by the test for the second break, etc.
    double onBreak = 0;
    for (j = 0; j < bDays / 2; j++) {
      if(day >= ldata[4 + dTimes + wDays + 2 * j] && day < ldata[4 + dTimes + wDays + 2 * j + 1]) {
        schoolTerm = 0;
        onBreak = 1;
      } else {
        if(onBreak == 0) {
          schoolTerm = 1;
        }
      }
    }
    
    if(schoolDay == 0) {
      outOfSchool = night;
    }
    else if(schoolWeek * schoolTerm == 0) {
      outOfSchool = noSchool;
    }
    else {
      outOfSchool = 1;
    }
    
    // Update new continuous variables
    v_new[0] = outOfSchool;
    v_new[1] = season;
    v_new[2] = sympProp;
    v_new[3] = schoolDay;
    v_new[4] = schoolWeek;
    v_new[5] = schoolTerm;
    
    return 0;"
  )
  
  return(pts_fun)
}


########### CLASSROOM QUARANTINE ############
########### CLASSROOM QUARANTINE ############
########### CLASSROOM QUARANTINE ############
########### CLASSROOM QUARANTINE ############
########### CLASSROOM QUARANTINE ############

pts_funScriptClassQ <- function(
  startDay = startofSimDay,               # start of simulation
  enn = N,                                # number of nodes in a trial
  numComp = length(compartments),         # number of compartments
  ncolV0 = ncol(v0),                      # number of continuous variables
  cosA = cosAmp,                          # Amplitude of seasonal variation in beta
  sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
  sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
  preDet = preDetectSuccess,             # Probability of detecting a pre-symptomatic individual
  sympDet = sympDetectSuccess,           # Probability of detecting a symptomatic individual
  postDet = postDetectSuccess,             # Probability of detecting a post-symptomatic individual
  asympDet = asympDetectSuccess,         # Probability of detecting an asymptomatic individual
  dsTimeMu = detectSuccessTimeMu,        # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeSig = detectSuccessTimeSig,      # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeLB = detectSuccessLowerBound,    # Lower bound on probability of detecting as time increases
  dsTimeUB = dsTimeUpperBound,            # Upper bound on probability of detecting as time increases
  qDays = quarantineDays,                 # Length in days of quarantine/isolation
  dTimes = maxDayTimes,                   # number of school day time start/stops
  wDays = maxWeekDays,                    # number of school week start/stops
  bDays = maxBreakDays,                   # number of school break start/stops
  noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
  night = nightFactor,                     # proportionate reduction in beta at night
  quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
) {
  pts_fun <- paste0(
    "const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const int Nn = ",
    enn,
    "; /* number of nodes in trial */
    const int Nc = ",
    numComp,
    "; /* number of compartments */
    const int Nv = ",
    ncolV0,
    "; /* number of continuous variables */
    const double cosAmp = ",
    cosA,
    "; /* seasonality factor */
    const double sympProp0 = ",
    sympProp0,
    "; /* Proportion of infections that are symptomatic */
    const double sympVar = ",
    sympVar,
    "; /* Variance of proportion of infections that are symptomatic */
    const double preDet = ",
    preDet,
    "; /* Probability of detecting a pre-symptomatic infection */
    const double sympDet = ",
    sympDet,
    "; /* Probability of detecting a symptomatic infection */
    const double postDet = ",
    postDet,
    "; /* Probability of detecting a post-symptomatic infection */
    const double asympDet = ",
    asympDet,
    "; /* Probability of detecting an asymptomatic infection */
    const double dsTimeMu = ",
    dsTimeMu,
    "; /* First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeSig = ",
    dsTimeSig,
    "; /* Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeLB = ",
    1 - dsTimeLB,
    "; /* Lower bound on probability of detecting as time increases */
    const double dsTimeUB = ",
    1 - dsTimeUB,
    "; /* Upper bound on probability of detecting as time increases */
    const double qDays = ",
    qDays,
    "; /* Number of days in quarantine */
    const int dTimes = ",
    dTimes,
    "; /* number of school day time start/stops */
    const int wDays = ",
    wDays,
    "; /* number of school week start/stops */
    const int bDays = ",
    bDays,
    "; /* number of school break start/stops */    
    const double noSchool = ",
    noSchool,
    "; /* Factor to reduce beta when outside of school */
    const double night = ",
    night,
    "; /* Factor to reduce beta at night */
    const double quarFactor = ",
    quarFactor,
    "; /* Factor to reduce beta in quarantine */
    
    const double pi = 3.141593;        
    
    const int tint = t;
    
    const double day = ceil(t / 24);
    
    const int dayint = day - 1;
    
    double outOfSchool = v[0];    // out-of-school factor
    double season = v[1];         // seasonality factor
    double sympProp = v[2];       // proportion of infections that are symptomatic
    double schoolDay = v[3];      // Is it during school hours
    double schoolWeek = v[4];     // Is it during week
    double schoolTerm = v[5];     // Is school in session (or on break)
    double noQuarantine = v[6];   // Is there a quarantine for this classroom (1 = no)
    int infectionTimer = v[7]; // Timer for tracking how long classroom has at least one infection
    double quarTimer = v[8];      // Timer for tracking how long quarantine has been going on
    double quarCounter = v[9];    // Counting number of quarantines for this classroom
    
    // Starting to update variables
    
    // Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1
    season = cosAmp * cos(pi * (day - 32 + startDay) / 182.5) + 1;
    
    // Generate random variable for proportion of new infectious who are symptomatic
    double rvSymp;
    rvSymp = (double)rand() / RAND_MAX;
    sympProp = sympProp0 + (rvSymp - .5) * 2 * sympVar;
    
    // Count number of active infections in node
    const double preSymp_i = u[2] + u[3];
    const double Symp_i = u[4] + u[5];
    const double postSymp_i = u[6] + u[7];
    const double aSymp_i = u[8] + u[9];
    double PnoDetectPre, PnoDetectS, PnoDetectPost, PnoDetectA, PTimeFactor, PDetect, rvDetect;
    int infectionTimerInt;
    
    // Update infectionTimer
    if(preSymp_i + Symp_i + postSymp_i + aSymp_i > 0 && quarTimer == 0) {
      infectionTimer++;
      
      infectionTimerInt = (int)infectionTimer;
      if(infectionTimerInt % 24 == 8) {
        
        // Try to detect infections
        PnoDetectPre = pow(1 - preDet, preSymp_i);
        PnoDetectS = pow(1 - sympDet,Symp_i);
        PnoDetectPost = pow(1 - postDet,postSymp_i);
        PnoDetectA = pow(1 - asympDet,aSymp_i);
        PTimeFactor = 1 + (dsTimeLB - dsTimeUB) * (1 / exp(dsTimeMu * (ceil(infectionTimer / 24) - dsTimeSig)) - 1 / exp(dsTimeMu * (- dsTimeSig)));
        PDetect = 1 - PnoDetectPre * PnoDetectS * PnoDetectPost * PnoDetectA * PTimeFactor;
        rvDetect = (double)rand() / RAND_MAX;
        
        if(rvDetect < PDetect) {
          noQuarantine = 0;
          quarTimer = 0;
          quarCounter++;
        }
      }
    } else {
      infectionTimer = 0;
    }
    
    // reset quarantine after qDays days
    if(noQuarantine == 0) {
      quarTimer++;
      if(quarTimer > qDays * 24.0) {
        noQuarantine = 1;
        quarTimer = 0;
      }
    }
    
    // School in session or not
    
    int j; // counter for dayTimes, weekDays, and breakDays
    
    // school day
    // Safety so the first school time isn't overwritten by the test for the second time, etc.
    double outSchool = 0;
    for(j = 0; j < dTimes / 2; j++) {
      if(tint % 24 < ldata[4 + 2 * j] || tint % 24 >= ldata[4 + 2 * j + 1]) {
        schoolDay = 0;
        outSchool = 1;
      } else {
        if(outSchool == 0) {
          schoolDay = 1;
        }
      }
    }
    
    // week day
    // Safety so the first school day isn't overwritten by the test for the second day, etc.
    double weekend = 0;
    for(j = 0; j < wDays / 2; j++){
      if(dayint % 7 < ldata[4 + dTimes + 2 * j] || dayint % 7  > ldata[4 + dTimes + 2 * j + 1]) {
        schoolWeek = 0;
        weekend = 1;
      } else {
      if(weekend == 0) {
          schoolWeek = 1;
        }
      }
    }
    
    // school break
    // Safety so the first school break isn't overwritten by the test for the second break, etc.
    double onBreak = 0;
    for (j = 0; j < bDays / 2; j++) {
      if(day >= ldata[4 + dTimes + wDays + 2 * j] && day < ldata[4 + dTimes + wDays + 2 * j + 1]) {
        schoolTerm = 0;
        onBreak = 1;
      } else {
        if(onBreak == 0) {
          schoolTerm = 1;
        }
      }
    }
    
    if(schoolDay == 0) {
      outOfSchool = night;
    }
    else if(noQuarantine == 0) {
      outOfSchool = quarFactor;
    }
    else if(schoolWeek * schoolTerm == 0) {
      outOfSchool = noSchool;
    }
    else {
      outOfSchool = 1;
    }
    
    // Update new continuous variables
    v_new[0] = outOfSchool;
    v_new[1] = season;
    v_new[2] = sympProp;
    v_new[3] = schoolDay;
    v_new[4] = schoolWeek;
    v_new[5] = schoolTerm;
    v_new[6] = noQuarantine;
    v_new[7] = infectionTimer;
    v_new[8] = quarTimer;
    v_new[9] = quarCounter;
    
    return 0;"
  )
  
  return(pts_fun)
}



########### GRADE QUARANTINE ############
########### GRADE QUARANTINE ############
########### GRADE QUARANTINE ############
########### GRADE QUARANTINE ############
########### GRADE QUARANTINE ############

pts_funScriptGradeQ <- function(
  startDay = startofSimDay,               # start of simulation
  enn = N,                                # number of nodes in a trial
  numComp = length(compartments),         # number of compartments
  ncolV0 = ncol(v0),                      # number of continuous variables
  cosA = cosAmp,                          # Amplitude of seasonal variation in beta
  sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
  sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
  preDet = preDetectSuccess,             # Probability of detecting a pre-symptomatic individual
  sympDet = sympDetectSuccess,           # Probability of detecting a symptomatic individual
  postDet = postDetectSuccess,             # Probability of detecting a post-symptomatic individual
  asympDet = asympDetectSuccess,         # Probability of detecting an asymptomatic individual
  dsTimeMu = detectSuccessTimeMu,        # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeSig = detectSuccessTimeSig,      # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeLB = detectSuccessLowerBound,    # Lower bound on probability of detecting as time increases
  dsTimeUB = dsTimeUpperBound,            # Upper bound on probability of detecting as time increases
  qDays = quarantineDays,                 # Length in days of quarantine/isolation
  dTimes = maxDayTimes,                   # number of school day time start/stops
  wDays = maxWeekDays,                    # number of school week start/stops
  bDays = maxBreakDays,                   # number of school break start/stops
  noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
  night = nightFactor,                     # proportionate reduction in beta at night
  quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
) {
  pts_fun <- paste0(
    "const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const int Nn = ",
    enn,
    "; /* number of nodes in trial */
    const int Nc = ",
    numComp,
    "; /* number of compartments */
    const int Nv = ",
    ncolV0,
    "; /* number of continuous variables */
    const double cosAmp = ",
    cosA,
    "; /* seasonality factor */
    const double sympProp0 = ",
    sympProp0,
    "; /* Proportion of infections that are symptomatic */
    const double sympVar = ",
    sympVar,
    "; /* Variance of proportion of infections that are symptomatic */
    const double preDet = ",
    preDet,
    "; /* Probability of detecting a pre-symptomatic infection */
    const double sympDet = ",
    sympDet,
    "; /* Probability of detecting a symptomatic infection */
    const double postDet = ",
    postDet,
    "; /* Probability of detecting a post-symptomatic infection */
    const double asympDet = ",
    asympDet,
    "; /* Probability of detecting an asymptomatic infection */
    const double dsTimeMu = ",
    dsTimeMu,
    "; /* First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeSig = ",
    dsTimeSig,
    "; /* Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeLB = ",
    1 - dsTimeLB,
    "; /* Lower bound on probability of detecting as time increases */
    const double dsTimeUB = ",
    1 - dsTimeUB,
    "; /* Upper bound on probability of detecting as time increases */
    const double qDays = ",
    qDays,
    "; /* Number of days in quarantine */
    const int dTimes = ",
    dTimes,
    "; /* number of school day time start/stops */
    const int wDays = ",
    wDays,
    "; /* number of school week start/stops */
    const int bDays = ",
    bDays,
    "; /* number of school break start/stops */
    const double noSchool = ",
    noSchool,
    "; /* Factor to reduce beta when outside of school */
    const double night = ",
    night,
    "; /* Factor to reduce beta at night */
    const double quarFactor = ",
    quarFactor,
    "; /* Factor to reduce beta in quarantine */
    
    
    // For quarType 2
    const int *u_0 = &u[-Nc*node]; // Very first node in simulation
    const double *v_0 = &v[-Nv*node]; // Very first entry in v0
    const int firstInGroup = ldata[0]; // first node in node group
    const int numberNodeGroup = ldata[1]; // number of nodes in node group
    
    
    const double pi = 3.141593;        
    
    const int tint = t;
    
    const double day = ceil(t / 24);
    
    const int dayint = day - 1;
    
    double outOfSchool = v[0];    // out-of-school factor
    double season = v[1];         // seasonality factor
    double sympProp = v[2];       // proportion of infections that are symptomatic
    double schoolDay = v[3];      // Is it during school hours
    double schoolWeek = v[4];     // Is it during week
    double schoolTerm = v[5];     // Is school in session (or on break)
    double noQuarantine = v_0[firstInGroup * Nv + 6]; // Is there a quarantine for this grade (1 = no)
    double quarTimer = v_0[firstInGroup * Nv + 8];  // Timer for tracking how long grade has at least one infection
    double quarCounter = v_0[firstInGroup * Nv + 9]; // Timer for tracking how long quarantine has been going on
    int infectionTimer = v_0[firstInGroup * Nv + 7]; // Counting number of quarantines for this grade
    
    // Starting to update variables
    
    // Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1
    season = cosAmp * cos(pi * (day - 32 + startDay) / 182.5) + 1;
    
    // Generate random variable for proportion of new infectious who are symptomatic
    double rvSymp;
    rvSymp = (double)rand() / RAND_MAX;
    sympProp = sympProp0 + (rvSymp - .5) * 2 * sympVar;
    
    
    // School in session or not
    
    int j; // counter for dayTimes, weekDays, and breakDays
    
    // school day
    // Safety so the first school time isn't overwritten by the test for the second time, etc.
    double outSchool = 0;
    for(j = 0; j < dTimes / 2; j++) {
      if(tint % 24 < ldata[4 + 2 * j] || tint % 24 >= ldata[4 + 2 * j + 1]) {
        schoolDay = 0;
        outSchool = 1;
      } else {
        if(outSchool == 0) {
          schoolDay = 1;
        }
      }
    }
    
    // week day
    // Safety so the first school day isn't overwritten by the test for the second day, etc.
    double weekend = 0;
    for(j = 0; j < wDays / 2; j++){
      if(dayint % 7 < ldata[4 + dTimes + 2 * j] || dayint % 7  > ldata[4 + dTimes + 2 * j + 1]) {
        schoolWeek = 0;
        weekend = 1;
      } else {
      if(weekend == 0) {
          schoolWeek = 1;
        }
      }
    }
    
    // school break
    // Safety so the first school break isn't overwritten by the test for the second break, etc.
    double onBreak = 0;
    for (j = 0; j < bDays / 2; j++) {
      if(day >= ldata[4 + dTimes + wDays + 2 * j] && day < ldata[4 + dTimes + wDays + 2 * j + 1]) {
        schoolTerm = 0;
        onBreak = 1;
      } else {
        if(onBreak == 0) {
          schoolTerm = 1;
        }
      }
    }
    
    // Grade quarantine
    
    // Count number of active infections in grade
    
    double preSymp_i = 0, Symp_i = 0, postSymp_i = 0, aSymp_i = 0;
        int m;
        for(m = firstInGroup; m < firstInGroup + numberNodeGroup; m++) {
          preSymp_i += u_0[m * Nc + 2] + u_0[m * Nc + 3];
          Symp_i += u_0[m * Nc + 4] + u_0[m * Nc + 5];
          postSymp_i += u_0[m * Nc + 6] + u_0[m * Nc + 7];
          aSymp_i = u_0[m * Nc + 8] + u_0[m * Nc + 9];
        }
    
    double PnoDetectPre, PnoDetectS, PnoDetectPost, PnoDetectA, PTimeFactor, PDetect = 0, rvDetect = 0;
        int infectionTimerInt;
  
      // Test all students in a grade, but only once per day, instead of testing all students in a grade for each class  
      
      // Update infectionTimer
      if(node == firstInGroup) {
        if(preSymp_i + Symp_i + postSymp_i + aSymp_i > 0 && quarTimer == 0) {
          infectionTimer++;
        
          infectionTimerInt = (int)infectionTimer;
          if(infectionTimerInt % 24 == 8) {
            
            // Try to detect infections
            PnoDetectPre = pow(1 - preDet, preSymp_i);
            PnoDetectS = pow(1 - sympDet,Symp_i);
            PnoDetectPost = pow(1 - postDet,postSymp_i);
            PnoDetectA = pow(1 - asympDet,aSymp_i);
            PTimeFactor = 1 + (dsTimeLB - dsTimeUB) * (1 / exp(dsTimeMu * (ceil(infectionTimer / 24) - dsTimeSig)) - 1 / exp(dsTimeMu * (- dsTimeSig)));
            PDetect = 1 - (PnoDetectPre * PnoDetectS * PnoDetectPost * PnoDetectA * PTimeFactor);
            rvDetect = (double)rand() / RAND_MAX;
          } 
          
          if(rvDetect < PDetect) {
            noQuarantine = 0;  // noQuarantine
            quarTimer = 0;  // quarTimer
            quarCounter++; // quarCounter
          }
        } // End of active infection, no quarantine logic
        // If no active infections or if quarantine
        else {
            infectionTimer = 0;  // infectionTimer
        }
      
      // reset quarantine after qDays days
      if(noQuarantine == 0) {
        quarTimer++;
        if(quarTimer > qDays * 24.0) {
          noQuarantine = 1;
          quarTimer = 0;
        }
      }    
    } // End of firstInGroup logic
    
  // Effect on beta
    if(schoolDay == 0) {
      outOfSchool = night;
    }
    else if(noQuarantine == 0) {
      outOfSchool = quarFactor;
    }
    else if(schoolWeek * schoolTerm == 0) {
      outOfSchool = noSchool;
    }
    else {
      outOfSchool = 1;
    }
    
    // Update new continuous variables
    v_new[0] = outOfSchool;
    v_new[1] = season;
    v_new[2] = sympProp;
    v_new[3] = schoolDay;
    v_new[4] = schoolWeek;
    v_new[5] = schoolTerm;
    v_new[6] = noQuarantine;
    v_new[7] = infectionTimer;
    v_new[8] = quarTimer;
    v_new[9] = quarCounter;
    
    return 0;"
  )
  
  return(pts_fun)
}



########### SCHOOL CLOSURE ############
########### SCHOOL CLOSURE ############
########### SCHOOL CLOSURE ############
########### SCHOOL CLOSURE ############
########### SCHOOL CLOSURE ############



pts_funScriptSchoolQ <- function(
  startDay = startofSimDay,               # start of simulation
  enn = N,                                # number of nodes in a trial
  numComp = length(compartments),         # number of compartments
  ncolV0 = ncol(v0),                      # number of continuous variables
  cosA = cosAmp,                          # Amplitude of seasonal variation in beta
  sympProp0 = symptomaticProp,            # Baseline high spreader proportion that gets re-randomized each time step
  sympVar = symptomaticVariance,          # Uniform variance of high spreader proportion
  preDet = preDetectSuccess,             # Probability of detecting a pre-symptomatic individual
  sympDet = sympDetectSuccess,           # Probability of detecting a symptomatic individual
  postDet = postDetectSuccess,             # Probability of detecting a post-symptomatic individual
  asympDet = asympDetectSuccess,         # Probability of detecting an asymptomatic individual
  dsTimeMu = detectSuccessTimeMu,        # First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeSig = detectSuccessTimeSig,      # Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform
  dsTimeLB = detectSuccessLowerBound,    # Lower bound on probability of detecting as time increases
  dsTimeUB = dsTimeUpperBound,            # Upper bound on probability of detecting as time increases
  qDays = quarantineDays,                 # Length in days of quarantine/isolation
  dTimes = maxDayTimes,                   # number of school day time start/stops
  wDays = maxWeekDays,                    # number of school week start/stops
  bDays = maxBreakDays,                   # number of school break start/stops
  noSchool = noSchoolFactor,              # proportionate reduction in beta out of school
  night = nightFactor,                     # proportionate reduction in beta at night
  quarFactor = quarantineFactor           # proportionate reduction in beta due to quarantine
) {
  pts_fun <- paste0(
    "const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const int Nn = ",
    enn,
    "; /* number of nodes in trial */
    const int Nc = ",
    numComp,
    "; /* number of compartments */
    const int Nv = ",
    ncolV0,
    "; /* number of continuous variables */
    const double cosAmp = ",
    cosA,
    "; /* seasonality factor */
    const double sympProp0 = ",
    sympProp0,
    "; /* Proportion of infections that are symptomatic */
    const double sympVar = ",
    sympVar,
    "; /* Variance of proportion of infections that are symptomatic */
    const double preDet = ",
    preDet,
    "; /* Probability of detecting a pre-symptomatic infection */
    const double sympDet = ",
    sympDet,
    "; /* Probability of detecting a symptomatic infection */
    const double postDet = ",
    postDet,
    "; /* Probability of detecting a post-symptomatic infection */
    const double asympDet = ",
    asympDet,
    "; /* Probability of detecting an asymptomatic infection */
    const double dsTimeMu = ",
    dsTimeMu,
    "; /* First parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeSig = ",
    dsTimeSig,
    "; /* Second parameter for function to increase probability of detecting as time goes on; can be linear, sigmoid, or uniform */
    const double dsTimeLB = ",
    1 - dsTimeLB,
    "; /* Lower bound on probability of detecting as time increases */
    const double dsTimeUB = ",
    1 - dsTimeUB,
    "; /* Upper bound on probability of detecting as time increases */
    const double qDays = ",
    qDays,
    "; /* Number of days in quarantine */
    const int dTimes = ",
    dTimes,
    "; /* number of school day time start/stops */
    const int wDays = ",
    wDays,
    "; /* number of school week start/stops */
    const int bDays = ",
    bDays,
    "; /* number of school break start/stops */
    const double noSchool = ",
    noSchool,
    "; /* Factor to reduce beta when outside of school */
    const double night = ",
    night,
    "; /* Factor to reduce beta at night */
    const double quarFactor = ",
    quarFactor,
    "; /* Factor to reduce beta in quarantine */
    
    
    // For quarType 2
    const int *u_0 = &u[-Nc*node]; // Very first node in simulation
    const double *v_0 = &v[-Nv*node]; // Very first entry in v0
    const int firstInTrial = ldata[2]; // first node in node group
    
    const double pi = 3.141593;        
    
    const int tint = t;
    
    const double day = ceil(t / 24);
    
    const int dayint = day - 1;
    
    double outOfSchool = v[0];    // out-of-school factor
    double season = v[1];         // seasonality factor
    double sympProp = v[2];       // proportion of infections that are symptomatic
    double schoolDay = v[3];      // Is it during school hours
    double schoolWeek = v[4];     // Is it during week
    double schoolTerm = v[5];     // Is school in session (or on break)
    double noQuarantine = v_0[firstInTrial * Nv + 6]; // Is there a quarantine for this grade (1 = no)
    double quarTimer = v_0[firstInTrial * Nv + 8];  // Timer for tracking how long grade has at least one infection
    double quarCounter = v_0[firstInTrial * Nv + 9]; // Timer for tracking how long quarantine has been going on
    int infectionTimer = v_0[firstInTrial * Nv + 7]; // Counting number of quarantines for this grade
    
    // Starting to update variables
    
    // Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1
    season = cosAmp * cos(pi * (day - 32 + startDay) / 182.5) + 1;
    
    // Generate random variable for proportion of new infectious who are symptomatic
    double rvSymp;
    rvSymp = (double)rand() / RAND_MAX;
    sympProp = sympProp0 + (rvSymp - .5) * 2 * sympVar;
    
    
    // School in session or not
    
    int j; // counter for dayTimes, weekDays, and breakDays
    
    // school day
    // Safety so the first school time isn't overwritten by the test for the second time, etc.
    double outSchool = 0;
    for(j = 0; j < dTimes / 2; j++) {
      if(tint % 24 < ldata[4 + 2 * j] || tint % 24 >= ldata[4 + 2 * j + 1]) {
        schoolDay = 0;
        outSchool = 1;
      } else {
        if(outSchool == 0) {
          schoolDay = 1;
        }
      }
    }
    
    // week day
    // Safety so the first school day isn't overwritten by the test for the second day, etc.
    double weekend = 0;
    for(j = 0; j < wDays / 2; j++){
      if(dayint % 7 < ldata[4 + dTimes + 2 * j] || dayint % 7  > ldata[4 + dTimes + 2 * j + 1]) {
        schoolWeek = 0;
        weekend = 1;
      } else {
      if(weekend == 0) {
          schoolWeek = 1;
        }
      }
    }
    
    // school break
    // Safety so the first school break isn't overwritten by the test for the second break, etc.
    double onBreak = 0;
    for (j = 0; j < bDays / 2; j++) {
      if(day >= ldata[4 + dTimes + wDays + 2 * j] && day < ldata[4 + dTimes + wDays + 2 * j + 1]) {
        schoolTerm = 0;
        onBreak = 1;
      } else {
        if(onBreak == 0) {
          schoolTerm = 1;
        }
      }
    }
    
    // Grade quarantine
    
    // Count number of active infections in school
    
    double preSymp_i = 0, Symp_i = 0, postSymp_i = 0, aSymp_i = 0;
        int m;
        for(m = firstInTrial; m < firstInTrial + Nn; m++) {
          preSymp_i += u_0[m * Nc + 2] + u_0[m * Nc + 3];
          Symp_i += u_0[m * Nc + 4] + u_0[m * Nc + 5];
          postSymp_i += u_0[m * Nc + 6] + u_0[m * Nc + 7];
          aSymp_i = u_0[m * Nc + 8] + u_0[m * Nc + 9];
        }
    
    double PnoDetectPre, PnoDetectS, PnoDetectPost, PnoDetectA, PTimeFactor, PDetect = 0, rvDetect = 0;
        int infectionTimerInt;
  
      // Test all students in a grade, but only once per day, instead of testing all students in a grade for each class  
      
      // Update infectionTimer
      if(node == firstInTrial) {
        if(preSymp_i + Symp_i + postSymp_i + aSymp_i > 0 && quarTimer == 0) {
          infectionTimer++;
        
          infectionTimerInt = (int)infectionTimer;
          if(infectionTimerInt % 24 == 8) {
            
            // Try to detect infections
            PnoDetectPre = pow(1 - preDet, preSymp_i);
            PnoDetectS = pow(1 - sympDet,Symp_i);
            PnoDetectPost = pow(1 - postDet,postSymp_i);
            PnoDetectA = pow(1 - asympDet,aSymp_i);
            PTimeFactor = 1 + (dsTimeLB - dsTimeUB) * (1 / exp(dsTimeMu * (ceil(infectionTimer / 24) - dsTimeSig)) - 1 / exp(dsTimeMu * (- dsTimeSig)));
            PDetect = 1 - (PnoDetectPre * PnoDetectS * PnoDetectPost * PnoDetectA * PTimeFactor);
            rvDetect = (double)rand() / RAND_MAX;
          } 
          
          if(rvDetect < PDetect) {
            noQuarantine = 0;  // noQuarantine
            quarTimer = 0;  // quarTimer
            quarCounter++; // quarCounter
          }
        } // End of active infection, no quarantine logic
        // If no active infections or if quarantine
        else {
            infectionTimer = 0;  // infectionTimer
        }
      
      // reset quarantine after qDays days
      if(noQuarantine == 0) {
        quarTimer++;
        if(quarTimer > qDays * 24.0) {
          noQuarantine = 1;
          quarTimer = 0;
        }
      }    
    } // End of firstInTrial logic
    
  // Effect on beta
    if(schoolDay == 0) {
      outOfSchool = night;
    }
    else if(noQuarantine == 0) {
      outOfSchool = quarFactor;
    }
    else if(schoolWeek * schoolTerm == 0) {
      outOfSchool = noSchool;
    }
    else {
      outOfSchool = 1;
    }
    
    // Update new continuous variables
    v_new[0] = outOfSchool;
    v_new[1] = season;
    v_new[2] = sympProp;
    v_new[3] = schoolDay;
    v_new[4] = schoolWeek;
    v_new[5] = schoolTerm;
    v_new[6] = noQuarantine;
    v_new[7] = infectionTimer;
    v_new[8] = quarTimer;
    v_new[9] = quarCounter;
    
    return 0;"
  )
  
  return(pts_fun)
}
