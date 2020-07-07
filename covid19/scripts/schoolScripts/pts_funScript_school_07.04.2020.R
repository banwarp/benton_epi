# pts_funScript_school_07.04.2020.R

# changes from pts_funScript_school_07.02.2020.R
# Changed quarantine decision process from 1 case leads to quarantine to user-defined number of cases leads to quarantine
# Removed false positives - too complex for right now

# changes from pts_funScript06.30.2020.R
# Removed Infection Timer
# Added false positives

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
#   moved constant parameters from v0/ldata to gdata/pts_function parameters if they didnt vary across trial
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


pts_funScriptSchool <- function(
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
  detVar = detectionProbabilityVar,       # Proportionate variable decrease in probability
  cThreshold = classroomThreshold,
  gThreshold = gradeThreshold,
  sThreshold = schoolThreshold,
  cClassrooms = countClassrooms,
  qDaysC = quarantineDaysClassroom,                 # Length in days of quarantine/isolation
  qDaysG = quarantineDaysGrade,                 # Length in days of quarantine/isolation
  qDaysS = quarantineDaysSchool,                 # Length in days of quarantine/isolation
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
    const double detVar = ",
    detVar,
    "; /* Proportionate variable decrease in probability */
    const double cThreshold = ",
    cThreshold,
    "; /* Number of cases before classroom quarantined */
    const double gThreshold = ",
    gThreshold,
    "; /* Number of cases before grade quarantined */
    const double sThreshold = ",
    sThreshold,
    "; /* Number of cases before school quarantined */
    const double cClassrooms = ",
    cClassrooms,
    "; /* Whether to count classrooms that are on quarantine when counting active detected cases; 1 = TRUE */
    const double qDaysC = ",
    qDaysC,
    "; /* Number of days in quarantine */
    const double qDaysG = ",
    qDaysG,
    "; /* Number of days in quarantine */
    const double qDaysS = ",
    qDaysS,
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
    
    const int *u_0 = &u[-Nc*node]; // Very first node in simulation
    const double *v_0 = &v[-Nv*node]; // Very first entry in v0
    const int firstInGroup = ldata[0]; // first node in node group
    const int numberNodeGroup = ldata[1]; // number of nodes in node group
    const int firstInTrial = ldata[2]; // first node in node group
    
    double outOfSchool = v[0];    // out-of-school factor
    double season = v[1];         // seasonality factor
    double sympProp = v[2];       // proportion of infections that are symptomatic
    double schoolDay = v[3];      // Is it during school hours
    double schoolWeek = v[4];     // Is it during week
    double schoolTerm = v[5];     // Is school in session (or on break)
    
    double noQuarC = v[10];   // Is there a quarantine for this classroom (1 = no)
    double qTimerC = v[11];      // Timer for tracking how long quarantine has been going on
    double qCounterC = v[12];    // Counting number of quarantines for this classroom

    double noQuarG = v_0[firstInGroup * Nv + 13]; // Is there a quarantine for this grade (1 = no)
    double qTimerG = v_0[firstInGroup * Nv + 14];  // Timer for tracking length of grade quarantine
    double qCounterG = v_0[firstInGroup * Nv + 15]; // Counter for number of grade quarantines
    
    double noQuarS = v_0[firstInTrial * Nv + 16]; // Is there a quarantine for school (1 = no)
    double qTimerS = v_0[firstInTrial * Nv + 17];  // Timer for tracking length of school closure
    double qCounterS = v_0[firstInTrial * Nv + 18]; // Counter for number of school closures
    
    // Updating variables
    
    // Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1
    season = cosAmp * cos(pi * (day - 32 + startDay) / 182.5) + 1;
    
    // Generate random variable for proportion of new infectious who are symptomatic
    double rvSymp;
    rvSymp = (double)rand() / RAND_MAX;
    sympProp = sympProp0 + (rvSymp - .5) * 2 * sympVar;
    
    // update quarantine variables because logic is not working
    v_new[10] = noQuarC;
    v_new[11] = qTimerC;
    v_new[12] = qCounterC;
    v_new[13] = noQuarG;
    v_new[14] = qTimerG;
    v_new[15] = qCounterG;
    v_new[16] = noQuarS;
    v_new[17] = qTimerS;
    v_new[18] = qCounterS;
    
    // Attempting to detect infections
    // Always test at 8am, for simplicity
    // Positive detection moves individual to isolation right away, hence dividing by _Period
    // Failure to detect does not move individual at higher rate, hence not dividing by _Period
    // modulo 7 to set up next time step
    
    double rvDet;
    if(tint % 24 == 7) {
      rvDet = (double)rand() / RAND_MAX;
      v_new[6] = preDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[7] = sympDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[8] = postDet * (1 + (rvDet - 1) * detVar);
      
      rvDet = (double)rand() / RAND_MAX;
      v_new[9] = asympDet * (1 + (rvDet - 1) * detVar);
    } else {
      v_new[6] = 0;
      v_new[7] = 0;
      v_new[8] = 0;
      v_new[9] = 0;
    }
    
    // Determine if quarantine is warranted
    
    // initialize counts by classroom, grade, school, also random variable for false positive
    double isoC = 0, isoG = 0, isoS = 0, rvFPos;
    int m;

    // Number of active detected cases in classroom
    isoC = u[14] + u[15];
  
    // Number of active detected cases in grade
    if(node == firstInGroup) {
      for(m = firstInGroup; m < firstInGroup + numberNodeGroup; m++) {
        // Option to only count classrooms not on classroom quarantine
        if(cClassrooms == 1 || v[m * Nv + 10] == 1) {
          isoG += u_0[m * Nc + 14] + u_0[m * Nc + 15];
        }
      }
    }
  
    // Number of active detected cases in school
    if(node == firstInTrial) {
      for(m = firstInTrial; m < firstInTrial + Nn; m++) {
        // Option to only count classrooms not on classroom quarantine or grade quarantine
        if(cClassrooms == 1 || v[m * Nv + 10] * v[m * Nv + 13] == 1) {      
          isoS += u_0[m * Nc + 14] + u_0[m * Nc + 15];
        }
      }
    }

    if(tint % 24 == 8 && (qTimerC + qTimerG + qTimerS) == 0) {
    
      // Start quarantine(s) if warranted
      if(isoC >= cThreshold){
        noQuarC = 0;  // noQuarantine - classroom
        qTimerC = 0;  // quarTimer - classroom
        qCounterC++; // quarCounter - classroom
      }
      if(isoG >= gThreshold) {
        noQuarG = 0;  // noQuarantine - grade
        qTimerG = 0;  // quarTimer - grade
        qCounterG++; // quarCounter - grade
      }
      if(isoS >= sThreshold) {
        noQuarS = 0;  // noQuarantine - school
        qTimerS = 0;  // quarTimer - school
        qCounterS++; // quarCounter - school
      }
    }
    
    // End classroom quarantine after qDaysC
    if(noQuarC == 0) {
      qTimerC++;
      if(qTimerC >= qDaysC * 24.0) {
        noQuarC = 1;
        qTimerC = 0;
      }
    }
    
    // End grade quarantine after qDaysG
    if(noQuarG == 0) {
      qTimerG++;
      if(qTimerG >= qDaysG * 24.0) {
        noQuarG = 1;
        qTimerG = 0;
      }
    }
    
    // End school quarantine after qDaysS
    if(noQuarS == 0) {
      qTimerS++;
      if(qTimerS >= qDaysS * 24.0) {
        noQuarS = 1;
        qTimerS = 0;
      }
    }
    
    // School in session or not
    
    int j; // counter for dayTimes, weekDays, and breakDays
    
    // school day
    // Safety so the first school time isnt overwritten by the test for the second time, etc.
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
    // Safety so the first school day isnt overwritten by the test for the second day, etc.
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
    // Safety so the first school break isnt overwritten by the test for the second break, etc.
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
  
  
    // Apply reduction(s) to beta
    if(schoolDay == 0) {
      outOfSchool = night;
    }
    else if(noQuarC == 0 || noQuarG == 0 || noQuarS == 0) {
      outOfSchool = quarFactor;
    }
    else if(schoolWeek * schoolTerm == 0) {
      outOfSchool = noSchool;
    }
    else {
      outOfSchool = 1;
    }
    
    // Update other new continuous variables
    v_new[0] = outOfSchool;
    v_new[1] = season;
    v_new[2] = sympProp;
    v_new[3] = schoolDay;
    v_new[4] = schoolWeek;
    v_new[5] = schoolTerm;
    v_new[10] = noQuarC;
    v_new[11] = qTimerC;
    v_new[12] = qCounterC;
    v_new[13] = noQuarG;
    v_new[14] = qTimerG;
    v_new[15] = qCounterG;
    v_new[16] = noQuarS;
    v_new[17] = qTimerS;
    v_new[18] = qCounterS;
    
    return 0;"
  )
  
  return(pts_fun)
}