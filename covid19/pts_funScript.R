# pts_funScript.R

pts_funScript <- function(
  cosAmp = .25,                   # amplitude of cosine wave for seasonality of beta
  startDay = startofSimDay,       # start day of simulation - number of days after 1/1/YEAR
  kbDay1 = kbDate1,               # day current physical distancing is lifted - phase 1
  kbDay2 = kbDate2,               # day current physical distancing is lifted - phase 2
  kbDay3 = kbDate3,               # day current physical distancing is lifted - phase 3
  parmCol = ncol(parms),          # number of columns of local parameters
  numComp = length(compartments), # number of compartments
  prevType = 0,                   # type of prevalence: 0 = count, 1 = proportion
  upDelay = 14,                   # number of days that beta-reduction measures are delayed after prevalence passes threshold
  downDelay = 28,                 # number of days that beta-reduction measures are relaxed after prevalence decreases past threshold
  switchOffPolicies = 0,          # logical variable to switch off policies (used for counterfactuals)
  switchOffDay = 90,              # day policies would be switched off (used for counterfactuals)
  phiMoveUp = .5,                 # rate at which phi converges up to phiFactor; [0,1] bigger is faster convergence
  phiMoveDown = .5                # rate at which phi converges down to phiFactor; [0,1] bigger is faster convergence
) {
  pts_fun <- paste0(
    "const double cosAmp = ",
    cosAmp,
    "; /* amplitude of cosine wave that gives seasonality to beta */
    const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const double kbDay1 = ",
    kbDay1,
    "; /* Day that current physical distancing is lifted - phase 1 */
    const double kbDay2 = ",
    kbDay2,
    "; /* Day that current physical distancing is lifted  - phase 2 */
    const double kbDay3 = ",
    kbDay3,
    "; /* Day that current physical distancing is lifted - phase 3 */
    const double prevType = ",
    prevType,
    "; /* type of prevalence: 0 = count, 1 = proportion */
    const double upDelay = ",
    upDelay,
    "; /* Delay before starting response to increased prevalence */
    const double downDelay = ",
    downDelay,
    "; /* Delay before stepping down response */
    const double switchOffPolicies = ",
    switchOffPolicies,
    "; /* Logical if policies will be switched off. 1 if they will be switched off; 0 if the won't be switched off */
    const double switchOffDay = ",
    switchOffDay,
    ";     /* Day that policies will be switched off */
    const double phiMoveUp = ",
    phiMoveUp,
    ";    /* rate at which phi converges up to phiFactor; [0,1] bigger is faster convergence */
    const double phiMoveDown = ",
    phiMoveDown,
    ";    /* rate at which phi converges down to phiFactor; [0,1] bigger is faster convergence */
    
    /* Initialize neighbors at the correct row. */
    const double *neighbors = &ldata[",
    parmCol,
    "];
    /* number of compartments: S E I R Im M. Include M for code purposes */
    const int Nc = ",
    numComp,
    "; 

    const double pi = 3.141593;        
    
    const double I_i = u[2];

    const double N_i = u[0] + u[1] + u[3] + u[4] + I_i;
    
    const double maxPrev1 = ldata[0];
    const double maxPrev2 = ldata[1];
    const double phiFactor1 = ldata[2];
    const double phiFactor2 = ldata[3];
    
    const double phi = v[0];
    const double prevUp01 = v[1];    /* variable capturing delay for moving up from lowest prev (0) to intermediate prev (1) */
    const double prevUp12 = v[2];    /* variable capturing delay for moving up from intermediate prev (1) to highest prev (2) */
    const double prevDown21 = v[3];  /* variable capturing delay for moving down from highest prev (2) to intermediate prev (1) */
    const double prevDown10 = v[4];  /* variable capturing delay for moving down from intermediate prev (1) to lowest prev (0) */
    const double kbPhase_0 = v[5];  /* logical: 1 means current physical distancing, 0 means physical distancing lifted */
    const double policy_0 = v[6];    /* logical: 1 means policies are used, 0 means not */
    const double prev_0 = v[8];      /* prevalence in previous time step */
    const double prevState_0 = v[9]; /* intervention state before minor intervention */
    
    // update variables in case they arent updated later
    v_new[1] = prevUp01;
    v_new[2] = prevUp12;
    v_new[3] = prevDown21;
    v_new[4] = prevDown10;
    v_new[5] = kbPhase_0;
    v_new[6] = policy_0;
    v_new[8] = prev_0;
    v_new[9] = prevState_0;
    
    /* Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1 */
    v_new[7] = cosAmp * cos(pi * (t - 32 + startDay) / 180) + 1;
    
    /* ***** Computing prevalence for each trial ***** */
    
    /* Initalize the compartment state vector to find number of individuals in neighbors */
    const int *u_0 = &u[-Nc*node];
    
    int j, k;
    double N_j = N_i, I_j = I_i, prev, kbPhase = kbPhase_0, policy=policy_0;
    
    j = (int)*neighbors++;
    while(j >= 0) {
      
      /* Add infected in node j to total */
      I_j += u_0[j * Nc + 2];
      
      /* Add pop in node j, compartment k to total */
      for (k = j * Nc; k < (j + 1) * Nc; k++)
        N_j += u_0[k];
      
      N_j -= u_0[((j + 1) * Nc) - 1]; /* subtracts M from total pop */
        
      neighbors++;
      
      j = (int)*neighbors++;
    
    }
    
    /* Compute prevalence across all nodes in trial */
    if(prevType == 0) {
     prev = I_j;
    } else {
     prev = I_j/N_j;
    }
    v_new[8] = prev; /* recording prevalence */
    
    
    
    /****************************************************/
    /* update phi with delays if policies are in effect */
    /****************************************************/
    
    /* Phase out current physical distancing */
    if(t == kbDay1) {
      kbPhase = 1;
      v_new[5] = kbPhase;
    } else if(t == kbDay2) {
      kbPhase = 2;
      v_new[5] = kbPhase;
    } else if(t == kbDay3) {
      kbPhase = 3;
      v_new[5] = kbPhase;
    }
    
    /* Switch off policies if indicated */
    if(switchOffPolicies == 1) {
      if(t == switchOffDay) {
      policy = 0;
      v_new[6] = policy;
      }
    }
    
    
    /* ***** Updating phi ***** */
    
    /* temporary variables storing phi targets */
    double phiBaseline = 1;
    double phiTarget1 = phiFactor1;
    double phiTarget2 = phiFactor2;
    
    // overrides lifting physical distancing in phases
    if(kbPhase_0 == 0) {
      phiBaseline = phiFactor2;
      phiTarget1 = phiFactor2;
    }
    else if(kbPhase_0 == 1) {
      phiBaseline = (phiFactor1+phiFactor2) / 2;
      phiTarget1 = (phiFactor1+phiFactor2) / 2;
    }
    else if(kbPhase_0 == 2) {
      phiBaseline = phiFactor1;
    }
    
    if(policy == 1) {
      
        /* ***** baseline phi - no intervention ***** */
        if(prev <= maxPrev1) {
          // past state A: major intervention
          if(prev_0 > maxPrev2) {
            v_new[0] = phi * (1 - phiMoveUp) + phiTarget2 * phiMoveUp; // Maintain major intervention
            v_new[3] = 1; // Start prevDown21 delay
            v_new[4] = 0; // Wait before starting delay of prevDown10
            // reset upDelays
            v_new[1] = 0;
            v_new[2] = 0;
          }
          // past state B: minor intervention
          else if(prev_0 > maxPrev2) {
            // Logical: already doing minor intervention
            if(prevDown21 > downDelay) {
              v_new[0] = phi * (1-phiMoveDown) + phiTarget1 * phiMoveDown; // Transitioning toward no intervention
              v_new[4] = 1; // start prevDown10 delay
            }
            // still doing major intervention
            else {
              v_new[0] = phi * (1-phiMoveUp) + phiTarget2 * phiMoveUp; // Maintaining major intervention
              v_new[3] = prevDown21 + 1;  // Increase step toward completing major intervention
              v_new[4] = 0; // Wait to start delay of prevDown10
            }
            // reset upDelays
            v_new[1] = 0;
            v_new[2] = 0;
          }
          // past state C: baseline
          else {
            // Logical: already doing minor intervention and completed downDelay
            if(prevDown10 > downDelay) {
              v_new[0] = phi * (1-phiMoveDown) + phiBaseline * phiMoveDown; // Transition toward no interventions
              v_new[9] = 0; // Setting previous state to baseline
            }
            // doing minor intervention but downDelay not complete
            else if(prevDown21 > downDelay) {
              v_new[0] = phi * (1-phiMoveDown) + phiTarget1 * phiMoveDown; // Transition toward minor intervention
              v_new[4] = prevDown10 + 1; // Increase step toward completing low downDelay
            }
            // still doing major intervention
            else {
              v_new[0] = phi * (1-phiMoveUp) + phiTarget2 * phiMoveUp; // Maintain major intervention
              v_new[3] = prevDown21 + 1; // Increase step toward completing high downDelay
              v_new[4] = 0; // Wait to start delay of prevDown10
            }
          }
        }
        
        /* ***** minor intervention ***** */
        else if(prev <= maxPrev2) {
          // past state A: baseline
          if(prev_0 <= maxPrev1) {
            if(prevState_0 == 0) {
              v_new[0] = phi * (1-phiMoveUp) + phiBaseline * phiMoveUp; // Maintaining no intervention
              v_new[1] = 1; // starting delay to minor intervention
              // resetting prevUp12
              v_new[2] = 0;
            }
            // // temporary wobble below low threshold while in minor intervention
            else if(prevState_0 == 1) {
              v_new[0] = phi * (1-phiMoveDown) + phiTarget1 * phiMoveDown; // Maintaining minor intervention
            }
            // temporary wobble below low threshold while in process of decreasing from major intervention
            else if(prevState_0 == 2) {
              v_new[0] = phi * (1-phiMoveDown) + phiTarget2 * phiMoveDown; // Maintaining major intervention
              v_new[3] = prevDown21 + 1; // Increase step toward completing high downDelay
            }
          }
          // past state B: major intervention
          else if(prev_0 > maxPrev2) {
              v_new[0] = phi * (1-phiMoveDown) + phiTarget2 * phiMoveDown; // Maintain major intervention
              v_new[3] = 1; // Start delay toward completing major intervention
              v_new[4] = 0; // reset downDelay10
          }
          // past state C: minor intervention
          else {
            // Moving up from baseline
            if(prevState_0 == 0) {
              if(prevUp01 > upDelay) {
                v_new[0] = phi * (1-phiMoveUp) + phiTarget1 * phiMoveUp; // Transition toward minor intervention
                v_new[9] = 1; // Setting previous state to minor intervention
              } else {
                v_new[0] = phi * (1-phiMoveUp) + phiBaseline * phiMoveUp; // Maintain no intervention
                v_new[1] = prevUp01 + 1; // Increase step toward minor intervention
              }
            } else {
              // Logical: already doing major intervention and completed downDelay
              if(prevDown21 > downDelay) {
                v_new[0] = phi * (1-phiMoveDown) + phiTarget1 * phiMoveDown; // Transition toward minor intervention
                v_new[9] = 1; // Setting previous state to minor intervention
              } else {
                v_new[0] = phi * (1-phiMoveUp) + phiTarget2 * phiMoveUp; // Maintain major intervention
                v_new[3] = prevDown21 + 1; // Increase step toward completing high downDelay
              }
              // resetting prevUps
              v_new[1] = 0;
              v_new[2] = 0;
            }
          }
        }
        
        /* ***** major intervention ***** */
        else {
          // past state A: baseline
          if(prev_0 <= maxPrev1) {
            v_new[0] = phi * (1-phiMoveUp) + phiBaseline * phiMoveUp; // Maintaining no intervention
            v_new[2] = 1; // Starting delay for major intervention (note: skips minor intervention)
            // Reset downDelays
            v_new[3] = 0;
            v_new[4] = 0;
          }
          // past state B: minor intervention
          else if(prev_0 <= maxPrev2) {
            if(prevUp01 > upDelay) {
              v_new[0] = phi * (1-phiMoveUp) + phiTarget1 * phiMoveUp; // Maintaining minor intervention
              v_new[1] = prevUp01 + 1;
            } else {
              v_new[0] = phi * (1-phiMoveUp) + phiBaseline * phiMoveUp; // Maintaining no intervention
            }
            v_new[2] = 1; // Starting delay for major intervention (note: may skip part of minor intervention)
            // resetting down delays
            v_new[3] = 0;
            v_new[4] = 0;
          }
          // past state C: major intervention
          else {
            if(prevUp12 > upDelay) {
              v_new[0] = phi * (1-phiMoveUp) + phiTarget2 * phiMoveUp; // Transitioning toward major intervention
              v_new[9] = 2; // Setting previous state to major intervention 
            } else {
              v_new[2] = prevUp12 + 1; // Increasing step toward start of major intervention
              if(prevUp01 > upDelay) {
                v_new[0] = phi * (1-phiMoveUp) + phiTarget1 * phiMoveUp; // Already doing minor intervention
              } else {
                v_new[0] = phi * (1-phiMoveUp) + phiBaseline * phiMoveUp; // Still in state of no intervention
                v_new[1] = prevUp01 + 1;
              }
            }
          }
        }
      // end of policy == 1 clause
      }
      else {
        v_new[0] = phi * (1-phiMoveDown) + phiBaseline * phiMoveDown;
      }
    
    v_new[11] = 2.25 / (12 * v_new[0]); // DELETE AFTER TROUBLESHOOTING //
    
    return 0;"
  )
  
  return(pts_fun)
}