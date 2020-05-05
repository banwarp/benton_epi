# pts_funScript5.2.2020.R

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

pts_funScript <- function(
  phiPhysicalDistancing = (gdata$beta/gdata$gamma)/parmList$RBaseline, # Ongoing baseline Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  phiNoAction = (gdata$beta/gdata$gamma)/parmList$RNoActions, # Ongoing Rt, reflecting no actions
  maxPrev1 = parmList$maxPrev1,                              # Maximum prevalence before instituting minor intervention
  maxPrev2 = parmList$maxPrev2,                              # Maximum prevalence before instituting major intervention
  phiFactor1 = (gdata$beta/gdata$gamma)/parmList$RTarget1,   # Target for the reduction in R under minor intervention
  phiFactor2 = (gdata$beta/gdata$gamma)/parmList$RTarget2,   # Target for the reduction in R under major intervention
  cosAmp = parmList$cosAmp,                                  # Amplitude of seasonal variation in beta
  startDay = startofSimDay,                                  # start of simulation
  kbDay1 = as.numeric(as.Date(parmList$kbDay1)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of first phase
  kbDay2 = as.numeric(as.Date(parmList$kbDay2)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of second phase
  kbDay3 = as.numeric(as.Date(parmList$kbDay3)) - as.numeric(as.Date("2020-01-01")) - startofSimDay, # date of third phase
  parmCol = ncol(lparms),                                    # number of columns in lparms
  numComp = length(compartments),                            # number of compartments
  prevType = ifelse(parmList$maxPrev1 < 1,1,0),              # prevalence type. 0 = count, 1 = proportion
  upDelay = parmList$upDelay,                                # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = parmList$downDelay,                            # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = parmList$phiMoveUp,                            # Rate at which phi increases when interventions are imposed
  phiMoveDown = parmList$phiMoveDown,                        # Rate at which phi decreases when interventions are lifted
  pdDecay = parmList$pdDecay,                                # Rate at which phi decreases toward 1 in the absence of interventions. Represents gradual relaxation of physical distancing
  switchOffPolicies = parmList$switchOffPolicies,            # Logical variable to switch off policies (used for counterfactuals)
  switchOffDay = as.numeric(as.Date(parmList$switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
) {
  pts_fun <- paste0(
    "const double phiBase = ",
    phiPhysicalDistancing,
    "; /* Ongoing baseline Rt, reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders */
    const double phiNoAction = ",
    phiNoAction,
    "; /* Ongoing baseline Rt, under no preventive actions */
    const double maxPrev1 = ",
    maxPrev1,
    "; /* Maximum prevalence before instituting minor intervention */
    const double maxPrev2 = ",
    maxPrev2,
    "; /* Maximum prevalence before instituting major intervention */
    const double phiFactor1 = ",
    phiFactor1,
    "; /* Target for the reduction in R under minor intervention */
    const double phiFactor2 = ",
    phiFactor2,
    "; /* Target for the reduction in R under major intervention */
    const double cosAmp = ",
    cosAmp,
    "; /* amplitude of cosine wave that gives seasonality to beta */
    const double startDay = ",
    startDay,
    "; /* Start date of simulation for seasonality calculation */
    const double kbDay1 = ",
    kbDay1,
    "; /* Day that current stay-at-home is lifted - phase 1 */
    const double kbDay2 = ",
    kbDay2,
    "; /* Day that current stay-at-home is lifted  - phase 2 */
    const double kbDay3 = ",
    kbDay3,
    "; /* Day that current stay-at-home is lifted - phase 3 */
    const double prevType = ",
    prevType,
    "; /* type of prevalence: 0 = count, 1 = proportion */
    const double upDelay = ",
    upDelay,
    "; /* Delay before starting response to increased prevalence */
    const double downDelay = ",
    downDelay,
    "; /* Delay before stepping down response */
    const double phiMoveUp = ",
    phiMoveUp,
    ";    /* rate at which phi converges up to phiFactor; [0,1] bigger is faster convergence */
    const double phiMoveDown = ",
    phiMoveDown,
    ";    /* rate at which phi converges down to phiFactor; [0,1] bigger is faster convergence */
    const double pdDecay = ",
    pdDecay,
    ";    /* rate at which phi decays to 1 as people gradually stop physical distancing */
    const double switchOffPolicies = ",
    switchOffPolicies,
    "; /* Logical if policies will be switched off. 1 if they will be switched off; 0 if the won't be switched off */
    const double switchOffDay = ",
    switchOffDay,
    ";     /* Day that policies will be switched off */
    
    const double RNaught = gdata[0]; // Basic reproduction number
    const double betaRandomizer = ldata[0]; // Trial randomizer for beta or R0
    
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
    
    const double phi = v[0];
    const double prevUp01 = v[1];    /* variable capturing delay for moving up from lowest prev (0) to intermediate prev (1) */
    const double prevUp12 = v[2];    /* variable capturing delay for moving up from intermediate prev (1) to highest prev (2) */
    const double prevDown21 = v[3];  /* variable capturing delay for moving down from highest prev (2) to intermediate prev (1) */
    const double prevDown10 = v[4];  /* variable capturing delay for moving down from intermediate prev (1) to lowest prev (0) */
    const double kbPhase_0 = v[5];  /* logical: 0 means current stay-at-home, 1-3 means stay-at-home lifted */
    const double policy_0 = v[6];    /* logical: 1 means policies are used, 0 means not */
    const double prev_0 = v[8];      /* prevalence in previous time step */
    const double prevState_0 = v[9]; /* intervention state before minor intervention */
    const double Rtee_0 = v[10];     /* current effective Rt */
    const double pdCounter_0 = v[12]; // Counter for pdDecay
    
    // update variables in case they arent updated later
    v_new[1] = prevUp01;
    v_new[2] = prevUp12;
    v_new[3] = prevDown21;
    v_new[4] = prevDown10;
    v_new[5] = kbPhase_0;
    v_new[6] = policy_0;
    v_new[8] = prev_0;
    v_new[9] = prevState_0;
    v_new[10] = Rtee_0;
    v_new[12] = pdCounter_0;
    
    /* Update seasonal effect: one cosine wave over course of year, peak on Feb 1, trough on Aug 1 */
    v_new[7] = cosAmp * cos(pi * (t - 32 + startDay) / 180) + 1;
    
    /* ***** Computing prevalence for each trial ***** */
    
    /* Initalize the compartment state vector to find number of individuals in neighbors */
    const int *u_0 = &u[-Nc*node];
    
    int j, k;
    double N_j = N_i, I_j = I_i, prev, kbPhase = kbPhase_0, policy=policy_0, pdCounter = pdCounter_0, pdB;
    
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
    
    /* Phase out current stay-at-home */
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
    double phiBaseline = phiBase;
    double phiTarget1 = phiFactor1;
    double phiTarget2 = phiFactor2;
    
    if(pdDecay > 0 && kbPhase == 3) {
      // If prevalence has been below maxPrev1 for long enough, start pdDecay
      if(prev <= maxPrev1) {
        // if downDelay has been completed
        if(prevDown10 > downDelay)
          v_new[12] = pdCounter_0 + 1; // Increases pdCounter by 1
          if(v_new[12] < pdDecay) {
            // decrease phiBaseline toward 1
            pdB = log(-v_new[12] + (pdDecay + 1)) /  log(pdDecay + 1);
          } else {
            // maintain phiBaseline at phiNoAction
            pdB = 0;
          }
          phiBaseline = pdB * phiBaseline + (1-pdB) * phiNoAction;
      }
    } else {pdCounter = 0;}
    
    // overrides lifting stay-at-home in phases, also overrides pdDecay
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
    
    // Dynamic changes in phi based on prevalence
    
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
          v_new[12] = 0; // Reset pdCounter
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
          v_new[12] = 0; // Reset pdCounter
        }
      // end of policy == 1 clause
      }
      else {
        v_new[0] = phi * (1-phiMoveDown) + phiBaseline * phiMoveDown;
        v_new[4] = prevDown10 + 1;
      }
    
    v_new[10] = RNaught * betaRandomizer * v_new[7] / v_new[0];
    
    return 0;"
  )
  
  return(pts_fun)
}