## Post Time Step Function Readme
#### Introduction
This set of R scripts was developed for local modeling of the COVID-19 disease in Benton County, Oregon, and other small geographies. The scripts use a modified SEIR model with partially-mixed, but otherwise homogeneous population nodes. This markdown document describes the Post Time Step Function used in the policyScripts (and also the AgeGroup and Hilow scripts) as of 06/26/2020.

##### Acknowledgements:
I'd like to thank the developers of the SimInf package for writing the code that these scripts use and for assisting me in adapting their code to meet the requirements of COVID-19 modeling.  

Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12  

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales. International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

[SimInf vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf)  
[SimInf technical documentation](https://cran.r-project.org/web/packages/SimInf/SimInf.pdf)  
[SimInf git repository](https://github.com/stewid/SimInf)

##### Conceptual schematic of SimInf
![Conceptual schematic of SimInf](images/simInfSchematic.png)
The conceptual schematic of the SimInf model that I use has three stages: Building the model, Running the model, and the Result.
##### Building the model:
- Specify the model with parameters for disease spread, population compartments, and transitions between compartments.
- Set the initial states in the compartments and the initial values of the continuous variables.
- (If the model used is not one of the built-in model) Write the post-time-step function that changes continuous variables after every time step.
- Generate a dataframe of events that can introduce or remove individuals to different compartments, transfer them between nodes, or shift them between compartments, at different times.
##### Run the model:
- Advance individuals through compartments using the transitions and disease parameters. Record the compartment counts and continous variable values for the disease trajectory.
- Use the events dataframe to move individuals in/out/between nodes and/or compartments.
- Advance the time step.
- Change continuous variables according to the post-time-step functions.
- Repeat for every time step.
##### Result:
- Trajectories of the compartments across the timespan.
- Evolution of continous variables across the timespan.
- With the result, you can plot trajectories for the different compartments and/or the evolution of the continuous variables.

### Post Time Step Function
The post time step function runs after each time step and updates the continuous variables stored in the data frame `v`. These continuous variables, which are initialized in `v0`, are used in the compartment transitions along with `gdata` and `ldata`, but unlike `gdata` and `ldata`, they can vary across time, either exogeneously to the disease conditions (such as a seasonality factor) or in response to disease conditions (such as a policy intervention). In addition to having a direct impact on transitions, `v` stores variables that are used in the post time step function, while changing them as needed, and for counters/delay trackers that change over time. The post time step function is written and implemented in C.  

The main purpose of the post time step function is to calculate the observed prevalence across all nodes in a trial, then use that observed prevalence to decide whether to impose or lift a policy intervention that acts by reducing the effective R0 in the transitions. This represents, for example, an imperfectly mixed county (i.e. multiple nodes), imposing a county-wide policy even if the infections are concentrated in portion of the county.

#### The components of the Post Time Step Function
- `v` is the data frame that stores continuous variables
- pts_function arguments are user-defined constant parameters used to building the post time step function
- `pts_funScript` is a C-script written in a .R file that is read as a string to the main function
- `pts_fun` is the string that contains the post-time-step function and is used in building and running the SimInf model

#### Initializing v0
The continuous variables are initialized in `v0`:
```
v0 = data.frame(
  phi = rep(phi0,NnumTrials),                # initial beta reduction factor (larger phi = more reduction)
  prevUp01 = rep(upDelay+1,NnumTrials),      # initialized variable for capturing delay in response when prevalence increases
  prevUp12 = rep(upDelay+1,NnumTrials),      # variable for capturing delay in response when prevalence increases
  prevDown21 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  prevDown10 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
  kbPhase = rep(kbPhase,NnumTrials),         # 0 means current physical distancing, 1-3 is phase of lifting.
  policy = rep(1,NnumTrials),                # logical: 1 means policies can take effect; 0 means they don't
  season = rep(1,NnumTrials),                # seasonality factor for beta
  observedPrev = rep(0,NnumTrials),          # prevalence tracker - observed from pI. Does not include I or uI by definition.
  previousState = rep(2,NnumTrials),         # previous state: 0 = baseline, 1 = low intervention, 2 = high intervention, used for low intervention logic
  pdCounter = rep(0,NnumTrials)              # counter for pdDecay
)
```

##### Explanation of each continuous variable
`phi` is the key variable that influences the rate of infection. Here is the transition from Susceptible to Exposed:
```
# Explanation of exposure rate:
# "S -> (betaI*(1/phi)*I+                                             # Exposures: Exposure rate among initial infectious
#        betaP*pI+                                                    # Exposure rate among post infectious; does not depend on phi
#        betaU*(1/phi)*uI)*                                           # Exposure rate among unknown infectious
#        season*betaRandomizer*                                       # Seasonal effect and randomizer
#        S/(S+E+I+uI+pI+R+Im)-> E",                                   # Mixing effect
```
The larger that `phi` is, the lower the infection rate. `phi` is updated in the post time step function based on the observed prevalence of COVID-19. For example, suppose the prevalence increases above a certain threshold. A user-defined delay is observed, simulating the lag time between an uptick in prevalence and awareness of that uptick. After the delay is completed, `phi` converges to its new value at an exponential rate. `phi` remains there until prevalence changes again, at which point another delay is observed, and the `phi` converges to its new value.

`prevUp01, prevUp12, prevDown21, prevDown10` are the delay trackers used when changing `phi`. `prevUp01` is set to 0 whenever the prevalence increases from its lowest level (level 0) up to the middle level (level 1). Then, each time step increases `prevUp01` until `prevUp01 = upDelay`, at which point `phi` begins to change. `prevUp12, prevDown21, prevDown10` are similar.  

`kbPhase` is a stored variable representing the initial state of reopening following guidance from the governor (Kate Brown). `kbPhase = 0` initializes the simulation in the state of full "stay-home, stay healthy". `kbPhase = 1` is phase 1 of reopening, etc. The phases act by inflating `phi` above what the observed prevalence would call for. `kbPhase` increases on user-defined dates to represent the phased reopening.  

`policy` is a stored variable that either allows policies to be used (`policy = 1`), or prevents policies from being used (`policy = 0`). `policy` can be changed from 1 to 0 on a user-defined date.  

`season` is the seasonal variation in transmissibility.  

`observedPrev` stores the prevalence from the previous time step. This is used when deciding how to change `phi`. For example, if the previous prevalence was below threshold 01, and the current prevalence is above the threshold, initialize `prevUp01 = 0`.  

`previousState` stores the state of policy intervention from the previous time step. It is used to smooth transitions from one state to another. For example, it is possible that the prevalence temporarily rises above threshold 01, but then falls back below the threshold. Using `previousState = 0` prevents superfluous re-initializations and movements in the other variables.  

`pdCounter` is a counter for the decay in physical distancing. The model allows for the assumption that in the absence of policy interventions, people's adherence to good physical distancing decays over a period of time as they "forget" about COVID-19. `pdCounter` controls how long it takes for adherence to physical distancing to decay to a user-defined baseline.  

#### Defining arguments for pts_function
The arguments for `pts_function` are user-defined in the main script. They are used as constant parameters when building `pts_fun`. Some parameters affect `phi`, others affect the policy interventions, and others are used for calculating prevalence:
```
pts_fun <- pts_funScript(
  phiPhysicalDistancing = (R0I+R0U)/RPhysicalDistancing, # Phi reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
  phiNoAction = (R0I+R0U)/RNoAction,                     # Phi reflecting no actions
  maxPrev1 = maxPrev1,                                   # Maximum prevalence before instituting minor intervention
  maxPrev2 = maxPrev2,                                   # Maximum prevalence before instituting major intervention
  phiFactor1 = (R0I+R0U)/RTarget1,                       # Target for the reduction in R0 under minor intervention
  phiFactor2 = (R0I+R0U)/RTarget2,                       # Target for the reduction in R0 under major intervention
  cosAmp = cosAmp,                                       # Amplitude of seasonal variation in beta
  startDay = startofSimDay,                              # start of simulation
  kbDay1 = kbDay1,                                       # date of first phase
  kbDay2 = kbDay2,                                       # date of second phase
  kbDay3 = kbDay3,                                       # date of third phase
  enn = N,                                               # number of nodes in a trial
  numComp = length(compartments),                        # number of compartments
  pICompt = which(compartments=="pI")-1,                 # which compartment is pI for prevalence monitoring
  prevType = ifelse(maxPrev1 < 1,1,0),                   # prevalence type. 0 = count, 1 = proportion
  upDelay = upDelay,                                     # Number of days after prevalence passes threshold until minor/major intervention
  downDelay = downDelay,                                 # Number of days after prevalence drops below threshold until intervention lifted
  phiMoveUp = phiMoveUp,                                 # Rate at which phi increases when interventions are imposed
  phiMoveDown = phiMoveDown,                             # Rate at which phi decreases when interventions are lifted
  pdDecay = pdDecay,                                     # Rate at which phi decreases toward phiNoAction in the absence of interventions. Represents gradual relaxation of physical distancing
  switchOffPolicies = switchOffPolicies,                 # Logical variable to switch off policies (used for counterfactuals)
  switchOffDay = as.numeric(as.Date(switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
)
```

##### Explanation of each pts_function argument
'phiPhysicalDistancing' is the level of `phi` that reflects no policy interventions, but does include individual physical distancing. The model allows for the assumption that individuals will maintain individual physical distancing, and that this physical distancing will reduce the effective R0. `phi` decays from `phiPhysicalDistancing` according to the `pdDecay` factor (physical distancing decay) and `pdCounter` from `v`.  

`phiNoAction` is the level of `phi` that reflects no actions (i.e. no policy interventions and no physical distancing). `phiNoAction` is the target to which `phiPhysicalDistancing` decays.  

`maxPrev1, maxPrev2` are the user-defined maximum acceptable prevalence levels before instituting moderate and major policy interventions (subject to delays), and also the prevalence levels at which interventions are lifted (again subject to delays). `maxPrev1, maxPrev2` can be counts or proportions; `prevType` is the control for this.  

`phiFactor1, phiFactor2` are the levels of `phi` that reflect moderate and major interventions, respectively.  

`cosAmp` is the amplitude of the seasonal variation in transmissibility. `pts_function` uses a shifted cosine function for seasonality that peaks on February 1 and troughs on August 1.  

`startDay, N, numComp, pICompt` are simulation parameters. `pICompt` is defined because the observed prevalence is restricted to observed cases, and not to unknown cases (`uI`) or cases before the are observed/not observed (`I`).  

`upDelay, downDelay` are the user-defined delays in days for imposing or lifting a policy intervention. They are used with `prevUp01`, etc. in `v`.  

`phiMoveUp, phiMoveDown` are the parameters in `(0,1]` that define the rate at which `phi` exponentially converges when increasing(decreasing) when imposing(lifting) a policy intervention. When `phiMoveUp` is close to 0, the convergence rate is slow, when it is close to 1, the convergence rate is fast.  

`pdDecay` is number of days over which `phi` decreases toward `phiNoAction` in the absence of policy interventions. If `pdDecay = -1`, there is no decay and physical distancing is always used (`phi = phiPhysicalDistancing`). The decay is not exponential, since that would result in the fastest change at the beginning. Instead the decay is based on a log, where the decay is slowest at the beginning and increases as `pdCounter` approaches `pdDecay`.  

`switchOffPolicies, switchOffDay` allow for switching off policies, and the date on which policies are switched off.

#### Walkthrough of pts_funScript
`pts_funScript` is the function that produces the long string which is `pts_fun`. It is written in C code in an R file, and uses `paste0` to transform the arguments into parameters. It is dense and not very readable, since it is neither formatted C code nor formatted R code, but a mishmash of both. This section of the readme walks through the code.

##### Initializing parameters
C requires initializing parameters, which takes up the top section of `pts_function`. The function uses `paste0` to transform the arguments (in R) to parameters (in C). I wrote this code in stages, so not all of the parameter names coincide with the arguments. I may eventually clean this up, but don't hold your breath.

```
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
    const int Nn = ",
    enn,
    "; /* number of nodes in trial */
    const int Nc = ",
    numComp,
    "; /* number of compartments */
    const int pIc =",
    pICompt,
    "; /* which compartment is pI */
````

##### Initializing other parameters/variables
`pts_funScript` uses a number of parameters and variables from the compartments dataframe `u` and the continuous variable dataframe `v`. Best practice in C code is to enumerate arrays using `enum`, as is done in the built-in SimInf C routines, but this occurs before the actual function, which I don't think is possible with a user-defined `pts_fun`. Therefore I used arguments (e.g. `pICompt`) or hard-coded numbers when selecting compartments in `u` and variables in `v`. As a reminder, C initializes array at 0, which simplifies loops, but is different from how R does it. Fortunately the rest of the code is one long string, and does not require additional `paste0` clauses.  

```
const double firstNode = ldata[0]; // first node in trial
    
const double pi = 3.141593;        // pi
    
const double phi = v[0];
const double prevUp01 = v[1];    /* variable capturing delay for moving up from lowest prev (0) to intermediate prev (1) */
const double prevUp12 = v[2];    /* variable capturing delay for moving up from intermediate prev (1) to highest prev (2) */
const double prevDown21 = v[3];  /* variable capturing delay for moving down from highest prev (2) to intermediate prev (1) */
const double prevDown10 = v[4];  /* variable capturing delay for moving down from intermediate prev (1) to lowest prev (0) */
const double kbPhase_0 = v[5];  /* logical: 0 means current stay-at-home, 1-3 means stay-at-home lifted */
const double policy_0 = v[6];    /* logical: 1 means policies are used, 0 means not */
const double prev_0 = v[8];      /* prevalence in previous time step */
const double prevState_0 = v[9]; /* intervention state before minor intervention */
const double pdCounter_0 = v[10]; // Counter for pdDecay
```
`firstNode` is the first node in the trial; it is constant across all groups in that trial. It is used to calculate the trial-wide prevalence.  

`v[0]` through `v[10]` create constants from the continous variables stored in `v`.  

I noticed when debugging that I ran into some problems if I didn't make sure that continous variables were explictly updated. This may be a problem with my coding, but regardless, the easiest fix was to make sure I updated all the continuous variables at least once. The SimInf implementation of `pts_fun` stores the updated continous variable in an array called `v_new`.
```
// update variables in case they arent updated later
v_new[1] = prevUp01;
v_new[2] = prevUp12;
v_new[3] = prevDown21;
v_new[4] = prevDown10;
v_new[5] = kbPhase_0;
v_new[6] = policy_0;
v_new[8] = prev_0;
v_new[9] = prevState_0;
v_new[10] = pdCounter_0;
```

#### Updating seasonal effect on transmissibility
```
/* Update seasonal effect: one cosine wave over course of year, peak on Feb 1 (day 32), trough on Aug 1 */
    v_new[7] = cosAmp * cos(pi * (t - 32 + startDay) / 182.5) + 1;
```

#### Computing observed prevalence for each trial
The main purpose for the post time step function is to compute the trial-wide observed prevalence and update any policy interventions. This snippet of code uses the pointer functionality of C, which I don't understand super well, but I was able to successfully model this off of `SISe_sp.c` from the SimInf package. Basically using an asterisk \* and an ampersand \& creates an array that allows for greater flexibility in specifying which element of the array is selected. The SimInf implementation stores `u` and `v` as single-dimensional arrays (as is common in C), so to point to the correct element, the code needs to use the `row * (number of columns) + index` approach that R tends to avoid by using dataframes and lists. However, it is very fast.  

```
/* Initalize the compartment state vector to find number of individuals in neighbors */
    const int *u_0 = &u[-Nc*node];
````
`Nc` is the number of compartments in the model, i.e. the number of columns.  

`node` is built-in to the `pts_fun` implementation as SimInf cycles over the array. `node` is the row.  

Therefore `&u[-Nc*node]` is the index of the first node in the whole array that is `u`, regardless of what the current node is. `u_0` is the array holding all the compartments for each node. This allows the code to get the value of `u_0` at any element by specifying the correct index.

Continuing the code snippet, starting by initializing variables:

```
int j = firstNode, k;
double N_j = 0, I_j = 0, prev, kbPhase = kbPhase_0, policy=policy_0, pdCounter = pdCounter_0, pdB;
    
for (j = firstNode; j < firstNode + Nn; j++) {
      
  /* Add infected in node j to total */
  I_j += u_0[j * Nc + pIc];  // pI is compartment 6
  
  /* Add pop in node j, compartment k to total */
  for (k = j * Nc; k < (j + 1) * Nc; k++)
    N_j += u_0[k];
  
  N_j -= u_0[((j + 1) * Nc) - 1]; /* subtracts M from total pop */
  }
  /* Compute prevalence across all nodes in trial */
  if(prevType == 0) {
   prev = I_j;
  } else {
   prev = I_j/N_j;
  }
  v_new[8] = prev; /* recording prevalence */
````

`firstNode` is the first node in the trial; `Nn` is the number of nodes in the trial. By starting `j=firstNode` and incrementing until `j < firstNode + Nn`, the code iterates over all the nodes in the trial that contains the current node.  

`I_j += u_0[j * Nc + pIc];` increments the number of observed infections by the number of `pI` in each node in the given trial.  `N_j += u_0[k]`, iterated over all compartments `Nc` in the node, computes the population of the trial. In practice, the population is constant, but it is possible include birth rates, etc. `N_j -= u_0[((j + 1) * Nc) - 1];` subtracts the number of people who have died to get the living population.  

Then the trial observed prevalence is computed, depending on if prevalence should be count or a proportion.

#### Updating fixed-time variables
`kbPhase` and `policy` are changed at fixed times. This snippet of code updates them. `t` is built in to the SimInf implementation.
```
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
```

