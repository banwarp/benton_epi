## COVID-19 Vignettes - Post Time Step Function Readme
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

#### Post Time Step Function
The post time step function runs after each time step and updates the continuous variables stored in the data frame `v`. These continuous variables, which are initialized in `v0`, are used in the compartment transitions along with `gdata` and `ldata`, but unlike `gdata` and `ldata`, they can vary across time, either exogeneously to the disease conditions (such as a seasonality factor) or in response to disease conditions (such as a policy intervention). The post time step function is written and implemented in C.  

##### Initializing v0
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

#### Explanation of each continuous variable
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

`prevUp01, prevUp12, prevDown21,prevDown10` are delay trac
