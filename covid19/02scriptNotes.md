## Description of script segments in the Covid-19 scripts
#### Introduction
This set of R scripts was developed for local modeling of the COVID-19 disease in Benton County, Oregon, and other small geographies. The scripts use a modified SEIR model with partially-mixed, but otherwise homogeneous population nodes. This markdown document describes the process by which the script pre- and post- processes data.

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
- With the result, you can plot trajectories for the different compartments and/or the evolutino of the continuous variables.

### Annotated code
(Not every line of code is included; just the ones that need explanation)
#### Libraries and subroutines; saving parameters
The script uses a number of packages in addition to SimInf:
`ggplot2`, `reshape2`, `SimInf`, `data.table`, `dplyr`, `Matrix`, `zoo`  

The script also relies on two subroutines: `eventFunctionsDATE.R`, which generates events for the simInf model, and `simInfPlottingFunctionDATE.R`, which aggregates the different nodes in each trial and plots the median and confidence spread for different compartments/continuous variables.  

Then the script saves the parameters for reference.

#### Setting additional parameters
Many of the user-defined parameters are used to set additional, internal parameters.

##### Date parameters
The date parameters are changed to numeric values, e.g. `startofSimDay <- as.numeric(as.Date(simDate)) - as.numeric(as.Date("2020-01-01"))`.

##### Trial and node parameters
In order to make sure that the script runs within a reasonable amount of time, an upper limit is placed on the number of trials and the number of nodes. The limits can be changed in the code.
```
numTrials <- min(100,numTrials) # capping number of trials at 100
N <- min(1000,N) # capping number of nodes per trial at 1000
```
Additional parameters are generated from `trialPop`, `N`, and `nodeGroupList`.
```
nodeTrialMat <- data.table(node=c(1:(NnumTrials)),
                           nodeGroup = rep(nodeGroupList,numTrials),
                           trial=rep(1:numTrials,each=N))
```
`nodeTrialMat` is used in subroutines when specific nodes needs to be linked to trials and/or node groups.  

`NList <- split(nodeTrialMat$node,nodeTrialMat$trial)` is created in order to use `lapply` on different trials.
```
nodeAndGroupList <- split(nodeTrialMat,nodeTrialMat$nodeGroup)
nodeAndGroupList <- lapply(nodeAndGroupList,function(x) split(x$node,x$trial))
```
`nodeAndGroupList` is created to use `lapply` on different node groups.  

`maxNodePop <- floor(trialPop/N)` is created for distributing populations across nodes within trials.

##### Compartments
The model is based on SEIR dynamics, with a couple of twists.  

First, as we become more aware of COVID-19, we are doing a better job monitoring for symptoms and testing. Therefore, after an initial infectious stage `I`, there is a probability that the infection is identified and the individual is moved to the post infectious `pI` compartment, whether they are much less infectious. Or, if the infection is not identified, they move to the `uI` compartment and remain infectious in the community. I chose a single chance to identify the infection rather than a recurring probability to simplify the transitions. This also leaves open the possibility of choosing three different infectious levels, one for initial infectious, one for post infectious, and one for unknown infectious.  

Second, some proportion of infectious individuals require hospitalization. These individuals have differential recovery rates compared to the non-hospitalized infectious compartment.  

Finally, some proportion of the recovered can again become susceptible after a period of time.

1. `S` Susceptible
2. `E` Exposed (pre-infectious)
3. `I` Initial infectious
4. `uI` Unknown infectious
5. `R` Recovered
6. `Im` Immune
7. `pI` Post infectious
8. `H` Hospitalized
9. `cumI` Cumulative infected. Used for plotting, not for disease dynamics
10. `M` Deceased

##### Global and local parameters
SimInf uses `gdata` and `ldata` for the model transitions. Most of these parameters are user defined, except for `betaI`, `betaU`, and `betaP`.  

`betaI` is calculated as a product of the user-defined `R0I` and `initInfectiousPeriod`. `betaU` and `betaP` are defined similarly.

`firstNode` is the first node in a trial. Used for getting the correct indices in various subroutines.
`betaRandomizer` is generated from a uniform distribution with mean 1 and a user-defined spread. Used to generate a distribution of transmissibility rates (`beta`) to add variability to the trials in the simulation.

###### Model transitions
Transitions beginning or ending with `@` represent entry/exit from/to the empty set.
1. `@ -> mu*(S+E+I+R+Im) -> S` Natural population growth 
2. `R -> reSuscepRate*tempImmPeriod*R -> S` Some proportion of recovered become susceptible again after a peiod of temporary immunity.
3. `S -> (betaI*(1/phi)*I+betaP*pI+betaU*(1/phi)*uI)*season*betaRandomizer*S/(S+E+I+R+Im)-> E` Transition from Susceptible to Exposed.
    - `betaI`, `betaP`, and `betaU` are the baseline transmissibility rates for `I`, `pI`, and `uI`, respectively.
    - `(1/phi)` scales beta according to the intensity of policy intervention and/or physical distancing. A more intense intervention or more physical distancing increases `phi` which decreases the effective transmissibility rate.
    - `season` is the seasonality factor, which peaks in February and troughs in August.
    - `betaRandomizer` creates a small distribution around `beta_` for each separate trial to increase variability between trials.
4. `E -> (1-isoRate)*exposedPeriod*E -> I + cumI` Transition from Exposed to Infectious. Includes `isoRate` to represent Exposed individuals being isolated before they become infectious. `cumI` tracts the total number of infected individuals.
5. `E -> (1-hospRateExp)*isoRate*exposedPeriod*E -> pI + cumI` Transition from Exposed to post-infectious. Includes `hospRate` to represent that some post-infectious individuals are hospitalized during their isolation.
6. `E -> hospRateExp*isoRate*exposedPeriod*E -> H + cumI` Transition from Exposed to Hospitalized. Note that we assume Hospitalized individuals have zero transmissibility. This can be changed in the code.
7  `I -> monitoringSuccess*(1-hospRatePost)*initInfectiousPeriod*I -> pI` Transition from infectious to post-infectious (non-hospitalized). This represents individuals identified through symptoms or testing after being infectious for a short period.
8. `I -> monitoringSuccess*hospRatePost*initInfectiousPeriod*I -> H` Transition from infectious to post-infectious (hospitalized).
9. `I -> (1-monitoringSuccess)*initInfectiousPeriod*I -> uI` Transition from infectious to unknown infectious. This represents individuals who are never identified as infectious and continue to circulate in the community. It is assumed that none of these individuals are hospitalized or die.
10. `pI -> (1-nonHospDeathRate)*postInfectiousPeriod*pI -> R` Recovery among post-infectious
11. `pI -> nonHospDeathRate*postInfectiousPeriod*pI -> M` Mortality among post-infectious (non-hospitalized)
12. `uI -> (1-nonHospDeathRate)*unknownInfectiousPeriod*uI -> R` Recovery among unknown-infectious
13. `uI -> nonHospDeathRate*unknownInfectiousPeriod*I -> M` Deaths among non-hospitalized COVID patients (unknown-infectious)
14. `H -> (1-hospDeathRate)*hospPeriod*H -> R` Recovery among hospitalized
15. `H -> hospDeathRate*H -> M` Deaths among hospitalized COVID patients
16. `R -> (1-reSuscepRate)*tempImmPeriod*R -> Im` Most Recovered become permanently Immune after a period of temporary immunity. Some fraction (`reSuscepRate`) become Susceptible again.
17. `S -> nu*S -> @`, `E -> nu*E -> @`, `I -> nu*I -> @`, `pI -> nu*pI -> @`, `uI -> nu*uI -> @`,`R -> nu*R -> @`, `Im -> nu*Im -> @` Natural death rate, not associated with COVID-19.  

##### Initial state of compartments
The initial susceptible population is randomly distributed across the nodes in each trial. Most nodes will have a population approximately equal to `.8*trialPop/N`. The distribution is skewed left, so that some nodes will have a much smaller population, and a decent number of nodes will have a population closer to `trialPop/N`. Is this necessary? Probably not, and it could be replaced with constant `trialPop/N` for each node.
```
S0 <- matrix(floor(maxNodePop/log(maxNodePop)*log(runif(NnumTrials,min=1,max=maxNodePop))),ncol=numTrials)
S0 <- matrix(apply(S0,2,function(x) x+round(((trialPop-sum(x))/length(x)),0)))
```
The code allows the user to select which node groups contain initial infectious individuals and the maximum number of nodes within those node groups that can have 1 or more initial infectious individuals. This allows the user to model scenarios where the initial infection is highly clustered or widely spread.
```
eligibleInodes <- nodeTrialMat[which(nodeTrialMat$nodeGroup %in% I0nodeGroups & nodeTrialMat$trial==1),"node"]$node
Inodes <- lapply(0:(numTrials-1),function(x) (x*N)+sample(eligibleInodes,ceiling(maxINodeProp*length(eligibleInodes))))
I0indices <- lapply(Inodes,function(x) if(length(x)>1){sample(x,rbinom(1,2*I0Pop,.5),replace=TRUE)}
                                              else{rep(x,I0Pop)})
I0indices <- unlist(lapply(I0indices,function(x) table(x)))
I0 <- matrix(0,nrow=NnumTrials)
I0[as.numeric(names(I0indices))] <- I0indices
```
The default is to initialize all compartments except `S`, `I`, and `cumI` to 0.
```
u0 <- data.frame(
  S = S0,
  E = rep(E0Pop,NnumTrials),
  I = I0, # random number of initial infectious, allows for 0, capped at maxIPop for each trial
  R = rep(R0Pop,NnumTrials),
  Im = rep(Im0Pop,NnumTrials),
  H = rep(H0Pop,NnumTrials),
  Is = rep(Is0Pop,NnumTrials),
  cumI = I0,
  M=rep(M0Pop,NnumTrials)
)
```

##### Initial state of continous variables
`phi` is initialized to `phi0`, which is a user defined parameter.  

Whenever prevalence crosses a certain threshold, a timer starts before an intervention is imposed or lifted. At the beginning of the simulation, the timers are initialized to the ready state.
```
prevUp01 = rep(upDelay+1,NnumTrials),      # initialized variable for capturing delay in response when prevalence increases
prevUp12 = rep(upDelay+1,NnumTrials),      # variable for capturing delay in response when prevalence increases
prevDown21 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
prevDown10 = rep(downDelay+1,NnumTrials),  # variable for capturing delay in response when prevalence decreases
```
In the state of Oregon, the stay-home stay-healthy orders were lifted in three phases. `kbPhase` initializes the phase.  

In order to explore counterfactuals, the user can choose to switch off policies on a certain date. `policy = 1` intializes the state to be able to use policies.  

`season = 1` centers the seasonal scaling of `beta` about `1`.  

The model records the daily `prevalence`, which is initialized to 0 for simplicity. It could be initialized to the prevalence in each trial, but that would not make a difference when the model runs.  

Various parameters, states, and variables depend on whether prevalence is increasing or decreasing. `previousState` records what category the `prevalence` was in during the previous time step: baseline, minor, or major intervention. For example, if `previousState` was major intervention, and the current prevalence is now below the major threshold, the post-time-step function will start the timer for lifting the major intervention.  

`RTee` tracks the effective reproduction number so that the user can plot it if they want to. It does not affect the disease dynamics or the model directly.

`pdCounter` is a timer variable for the decay of physical distancing. The model assumes that if prevalence is low (below the minor threshold), people will eventually stop physically distancing. `pdCounter` is used to track the decay at each time step.

##### Post-time-step function
A whole .md could be (and will be) written about the post-time-step function (pts_function). In brief, the pts_function tells SimInf how to change continous variables (with the option of other effects) after every time step. The pts_function is written in C, which is nothing if not verbose. The function is stored as a single character and passed to the SimInf routine. To save space, I wrote the pts_function in a separate script and call pts_funScript to generate it. The parameters are passed through from user-defined parameters or paramters already described in this .md.
```
pts_fun <- pts_funScript(
    phiPhysicalDistancing = (gdata$beta/gdata$infectiousPeriod)/RPhysicalDistancing, # Phi reflecting that physical distancing and contact tracing will reduce R0 even without stay-at-home orders
    phiNoAction = (gdata$beta/gdata$infectiousPeriod)/RNoAction,                     # Phi reflecting no actions
    maxPrev1 = maxPrev1,                                                             # Maximum prevalence before instituting minor intervention
    maxPrev2 = maxPrev2,                                                             # Maximum prevalence before instituting major intervention
    phiFactor1 = (gdata$beta/gdata$infectiousPeriod)/RTarget1,                       # Target for the reduction in R0 under minor intervention
    phiFactor2 = (gdata$beta/gdata$infectiousPeriod)/RTarget2,                       # Target for the reduction in R0 under major intervention
    cosAmp = cosAmp,                                                                 # Amplitude of seasonal variation in beta
    startDay = startofSimDay,                                                        # start of simulation
    kbDay1 = kbDay1,                                                                 # date of first phase
    kbDay2 = kbDay2,                                                                 # date of second phase
    kbDay3 = kbDay3,                                                                 # date of third phase
    enn = N,                                                                         # number of nodes in a trial
    numComp = length(compartments),                                                  # number of compartments
    prevType = ifelse(maxPrev1 < 1,1,0),                                             # prevalence type. 0 = count, 1 = proportion
    upDelay = upDelay,                                                               # Number of days after prevalence passes threshold until minor/major intervention
    downDelay = downDelay,                                                           # Number of days after prevalence drops below threshold until intervention lifted
    phiMoveUp = phiMoveUp,                                                           # Rate at which phi increases when interventions are imposed
    phiMoveDown = phiMoveDown,                                                       # Rate at which phi decreases when interventions are lifted
    pdDecay = pdDecay,                                                               # Rate at which phi decreases toward 1 in the absence of interventions. Represents gradual relaxation of physical distancing
    switchOffPolicies = switchOffPolicies,                                           # Logical variable to switch off policies (used for counterfactuals)
    switchOffDay = as.numeric(as.Date(switchOffDay)) - as.numeric(as.Date("2020-01-01")) - startofSimDay # day policies would be switched off (used for counterfactuals)
  )
  ```

##### Events
Events specify time, compartment, node, number, and type (entry, exit, transfer, shift). For full details, refer to the [SimInf vignette](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf).
In order to use events, the user must define a matrix(ces) for SimInf that tells the routine which compartments the events happen to. `E` is for events that primarily affect transfers between nodes; `N` is for events that primarily affect transfer between compartments.
```
# Events matrix
# leaves off Is, H, cumI, and M from last column
E <- cbind(diag(length(compartments)),c(rep(1,length(compartments)-4),rep(0,4)))
dimnames(E) <- list(compartments,c(1:ncol(E)))

# Shift matrix
Nmat <- matrix(c(1,rep(0,length(compartments)-1)),ncol=1)
dimnames(Nmat) <- list(compartments,1)
```
The columns of `E` (`N`) correspond to different types of events, and the populated rows specify which compartments are impacted in each event.

##### Parachute events
The model allows for infectious individuals to parachute into the populations. The user can define how many parachute events there are, how many infectious individuals parachute in each time, and which node groups can receive a parachuter. In order to try to approximate reality, I sampled parachute events from a truncated chi square distribution with degree 4. If a uniform distribution is preferred, this can be changed in the code.  

Each trial gets its own set of randomly distributed parachute events. Once the number of events for a trial is set, the specific node(s) is selected from among the eligible node groups (listed in `parachuteNodeGroups`). Then the `parachuteEvents` data frame is created.

```
# number of events in each trial
  numParachuteList <- lapply(1:numTrials,function(x) rpois(1,parachuteRate*max(tspan)))
  parachuteNTM <- nodeTrialMat[nodeTrialMat$nodeGroup %in% parachuteNodeGroups,]
  pNList <- split(parachuteNTM$node,parachuteNTM$trial)
  
  if(max(unlist(numParachuteList)) > 0){
    # parachute events data frame
    parachuteEvents <- data.frame(
      event = "enter",
      time = unlist(lapply(numParachuteList, function(x) sample(parachuteDist,x,replace=TRUE))),
      node = unlist(lapply(c(1:numTrials),function(x) sample(pNList[[x]],numParachuteList[[x]],replace=TRUE))),
      dest = 0,
      n = parachuteNum,
      proportion = 0,
      select = which(compartments == "I"),
      shift = 0
    )
  }
```

##### Transfer Events
In order to model partial mixing, I split each trial into a set of nodes, classified by node groups. The default is for there to be 1 node group, but the user can define a list of node groups. I use transfers between nodes to represent the broader mixing of the population. I assume that transfers within a node group are more common than transfer between node groups. This represents that people tend to mix with the same set of people in their node, occasionally mixing with others in their node group, and rarely mixing with still others in other node groups. For example, people in one suburb, occasionally mixing with people in their neighboring city, rarely mixing with people in other cities.  

As with parachute events, the number and date of in-node and out-node transfers is randomly distributed and different in each trial. The distribution is uniform with the number of events following a Poisson distribution. The script uses the `transferFunction` subroutine to populate the `transferEvents` data frame.
```
numInTransferList <- lapply(1:numGroups,function(x) 
  lapply(1:numTrials, function(y) rpois(1,inGroupTransferRate[x]*inGroupTransferNodeNum[x]*max(tspan))))
  
inTransferEventsList <- lapply(1:numGroups, transferFunction, ...)
# See script for complete code

numOutTransferList <- lapply(1:numGroups,function(x) 
    lapply(1:numTrials, function(y) rpois(1,outGroupTransferRate[x]*outGroupTransferNodeNum[x]*max(tspan))))
    
  outTransferEventsList <- lapply(1:numGroups, transferFunction...)
# See script for complete code
```
In order to keep a lid on transfers, I cap the number of people who can transfer based on the smaller population of the origination and destination node. I assume that a small node would not receive a huge in-transfer. Maybe this doesn't make much difference - I haven't explored it yet.
```
if(nrow(transferEvents)>0){
  Sgrid <- data.frame(node=1:length(u0$S),nodePop = u0$S)
  transferEvents <- merge(transferEvents,Sgrid,by="node",all.x=TRUE)
  names(Sgrid) <- c("dest","destPop")
  transferEvents <- merge(transferEvents,Sgrid,by="dest",all.x=TRUE)
  transferEvents$maxProp <- transferEvents$destPop/transferEvents$nodePop
  transferEvents[transferEvents$maxProp<1,"proportion"] <- transferEvents[transferEvents$maxProp<1,"maxProp"]*transferEvents[transferEvents$maxProp<1,"proportion"]
  transferEvents <- transferEvents[,c("event","time","node","dest","n","proportion","select","shift")]
}
```

##### Super-spreader events
There is a small but nonzero probability that any given infected person becomes a super-spreader, i.e. does not conform to the parameterized `R0` in terms of number of generated infections. To model this, I allow for zero, one, or more super-spreader events. The super-spreader parameters are in lists to allow for more than one. For each event, the user selects the date. This could be randomized using `superDate = sample(maxT,n)`, where `n` is the number of events. Each event also has a spread of days when the infections take place. The event transitions a given number of individuals from the Susceptible compartment to the Exposed compartment. As with parachuters, the newly Exposed can be assigned to certain node groups and a quantity of nodes within those groups, representing that the super-spreader does not come in contact with everyone in the population.  

Note: There is an unresolved problem with Super-spreader events that sometimes occurs when running the code. The SimInf package wisely expects that the population transfered is no larger than the population of the node/compartment. In order to randomize the super-spreader event, this script cannot project in advance the population of the node/compartment, and therefore runs a non-negligible risk of trying to transfer more individuals from I to E then exist in I. When this happens the function returns an error. You can always call the function again and hope it doesn't happen, or you can change the code to use proportions instead of numbers, but then you can't specify exactly how many individuals transition from I to E, only the proportion. This is not a problem with SimInf, but ideas on the best way to resolve it would be appreciated.
```
eligibleSuperNodes <- lapply(superNodeGroups,function(x) nodeTrialMat[which(nodeTrialMat$nodeGroup %in% x & nodeTrialMat$trial==1),"node"]$node)

superNodeList <- lapply(1:length(superNodes),
                        function(x) if(length(eligibleSuperNodes[[x]])>1)
                          {sort(sample(eligibleSuperNodes[[x]],min(length(eligibleSuperNodes[[x]]),superNodes[x])))}
                          else{eligibleSuperNodes[[x]]})

superEventList <- lapply(1:length(superInfections), superFunction,...)
# See script for complete code
```

##### Mass entry events
I developed this model for local projections of a county whose population changes by 25% multiple times a year. Therefore I wanted to be able to model a mass entry event to understand how the disease dynamics may change under a shock to the system like a 25% increase in the susceptible population and the introduction of a number of new infections. As in the super-spreader event the user can define which node groups and how many nodes receive the mass entry. I have not written in a mass exit event or multiple mass entry events, but these can be added to the code by using the super-spreader/mass entry events as templates.  

###### Mass entry parameters
The code allows for a single mass entry event. The mass entry event can be turned off by setting `massMassEntryNodes=0` in the user-defined parameters.  

`maxMassEntryNodes <- min(length(which(nodeGroupList %in% massEntryNodeGroups)),maxMassEntryNodes)` is user-defined, but also limited by the number of available nodes in the set of node groups chosen to host the mass entry event.  

`mRProp` through `mMProp` are created to distributed the mass entry individuals among the different departments.  `massEntryPropTable` gathers these proportions together into a single table.

```
eligibleMassEntryNodes <- nodeTrialMat[which(nodeTrialMat$nodeGroup %in% massEntryNodeGroups & nodeTrialMat$trial==1),"node"]$node
    massEntryNodeList <- sort(sample(eligibleMassEntryNodes,min(length(eligibleMassEntryNodes),maxMassEntryNodes)))
    
# exclude isolated, hospitalized, cumulative, and deceased by length(compartments)-4
massEntryEventList <- lapply(c(1:(length(compartments)-4)),massEntryEventFunction,...)
# See script for complete code
# Add introduced infections to cumulative total
massEntryCumI <- massEntryEvents[which(massEntryEvents$select == which(compartments=="I")),]
massEntryCumI$select <- which(compartments=="cumI")
massEntryEvents <- rbind(massEntryEvents,massEntryCumI)
```

##### Building and running the SimInf model
The SimInf model requires just one events data frame:
```
allEvents <- data.frame(NULL)
if(max(unlist(numParachuteList)) > 0) { allEvents <- rbind(allEvents,parachuteEvents) }
if(nrow(transferEvents)>0) { allEvents <- rbind(allEvents,transferEvents) }
if(maxMassEntryNodes > 0) { allEvents <- rbind(allEvents,massEntryEvents) }
if(length(superInfections)>0) { allEvents <- rbind(allEvents,superEvents)}
```

The `mparse` command is built into SimInf and allows the user to define their own model based on the defined transitions. All of the arguments are passed to the SimInf routine, then the model is run.
```
model <- mparse(transitions = transitions, compartments = compartments,
events = allEvents, E=E,N=Nmat,
gdata = gdata, ldata=ldata,
u0 = u0, v0=v0, tspan = tspan, pts_fun = pts_fun)

result <- run(model)
```

##### Output and plots
The output of `run()` is a model object that contains trajectories of the compartments, evolution of the continuous variables, and meta data about the arguments. To extract the trajectories use `trajectory(result)`. You can also directly plot the model using SimInf's implementation of `plot`. See the [vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf) for more details. However, the built in plotter treats each node as its own trial, but I need to sum compartments across all nodes in the trial. Therefore I built a plotting subroutine that sums all the nodes in a given trial (or a given nodeGroup), creates a credibile spread (aka confidence interval) from the different trials, and plots the median and spread of the simulation.
```
trajPlotInfections <- simInfPlottingFunction(
  result = result,                           # model result
  table = "U",                               # which table: U or V
  compts= c("I_Is_H"),                       # compartments that will be plotted
  groups = NULL,                             # List of groups to aggregate and plot
  sumGroups = TRUE,                          # Automatically sums over all groups
  uNames = names(u0),                        # list of compartments
  vNames = NULL,                             # list of variables
  rollM = rollM,                             # number of days for rolling mean
  confIntv = confIntv,                       # confidence interval for plotting spread
  nTM = nodeTrialMat,                        # node-Trial matrix
  tS = tspan,                                # length of simulation
  enn = N,                                   # number of nodes per trial
  nT = numTrials,                            # number of trials in simulation
  startDate = startofSimDay,                 # start date of simulation
  dateBreaks = dateBreaks,                   # plot parameter: Date axis format
  titleString = paste0("Active infections in ",titleString," County"),        # plot parameter: Title of plot
  xString = "Date",                          # plot parameter: Title of x axis
  yString = "Number of infections",          # plot parameter: Title of y axis
  lString = "Compartment"                    # plot parameter: Title of legend
)
```
Details of `simInfPlottingFunction` will be explained in a different .md.  

The overall function returns the model, so that the user can extract the different elements of the model, plot additional trajectories, etc.

```return(result)```

### Further reading
Read [03COVID-19 Vignettes.md](03covid19Vignettes.md) for examples of how to use the code.
