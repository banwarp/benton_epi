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

##### Mass entry parameters
The code allows for a single mass entry event. The mass entry event can be turned off by setting `massMassEntryNodes=0` in the user-defined parameters.  

`maxMassEntryNodes <- min(length(which(nodeGroupList %in% massEntryNodeGroups)),maxMassEntryNodes)` is user-defined, but also limited by the number of available nodes in the set of node groups chosen to host the mass entry event.  

`mRProp` through `mMProp` are created to distributed the mass entry individuals among the different departments.  `massEntryPropTable` gathers these proportions together into a single table.

##### Compartments
1. `S` Susceptible
2. `E` Exposed (pre-infectious)
3. `I` Infectious
4. `R` Recovered
5. `Im` Immune
6. `H` Hospitalized
7. `Is` Isolated
8. `cumI` Cumulative infected. Used for plotting, not for disease dynamics
9. `M` Deceased

##### Global and local parameters
SimInf uses `gdata` and `ldata` for the model transitions. Most of these parameters are user defined, except for `beta`, `isoRate`, and `betaIsolated`.  

`beta` is calculated as a product of the user-defined `R0` and `infectiousPeriod`.  

`isoRate` is user defined. However, assuming that all hospitalized individuals are isolated, `isoRate` is overwritten by the maximum of `isoRate` and `hospRate`.  

`betaIsolated` is the transmissibility of the virus among isolated individuals, and it is the product of the user-defined `RIsolated` and `isoPeriod`, similar to `beta`.  

`firstNode` is the first node in a trial. Used for getting the correct indices in various subroutines.
`betaRandomizer` is generated from a uniform distribution with mean 1 and a user-defined spread. Used to generate a distribution of transmissibility rates (`beta`) to add variability to the trials in the simulation.

###### Model transitions
Transitions beginning or ending with `@` represent entry/exit from/to the empty set.
1. `@ -> mu*(S+E+I+R+Im) -> S' Natural population growth 
2. `R -> reSuscepRate*tempImmPeriod*R -> S` Some proportion of recovered become susceptible again after a peiod of temporary immunity.
3. `S -> nu*S -> @`, `E -> nu*E -> @`, `I -> nu*I -> @`, `Is -> nu*Is -> @`, `R -> nu*R -> @`, `Im -> nu*Im -> @` Natural death rate, not associated with COVID-19.
4. `S -> ((1/phi)*beta*betaRandomizer*season*I+betaIsolated*Is)*S/(S+E+I+R+Im)-> E` Transition from Susceptible to Exposed.
 - `beta` is the baseline transmissibility rate.
 - `(1/phi)` scales beta according to the intensity of policy intervention and/or physical distancing. A more intense intervention or more physical distancing increases `phi` which decreases the effective transmissibility rate.
 - `season` is the seasonality factor, which peaks in February and troughs in August.
 - `betaRandomizer` creates a small distribution around `beta` for each separate trial to increase variability between trials.
 - `betaIsolated` is the transmissibility rate for isolated individuals. It is close to zero but positive to represent incomplete isolation/quarantine or slightly delayed isolation/quarantine.
5. `E -> (1-isoRate)*exposedPeriod*E -> I + cumI` Transition from Exposed to Infectious. Includes `isoRate` to represent Exposed individuals being isolated before they become infectious. `cumI` tracts the total number of infected individuals.
6. `E -> (isoRate-hospRate)*exposedPeriod*E -> Is + cumI` Transition from Exposed to Isolated. Isolated individuals have a much lower transmissiblity factor, which slows the epidemic spread. Includes `hospRate` to represent that some Isolated individuals are hospitalized during their isolation.
7. `E -> hospRate*exposedPeriod*E -> H` Transition from Exposed to Hospitalized. Note that we assume Hospitalized individuals have zero transmissibility. This can be changed in the code if (for example) you want to represent some infection before hospitalization.
8. `I -> (1-nonHospDeathRate)*infectiousPeriod*I -> R` Recovery from being Infectious, less COVID-19 deaths among the Infectious (defined to excluded hospitalized COVID-19 patients).
9. `I -> nonHospDeathRate*infectiousPeriod*I -> M` Deaths among non-hospitalized COVID patients
10. `H -> (1-hospDeathRate)*hospPeriod*H -> R` Recovery among Hospitalized
11. `H -> hospDeathRate*H -> M` Deaths among hospitalized COVID patients. Note that we assume the death rate among Isolated = 0, since they would have been hospitalized. This can be changed in the code.
12. `Is -> isoPeriod*Is -> R` Recovery from COVID among isolated, assumes non-Hospitalized Isolated all recover.
13. `R -> (1-reSuscepRate)*tempImmPeriod*R -> Im` Most Recovered become permanently immune after a period of temporary immunity. Some fraction (`reSuscepRate`) become Susceptible again.  
