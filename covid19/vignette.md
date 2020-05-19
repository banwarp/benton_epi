## A basic vignette to describe how to use the covid19 set of R scripts
### Introduction
This set of R scripts was developed for local modeling of the COVID-19 disease in Benton County, Oregon, and other small geographies. The scripts use a modified SEIR model with partially-mixed, but otherwise homogeneous population nodes. This vignette describes how to use the scripts.

### Acknowledgements:
I'd like to thank the developers of the SimInf package for writing the code that these scripts use and for assisting me in adapting their code to meet the requirements of COVID-19 modeling.  

Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12  

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales. International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

[SimInf vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf)  
[SimInf technical documentation](https://cran.r-project.org/web/packages/SimInf/SimInf.pdf)  
[SimInf git repository](https://github.com/stewid/SimInf)

### Conceptual schematic of SimInf
![Conceptual schematic of SimInf](images/simInfSchematic.png)
The conceptual schematic of the SimInf model that I use has three stages: Building the model, Running the model, and the Result.
#### Building the model:
- Specify the model with parameters for disease spread, population compartments, and transitions between compartments.
- Set the initial states in the compartments and the initial values of the continuous variables.
- (If the model used is not one of the built-in model) Write the post-time-step function that changes continuous variables after every time step.
- Generate a dataframe of events that can introduce or remove individuals to different compartments, transfer them between nodes, or shift them between compartments, at different times.
#### Run the model:
- Advance individuals through compartments using the transitions and disease parameters. Record the compartment counts and continous variable values for the disease trajectory.
- Use the events dataframe to move individuals in/out/between nodes and/or compartments.
- Advance the time step.
- Change continuous variables according to the post-time-step functions.
- Repeat for every time step.
#### Result:
- Trajectories of the compartments across the timespan.
- Evolution of continous variables across the timespan.
- With the result, you can plot trajectories for the different compartments and/or the evolutino of the continuous variables.

### covid19_SimInf_5.16.2020.R (or updated date)
The covid19_SimInf_DATE.R script pre- and post-processes data and sets parameters to use in the SimInf model. The script has a number of pre-processing steps and associated parameters, which are discussed below:
#### Simulation parameters
The script relies on subroutines that need to be accessed. Set the folder where they are stored using the `folderpath` argument. This folderpath is also where the exported graphs will be saved. The `simID` argument is the heading for the different graphs.  

The arguments `simDate` and `maxT` set the date parameters - when the simulation starts and how long it runs (in days).  

The argument `numTrials` sets the number of trials in the simulation. The script limits the number of trials to a maximum of 100. This limit can be changed in the code itself.
#### Population parameters
The model assumes a homogenous, partially-mixed, mostly closed population. Homogenous because every individual has the same disease dynamics (chance of infection/dying, length of infectious period, etc.). Partially-mixed because the total population is divided into nodes. The disease evolves independently in different nodes, with occasional transfers between nodes to expose other nodes to disease. Mostly closed because no individuals leave the population except through death, and no individuals enter the population except for the occasional parachuting of infectious individuals or planned entry events.  

The `trialPop` argument is the total population of initially susceptible individuals at time `t`.  

`N` is the number of nodes. The node population is randomly generated with a mean population of approximately `0.8*trialPop/N`. The distribution is slightly skewed left. The distribution of the trial population can be changed in the code itself.  

`nodeGroupList` is an optional parameter to split the nodes into different groups. This represents subpopulations that are more likely to mix within the group compared to across groups. For example, distinct cities can be grouped into different nodeGroups. `nodeGroupList` is either `NULL` or a vector of length `N` with entries equal to the node group for each node.  

`I0Pop` through `M0Pop` are the initial populations of the non-susceptible compartments. All of these compartments except `I0Pop` are set to 0 by default.  

`maxINodeProp` sets an upper limit on how many nodes can have a positive initial infectious population, to represent that the disease is not evenly distributed. `I0nodeGroups` allows the distribution of initial infectious to be restricted to certain node groups.
#### Disease dynamics parameters
The disease dynamic parameters are split into three groups: global data (`gdata`), local data (`ldata`), and continuous variables (`v0`). Global data applies to all nodes and trials in the simulation. Local data can be specified to the individual trial, nodeGroup, or node. Neither global nor local data change during the simulation. Continuous variables can be specified locally, and can change during the simulation. You can change where parameters are specific in the code itself. For example, the default is for the basic reproduction number to be a global parameter, but you could make it a local parameter or a continuous variable. 
##### Global data parameters
```
R0 = 2.9                                # Basic reproduction number
infectiousPeriod = 1/8                  # Reciprocal of infectious period
exposedPeriod = 1/4                     # Reciprocal of exposed period
isoRate = .125                          # Proportion of exposed who are identified and isolated before they become infectious
isoPeriod = 1/10                        # Reciprocal of length of isolation period
RIsolated = .125                        # Reproduction number of isolated infectious
hospRate = .033                         # Proportion of infectious/isolated that are hospitalized
hospPeriod = 1/14                       # Reciprocal of length of hospitalization period
nonHospDeathRate = 0                    # Non-hospitalized fatality rate
hospDeathRate = .125                    # Hospitalized fatality rate NOT CASE FATALITY RATE
reSuscepRate = .1                       # Proportion of recovereds who eventually become susceptible again
tempImmPeriod = 1/100                   # Reciprocal of temporary immunity period, after which R becomes Im or S
mu = 0                                  # Natural birth/susceptible immigration rate
nu = 0                                  # Natural non-Covid death rate
```
##### Local data parameters
` R0Spread = .1` is a parameter to generate a uniform distribution of R0 across different trials. `R0` is scaled by `[1-R0Spread,1+R0Spread].`
##### Continuous disease dynamics parameters
`phi` is a factor that scales `R0` to produce an effective reproduction number at different time steps. `phi` responds to changes in the trial-wide prevalence of the disease. The mechanism is as follows:

