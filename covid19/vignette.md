## A basic vignette to describe how to use the covid19 set of R scripts
### Introduction
This set of R scripts was developed for local modeling of the COVID-19 disease in Benton County, Oregon, and other small geographies. The scripts use a modified SEIR model with partially-mixed, but otherwise homogeneous population nodes. This vignette describes how to use the scripts.

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

