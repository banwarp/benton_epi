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
In order to make sure that the script runs within a reasonable amount of time, an upper limit is placed on the number of tirals and the number of nodes. The limits can be changed in the code.
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
