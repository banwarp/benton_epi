## COVID-19 Vignettes
#### Introduction
This set of R scripts was developed for local modeling of the COVID-19 disease in Benton County, Oregon, and other small geographies. The scripts use a modified SEIR model with partially-mixed, but otherwise homogeneous population nodes. This markdown document presents a set of vignettes and explains the output.

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

#### Vignette 1: The basics
All of the arguments for the COVID-19 script (`covidWrapper`) are defaulted, so you can just run the script with no arguments:
```
r1 <- covidWrapper()
```

The output for `r1` is the model after it has been run:
```
Model: SimInf_model
Number of nodes: 50000
Number of transitions: 18
Number of scheduled events: 262217

Global data
-----------
 Parameter        Value     
 beta             0.36250000
 infectiousPeriod 0.12500000
 exposedPeriod    0.25000000
 isoRate          0.12500000
 isoPeriod        0.10000000
 hospRate         0.03300000
 hospPeriod       0.07142857
 betaIsolated     0.01250000
 nonHospDeathRate 0.00000000
 hospDeathRate    0.12500000
 reSuscepRate     0.10000000
 tempImmPeriod    0.01000000
 mu               0.00000000
 nu               0.00000000

Local data
----------
                    Min.  1st Qu.   Median     Mean
 firstNode      0.00e+00 1.24e+04 2.48e+04 2.48e+04
 betaRandomizer 9.01e-01 9.48e-01 1.01e+00 1.00e+00
                 3rd Qu.     Max.
 firstNode      3.71e+04 4.95e+04
 betaRandomizer 1.05e+00 1.10e+00

Continuous state variables
--------------------------
                   Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
 phi           1.00e+00 1.20e+00 1.44e+00 1.65e+00 1.93e+00 3.22e+00
 prevUp01      0.00e+00 0.00e+00 0.00e+00 1.30e+00 1.00e+00 1.10e+01
 prevUp12      0.00e+00 0.00e+00 0.00e+00 2.56e-02 0.00e+00 1.10e+01
 prevDown21    0.00e+00 2.90e+01 2.90e+01 2.85e+01 3.00e+01 3.00e+01
 prevDown10    0.00e+00 2.90e+01 2.90e+01 2.66e+01 2.90e+01 2.90e+01
 kbPhase       0.00e+00 3.00e+00 3.00e+00 2.60e+00 3.00e+00 3.00e+00
 policy        1.00e+00 1.00e+00 1.00e+00 1.00e+00 1.00e+00 1.00e+00
 season        7.50e-01 8.26e-01 9.96e-01 1.00e+00 1.17e+00 1.25e+00
 prevalence    0.00e+00 4.59e-04 8.08e-04 7.95e-04 1.07e-03 2.53e-03
 previousState 0.00e+00 0.00e+00 0.00e+00 3.66e-01 1.00e+00 2.00e+00
 RTee          7.12e-01 1.51e+00 2.02e+00 1.98e+00 2.43e+00 3.93e+00
 pdCounter     0.00e+00 0.00e+00 2.00e+00 1.95e+01 2.70e+01 2.91e+02

Compartments
------------
          Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
 S    0.00e+00 9.20e+01 1.25e+02 1.28e+02 1.56e+02 5.45e+02
 E    0.00e+00 0.00e+00 0.00e+00 5.80e-02 0.00e+00 4.40e+01
 I    0.00e+00 0.00e+00 0.00e+00 1.06e-01 0.00e+00 6.80e+01
 R    0.00e+00 0.00e+00 0.00e+00 1.03e+00 0.00e+00 2.06e+02
 Im   0.00e+00 0.00e+00 0.00e+00 1.29e+00 0.00e+00 1.81e+02
 H    0.00e+00 0.00e+00 0.00e+00 2.58e-03 0.00e+00 6.00e+00
 Is   0.00e+00 0.00e+00 0.00e+00 1.29e-02 0.00e+00 1.10e+01
 cumI 0.00e+00 0.00e+00 0.00e+00 2.53e+00 0.00e+00 3.07e+02
 M    0.00e+00 0.00e+00 0.00e+00 5.01e-02 0.00e+00 1.10e+01
 ```
From this model, you can extract compartment trajectories, continous variables, the events data frame, etc. Refer to the [SimInf Vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf) and the [SimInf technical documentation](https://cran.r-project.org/web/packages/SimInf/SimInf.pdf) for details.  

In addition to returning the model, the script also plots 5 trajectories and saves them as .png files.
```
trajPlotInfections # Plot of daily active infections
trajPlotPhi        # Plot of phi, the intervention intensity continuous variable
trajPlotnewI       # Plot of daily new infections
trajPlotHosp       # Plot of daily hospitalizations
trajPlotDeaths     # Plot of cumulative deaths
```
Here is the plot of daily active infections:
![trajPlotInfections](images/trajPlotIGeneric.png)
