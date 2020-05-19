This folder contains scripts for the SimInf implementation of COVID-19 modeling.

SimInf is a powerful R package that implements disease dynamics in C using a stochastic approach. It runs very quickly and has the ability to incorporate user-defined transitions. The disease dynamics are built on standard SIR, SEIR, and similar models.

The COVID-19 scripts in this repository pre- and post-process data to adapt the SimInf functions to disease dynamics characteristic of the COVID-19 pandemic.

At this point, the scripts are not packaged, because I have not done that before and I am not sure the best way to do it. The scripts are in their raw form and you are welcome to download and adapt them to your purposes.

I recommend that you begin by reviewing the file BentonCountycovidProjections2020.05.12-0580.pdf. It discusses the SimInf model and the different components of the COVID-19 script.

Briefly:
The model is a modified SEIR. Some proportion of Recovered lose their Immunity and return to Susceptible. Some proportion of Exposed are identified through contact tracing and are Isolated. Some proportion of Infectious are Hospitalized. Some proportion of Hospitalized (optional: also Infectious) die of COVID-19. I call this model SEIR-(S,Is,Im,H,M).
The script sets initial conditions, disease dynamic parameters, and a series of events. The events include parachuting infected individuals, routine transfers between sub-populations, and optional mass events (e.g. university students returning to a college town or super spreader).
When the model runs, the script used the real-time prevalence to automatically introduce and lift policy interventions (major and minor), and allows for continual or temporary physical distancing. In Oregon, the current stay-at-home orders will be lifted in 3 phases, so these events are included. Also, in order to explore counterfactuals, there is an option to "switch off" policy interventions after a certain date.
Each simulation produces a user-defined number of trials to generate a confidence or credibility spread.
The SimInf package produces trajectories of the different compartments (e.g. Infectious or Recovered) and of the continuous variables (e.g. the intensity of the policy intervention). The script can plot these different trajectories with their credibility spreads.

To use the scripts:
1. Download them to your machine, then open covid19_SimInf_5.12.2020.R and review the script.
2. Change any folder settings as needed in the script: search for "setwd(".
3. Open covidExampleScript_5.12.2020.R for an example of how to produce projections for all Oregon counties under three scenarios.
4. Select your parameters. There are 70+ adjustable parameters in the function. (Not counting the preset parameters in the rest of the code, which can be changed as you tinker). The parameters are explained in the comments. Refer to the pdf file BentonCountycovidProjections2020.05.12-0580.pdf for more detailed descriptions.
5. Unless you are modeling a community with a single large entry event (e.g. university students returning on/around a certain date), you probably don't want to include the studentEvent. Simply set maxStudentNodes = 0 and this event will not be created.
6. You can include zero, one, or more super-spreader events. The parameters for the super spreader events are lists to accommodate more than one event.
7. You have the option to group the nodes to represent segmented populations where movement is higher within node groups compared to across node groups.
8. Run the script. The script automatically produces plots of daily active infections (I), intervention measures (phi), daily new infections (newI), daily hospital bed use (H), and cumulative deaths (M), and saves them to a folder. You can also use the code for these plots to plot other compartments or continuous variables. Read through simInfPlottingFunction5.14.2020.R for details. The function itself returns the modeled result, from which you can extract trajectories and other data and make other plots.

You are welcome to leave comments/questions on github or email me at peter.banwarth@co.benton.or.us.

One consequence of this semi-mixed approach is that the number of nodes and the transfer parameters strongly impact the disease dynamics. Highly segmented (many nodes/low transfers) populations contain mini-epidemics, leading to lower overall infections, while highly-mixed populations (few nodes/high transfers) experience wide-spread epidemics.

Citation:
Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
https://github.com/stewid/SimInf