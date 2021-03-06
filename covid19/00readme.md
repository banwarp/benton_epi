These folders contains scripts for the SimInf implementation of COVID-19 modeling.

As of 06/26/2020, there are four sets of SimInf scripts:
1. Policy scripts - the original implementation of the SimInf that includes policy interventions to reduce R0 in response to prevalence, and includes parachuting infections, super-spreader events, and the option for a mass-entry event.
2. High-low scripts - a modified version of the basic scripts that allows for infectious individuals to be high-spreaders or low-spreaders
3. Age group scripts - a modified version of the high-low scripts that uses four age groups: 0-17, 18-29, 30-64, and 65+. Each age group can have a different R0 and different levels of complexity of COVID-19 illness (i.e. hospitalization and fatality rates).
4. School scripts - a reworked version of the original implementation for modeling disease dynamics in schools. This version does not have parachuting infections or mass-entry. Instead of policy responses, this version uses different quarantine schemes to control spread of disease: Individual isolation, classroom quarantine, grade quarantine, or school closure. The school scripts are different enough from the others that they have their own set of readmes in the scripts/schoolScripts folder.

SimInf is a powerful R package that implements disease dynamics in C using a stochastic approach. It runs very quickly and has the ability to incorporate user-defined transitions. The disease dynamics are built on standard SIR, SEIR, and similar models.

The COVID-19 scripts in this repository pre- and post-process data to adapt the SimInf functions to disease dynamics characteristic of the COVID-19 pandemic.

At this point, the scripts are not packaged, because I have not done that before and I am not sure the best way to do it. The scripts are in their raw form and you are welcome to download and adapt them to your purposes.

I recommend that you begin by reviewing the Policy Scripts folder, as this has the most comprehensive set of readmes for the scripts. If you review these readmes thoroughly, I think the other scripts (hilow, ageGroup, and school) will be clear enough from the commented code.

Briefly:
The model is a modified SEIR. Some proportion of Recovered lose their Immunity and return to Susceptible. Some proportion of Exposed are identified through contact tracing and are Isolated. Some proportion of Infectious are Hospitalized. Some proportion of Hospitalized (optional: also Infectious) die of COVID-19. I call this model SEIR-(S,Is,Im,H,M).
The script sets initial conditions, disease dynamic parameters, and a series of events. The events include parachuting infected individuals, routine transfers between sub-populations, and optional mass events (e.g. university students returning to a college town or super spreader).
When the model runs, the script used the real-time prevalence to automatically introduce and lift policy interventions (major and minor), and allows for continual or temporary physical distancing. In Oregon, the current stay-at-home orders will be lifted in 3 phases, so these events are included. Also, in order to explore counterfactuals, there is an option to "switch off" policy interventions after a certain date.
Each simulation produces a user-defined number of trials to generate a confidence or credibility spread.
The SimInf package produces trajectories of the different compartments (e.g. Infectious or Recovered) and of the continuous variables (e.g. the intensity of the policy intervention). The script can plot these different trajectories with their credibility spreads.

You are welcome to leave comments/questions on github or email me at peter.banwarth@co.benton.or.us.

Citation:
Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic
Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales.
International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf
https://cran.r-project.org/web/packages/SimInf/SimInf.pdf
https://github.com/stewid/SimInf
