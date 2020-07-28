### Description of parameter values for Benton County Covid-19 Simulator - Education Branch

Choosing the right parameter values is critical for producing reasonable simulations of disease dynamics of COVID-19 in schools. Please review these parameters and make suggestions about how they can be improved. The Benton County Covid-19 Simulator - Education Branch is a modified SEIR compartmental model

#### Last updated: 07/28/2020

### Disease dynamics parameters
- Baseline Reproductive Number: R0. 
'''
R0Base = 2.5,                            # Baseline R0
  R0preFrac = .4,                          # Presymptomatic R0 fraction (what proportion of new cases come from presymptomatic)
  R0sympFrac = .55,                          # Symptomatic R0 fraction
  R0postFrac = .05,                         # Postsymptomatic R0 fraction
  R0asympFrac = .4,                        # Asymptomatic R0 fraction
  studentTeacherDiff = rep(1,4),           # Differential R0 for stationary-stationary, transitory-transitory, and stationary-transitory interactions
  preSympPeriod = 1/3,                     # Reciprocal of presymptomatic period
  symptomaticPeriod = 1/4,                 # Reciprocal of symptomatic period
  postSympPeriod = 1/6,                    # Reciprocal of post-symptomatic period
  aSympPeriod = 1/6,                       # Reciprocal of asymptomatic period
  isoPeriod = 1/10,                        # Reciprocal of isolation period
  reSuscepRate = .1,                       # Proportion of recovereds who eventually become susceptible again
  tempImmPeriod = 1/100,                   # Reciprocal of temporary immunity period, after which R becomes Im or S
'''


##### Acknowledgements:
I'd like to thank the developers of the SimInf package for writing the code that these scripts use and for assisting me in adapting their code to meet the requirements of COVID-19 modeling.  

Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12  

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales. International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

[SimInf vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf)  
[SimInf technical documentation](https://cran.r-project.org/web/packages/SimInf/SimInf.pdf)  
[SimInf git repository](https://github.com/stewid/SimInf)
