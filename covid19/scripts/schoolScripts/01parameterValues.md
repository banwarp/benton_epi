### Description of parameter values for Benton County Covid-19 Simulator - Education Branch

Choosing the right parameter values is critical for producing reasonable simulations of disease dynamics of COVID-19 in schools. Please review these parameters and make suggestions about how they can be improved. The Benton County Covid-19 Simulator - Education Branch is a modified SEIR compartmental model

#### Last updated: 07/30/2020

### Disease dynamics parameters
- Estimated proportion of community cases that are known to public health = 0.2
  - Source: NEED SOURCE

- Baseline Reproductive Number: R0.
  - R0 = 2.5
    - Source: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box, Scenario 5
  - Alternate base R0 = 2.5 * .021 / .016 = 3.28; Source: OHA modeling used .021 for the pre-intervention beta instead of the .016 that IDM used; therefore 3.28 inflates the base 2.5 R0 by the same factor.
  
- Relative infectiousness of presymptomatic or asymptomatic infections
  - Relative infectiousness of asymptomatic infection = 0.75
  - Relative infectiousness of pre-symptomatic infection is assumed to be the same.
    - Source: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box, Scenario 5
  
- Transforming relative infectiousness into transmission probabilities
  - Proportionate transmissibility of pre- and asymptomatic infections = 0.75 / (0.75 + 1 + 0.1) = 0.4
  - 0.75 is relative infectiousness of pre- and asymptomatic infections; 1 is relative infectiousness of symptomatic infection; 0.1 is relative infectiousness of post-symptomatic infection
    - Source:  https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box, Scenario 5
    - Note: Post-symptomatic infectiousness = 0.1 does not come from literature; it is a baseline guess.
  - Proportionate transmissibility of symptomatic infection = 1 /  (0.75 + 1 + 0.1) = 0.55
  - Proportionate transmissibility of post-symptomatic infection = 0.1 /  (0.75 + 1 + 0.1) = 0.05
  
- Transforming proportionate transmissibility into differential reproductive numbers
  - R0 for pre-symptomatic and asymptomatic infections = 0.4 * 2.5 = 1
  - R0 for symptomatic infections - 0.55 * 2.5 = 1.375
  - R0 for post-symptomatic infections = .05 * 2.5 = 0.125

- Duration parameters
  - I settled on these duration parameters from a combination of the Covasim duration parameters and some educated guesses.
  - Length of pre-symptomatic period = 2 days (directly from Covasim)
    - Source: https://covid.idmod.org/data/Covasim_model_report.pdf
  - Length of symptomatic period = 4 days (Covasim recovery time for mild cases = 6-10 days from symptom onset. I assume that symptoms decrease (but not to zero) by day 4, moving the individual to the post-symptomatic compartment, where symptoms continue to decrease).
    - Source: https://covid.idmod.org/data/Covasim_model_report.pdf
  - Length of post-symptomatic period = 6 days (Covasim recovery time for mild cases = 6-10 days from symptom onset. I assume that symptoms have decreased significantly after 4 days, then continue to decrease for another 6 days).
    - Source: https://covid.idmod.org/data/Covasim_model_report.pdf
  - Length of asymptomatic period (after pre-symptomatic period) = 6 days (Covasim recovery time for asymptomatic = 8 days; subtract 2 days for the presymptomatic period)
    - Source: https://covid.idmod.org/data/Covasim_model_report.pdf
  - Length of period of temporary immunity (after recovery) = 100 days (Just a baseline guess to allow for some reinfection)
  - Probability of moving into permanent immunity status after the temporary immunity compartment = 0.9 (Just a baseline guess to allow for some reinfection)

- Probability of symptomatic infections
  - Probability of symptomatic infections in middle and high school = 0.43
    - Source: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.29.2001352#html_fulltext
  - Probability of symptomatic infections in elementary school = 0.43
    - Source: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.29.2001352#html_fulltext; not direct support, this just is a Bayesian starting point
  - Probability of symptomatic infections in adults (teachers/staff) = 0.65
    - Source: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html#box, Scenario 5

- Differential transmission rates (transmission from student-student; student-teacher; teacher-student; teacher-teacher)
  - Middle and high-school students are as likely as adults to transmit; young children are 1/3 as likely to transmit:
  - Student-student and student-teacher transmission rates factor for elementary = 0.33
  - Student-student and student-teacher transmission rates factor for middle/high school = 1
  - Teacher-teacher transmission rate factor = 1
  - Teacher-student transmission rate factor for elementary = .5, based on moderate levels of distancing and ACE2 reduction
  - Teacher-student transmission rate factor for middle/high = .5, based on high levels of distancing
    - Source: https://wwwnc.cdc.gov/eid/article/26/10/20-1315_article
  - Alternative: Children are 44% as likely as adults to transmit the coronavirus to household contacts:
  - Student-student and student-teacher transmission rates factor = 0.44
  - Teacher-teacher transmission rate factor = 1
  - Teacher- student transmission rate factor = ?
    - Source: https://www.medrxiv.org/content/10.1101/2020.05.20.20108126v1

- Seasonality of the coronavirus
  - Assume that the biological infectiousness of the coronavirus intrisically varies with temperature/humidity, according to seasonality like influenza.
  - The seasonality factor for the coronavirus = 0.125 (represents a 28% higher R0 in winter versus summer)
    - Stephen M Kissler, Christine Tedijanto, Marc Lipsitch, and Yonatan Grad. n.d. “Social Distancing Strategies for Curbing the COVID-19 Epidemic | medRxiv.” Accessed July 9, 2020. https://www.medrxiv.org/content/10.1101/2020.03.22.20041079v1.
 
 - Probability of detecting a coronavirus infection
   - Assume that the probability of detecting an infection is different depending on pre-, post, a- or symptomatic infection
   - Assume that detection efforts are limited to symptom monitoring (see list of child symptoms: https://www.cdc.gov/coronavirus/2019-ncov/symptoms-testing/symptoms.html)
   - These probabilities have to take into account the limited nature of the model; there is no scale of symptom intensity for an individual in the symptomatic compartment; some individuals may have a symptomatic level of transmission probability but unnoticeable mild symptoms; others may have more obvious symptoms; therefore the probabilities are discounted. Likewise, some individuals who are shedding at an asymptomatic rate may still have noticeable symptoms and are detected.
   - Probability of detecting COVID-19 in pre-symptomatic individual = 0.1
   - Probability of detecting COVID-19 in symptomatic individual = 0.8
   - Probability of detecting COVID-19 in post-symptomatic individual = 0.1
   - Probability of detecting COVID-19 in a-symptomatic individual = 0.1
     - Source: Just an educated guess.
   
- Factors that reduce R0 across the board
  - Assume that outside of school, people stay mostly isolated; which reduces mixing, which reduces epidemic pressure, and therefore reduces R0 for every group.
  - Proportionate reduction in baseline R0 during out-of-school hours (afternoons, weekends) = 0.2
  - Proportionate reduction in baseline R0 during nights = 0.05 (Not 0 because there could be some mixing during the night)
  - Proportionate reduction in baseline R0 during classroom/grade/school quarantines = 0.1 (Assumes lower than normal afternoon/weekend because people are being more cautious).
    - Source: Just an educated guess.
  - Proportionate reduction in baseline R0 for high-compliance with face coverings and physical distancing in schools = 0.2
    - Source: 3.1%/17.4% ~= 0.2; Chu, Derek K., Elie A. Akl, Stephanie Duda, Karla Solo, Sally Yaacoub, Holger J. Schünemann, Derek K. Chu, et al. 2020. “Physical Distancing, Face Masks, and Eye Protection to Prevent Person-to-Person Transmission of SARS-CoV-2 and COVID-19: A Systematic Review and Meta-Analysis.” The Lancet 395 (10242): 1973–87. https://doi.org/10.1016/S0140-6736(20)31142-9. 
  - Proportionate reduction in baseline R0 for adequate ventilation in schools/classrooms = 0.33
    - Source = 13.5/37 ~= .33; Table E.1 from https://apps.who.int/iris/bitstream/handle/10665/44167/9789241547857_eng.pdf?sequence=1
    



##### Acknowledgements:
I'd like to thank the developers of the SimInf package for writing the code that these scripts use and for assisting me in adapting their code to meet the requirements of COVID-19 modeling.  

Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R Package for Data-Driven Stochastic Disease Spread Simulations. Journal of Statistical Software, 91(12), 1--42. doi: 10.18637/jss.v091.i12  

Bauer P, Engblom S, Widgren S (2016) Fast event-based epidemiological simulations on national scales. International Journal of High Performance Computing Applications, 30(4), 438--453. doi: 10.1177/1094342016635723

[SimInf vignettes](https://cran.r-project.org/web/packages/SimInf/vignettes/SimInf.pdf)  
[SimInf technical documentation](https://cran.r-project.org/web/packages/SimInf/SimInf.pdf)  
[SimInf git repository](https://github.com/stewid/SimInf)
