# transitionsScriptAgeGroups06.19.2020.R

# script to hold the transitions for the 4 age group SimInf since it is so big now.

transitionsScript <- function() {
  transitions <-  c(
    "Rk -> reSuscepRate*tempImmPeriod*Rk -> Sk",                           # Recovereds who become susceptible again
    paste0("Sk -> ",
           "(kkf*(betaI*(1/phi)*(spreadRatio*hIk+(1/spreadRatio)*lIk)+",
           "betaU*(1/phi)*(spreadRatio*uhIk+(1/spreadRatio)*ulIk)+",
           "betaP*pIk)+",
           "kyf*(betaI*(1/phi)*(spreadRatio*hIy+(1/spreadRatio)*lIy)+",
           "betaU*(1/phi)*(spreadRatio*uhIy+(1/spreadRatio)*ulIy)+",
           "betaP*pIy)+",
           "kaf*(betaI*(1/phi)*(spreadRatio*hIa+(1/spreadRatio)*lIa)+",
           "betaU*(1/phi)*(spreadRatio*uhIa+(1/spreadRatio)*ulIa)+",
           "betaP*pIa)+",
           "ksf*(betaI*(1/phi)*(spreadRatio*hIs+(1/spreadRatio)*lIs)+",
           "betaU*(1/phi)*(spreadRatio*uhIs+(1/spreadRatio)*ulIs)+",
           "betaP*pIs))*",
           "season*betaRandomizer*Sk/",
           "(Sk+Sy+Sa+Ss+Ek+Ey+Ea+Es+hIk+hIy+hIa+hIs+lIk+lIy+lIa+lIs+uhIk+uhIy+uhIa+uhIs+ulIk+ulIy+ulIa+ulIs+Rk+Ry+Ra+Rs+Imk+Imy+Ima+Ims)",
           " -> Ek"
    ),
    "Ek -> highSpreadProp*(1-isoK)*exposedPeriod*Ek -> hIk + cumIk",     # Development of disease among non-isolated exposed high spreaders
    "Ek -> (1-highSpreadProp)*(1-isoK)*exposedPeriod*Ek -> lIk + cumIk", # Development of disease among non-isolated exposed low spreaders
    "Ek -> (1-hospExpK)*isoK*exposedPeriod*Ek -> pIk + cumIk",        # Isolation of exposed (non-hospitalized)
    "Ek -> hospExpK*isoK*exposedPeriod*Ek -> Hk + cumIk",                       # New hospitalizations from exposed
    "hIk -> mSk*(1-hospPostK)*initInfectiousPeriod*hIk -> pIk",              # Transition from high spreader Infectious to Post-Infectious because of symptom monitoring
    "lIk -> mSk*(1-hospPostK)*initInfectiousPeriod*lIk -> pIk",              # Transition from low spreader Infectious to Post-Infectious because of symptom monitoring
    "hIk -> mSk*hospPostK*initInfectiousPeriod*hIk -> Hk",
    "lIk -> mSk*hospPostK*initInfectiousPeriod*lIk -> Hk",
    "hIk -> (1-mSk)*initInfectiousPeriod*hIk -> uhIk",              # Transition to high spread unknown-infectious if not identified by symptoms/testing
    "lIk -> (1-mSk)*initInfectiousPeriod*lIk -> ulIk",              # Transition to low spread unknown-infectious if not identified by symptoms/testing
    "pIk -> (1-nonHospDeathK)*postInfectiousPeriod*pIk -> Rk",          # Recovery among post-infectious
    "pIk -> nonHospDeathK*postInfectiousPeriod*pIk -> Mk",              # Deaths among non-hospitalized COVID patients (post-infectious)
    "uhIk -> (1-nonHospDeathK)*unknownInfectiousPeriod*uhIk -> Rk",     # Recovery among high spread unknown-infectious
    "ulIk -> (1-nonHospDeathK)*unknownInfectiousPeriod*ulIk -> Rk",     # Recovery among low spread unknown-infectious
    "uhIk -> nonHospDeathK*unknownInfectiousPeriod*uhIk -> Mk",         # Deaths among non-hospitalized COVID patients (high spread unknown-infectious)
    "ulIk -> nonHospDeathK*unknownInfectiousPeriod*ulIk -> Mk",         # Deaths among non-hospitalized COVID patients (low spread unknown-infectious)
    "Hk -> (1-hospDeathK)*hospPeriod*Hk -> Rk",                         # Recovery among hospitalized
    "Hk -> hospDeathK*Hk -> Mk",                                        # Deaths among hospitalized COVID patients
    "Rk -> (1-reSuscepRate)*tempImmPeriod*Rk -> Imk",                      # Development of immunity
    
    "Ry -> reSuscepRate*tempImmPeriod*Ry -> Sy",                           # Recovereds who become susceptible again
    paste0("Sy -> ",
           "(kyf*(betaI*(1/phi)*(spreadRatio*hIk+(1/spreadRatio)*lIk)+",
           "betaU*(1/phi)*(spreadRatio*uhIk+(1/spreadRatio)*ulIk)+",
           "betaP*pIk)+",
           "yyf*(betaI*(1/phi)*(spreadRatio*hIy+(1/spreadRatio)*lIy)+",
           "betaU*(1/phi)*(spreadRatio*uhIy+(1/spreadRatio)*ulIy)+",
           "betaP*pIy)+",
           "yaf*(betaI*(1/phi)*(spreadRatio*hIa+(1/spreadRatio)*lIa)+",
           "betaU*(1/phi)*(spreadRatio*uhIa+(1/spreadRatio)*ulIa)+",
           "betaP*pIa)+",
           "ysf*(betaI*(1/phi)*(spreadRatio*hIs+(1/spreadRatio)*lIs)+",
           "betaU*(1/phi)*(spreadRatio*uhIs+(1/spreadRatio)*ulIs)+",
           "betaP*pIs))*",
           "season*betaRandomizer*Sy/",
           "(Sk+Sy+Sa+Ss+Ek+Ey+Ea+Es+hIk+hIy+hIa+hIs+lIk+lIy+lIa+lIs+uhIk+uhIy+uhIa+uhIs+ulIk+ulIy+ulIa+ulIs+Rk+Ry+Ra+Rs+Imk+Imy+Ima+Ims)",
           " -> Ey"
    ),
    "Ey -> highSpreadProp*(1-isoY)*exposedPeriod*Ey -> hIy + cumIy",     # Development of disease among non-isolated exposed high spreaders
    "Ey -> (1-highSpreadProp)*(1-isoY)*exposedPeriod*Ey -> lIy + cumIy", # Development of disease among non-isolated exposed low spreaders
    "Ey -> (1-hospExpY)*isoY*exposedPeriod*Ey -> pIy + cumIy",        # Isolation of exposed (non-hospitalized)
    "Ey -> hospExpY*isoY*exposedPeriod*Ey -> Hy + cumIy",                       # New hospitalizations from exposed
    "hIy -> mSy*(1-hospPostY)*initInfectiousPeriod*hIy -> pIy",              # Transition from high spreader Infectious to Post-Infectious because of symptom monitoring
    "lIy -> mSy*(1-hospPostY)*initInfectiousPeriod*lIy -> pIy",              # Transition from low spreader Infectious to Post-Infectious because of symptom monitoring
    "hIy -> mSy*hospPostY*initInfectiousPeriod*hIy -> Hy",
    "lIy -> mSy*hospPostY*initInfectiousPeriod*lIy -> Hy",
    "hIy -> (1-mSy)*initInfectiousPeriod*hIy -> uhIy",              # Transition to high spread unknown-infectious if not identified by symptoms/testing
    "lIy -> (1-mSy)*initInfectiousPeriod*lIy -> ulIy",              # Transition to low spread unknown-infectious if not identified by symptoms/testing
    "pIy -> (1-nonHospDeathY)*postInfectiousPeriod*pIy -> Ry",          # Recovery among post-infectious
    "pIy -> nonHospDeathY*postInfectiousPeriod*pIy -> My",              # Deaths among non-hospitalized COVID patients (post-infectious)
    "uhIy -> (1-nonHospDeathY)*unknownInfectiousPeriod*uhIy -> Ry",     # Recovery among high spread unknown-infectious
    "ulIy -> (1-nonHospDeathY)*unknownInfectiousPeriod*ulIy -> Ry",     # Recovery among low spread unknown-infectious
    "uhIy -> nonHospDeathY*unknownInfectiousPeriod*uhIy -> My",         # Deaths among non-hospitalized COVID patients (high spread unknown-infectious)
    "ulIy -> nonHospDeathY*unknownInfectiousPeriod*ulIy -> My",         # Deaths among non-hospitalized COVID patients (low spread unknown-infectious)
    "Hy -> (1-hospDeathY)*hospPeriod*Hy -> Ry",                         # Recovery among hospitalized
    "Hy -> hospDeathY*Hy -> My",                                        # Deaths among hospitalized COVID patients
    "Ry -> (1-reSuscepRate)*tempImmPeriod*Ry -> Imy",                      # Development of immunity
    
    "Ra -> reSuscepRate*tempImmPeriod*Ra -> Sa",                           # Recovereds who become susceptible again
    paste0("Sa -> ",
           "(kaf*(betaI*(1/phi)*(spreadRatio*hIk+(1/spreadRatio)*lIk)+",
           "betaU*(1/phi)*(spreadRatio*uhIk+(1/spreadRatio)*ulIk)+",
           "betaP*pIk)+",
           "yaf*(betaI*(1/phi)*(spreadRatio*hIy+(1/spreadRatio)*lIy)+",
           "betaU*(1/phi)*(spreadRatio*uhIy+(1/spreadRatio)*ulIy)+",
           "betaP*pIy)+",
           "aaf*(betaI*(1/phi)*(spreadRatio*hIa+(1/spreadRatio)*lIa)+",
           "betaU*(1/phi)*(spreadRatio*uhIa+(1/spreadRatio)*ulIa)+",
           "betaP*pIa)+",
           "asf*(betaI*(1/phi)*(spreadRatio*hIs+(1/spreadRatio)*lIs)+",
           "betaU*(1/phi)*(spreadRatio*uhIs+(1/spreadRatio)*ulIs)+",
           "betaP*pIs))*",
           "season*betaRandomizer*Sa/",
           "(Sk+Sy+Sa+Ss+Ek+Ey+Ea+Es+hIk+hIy+hIa+hIs+lIk+lIy+lIa+lIs+uhIk+uhIy+uhIa+uhIs+ulIk+ulIy+ulIa+ulIs+Rk+Ry+Ra+Rs+Imk+Imy+Ima+Ims)",
           " -> Ea"
    ),
    "Ea -> highSpreadProp*(1-isoA)*exposedPeriod*Ea -> hIa + cumIa",     # Development of disease among non-isolated exposed high spreaders
    "Ea -> (1-highSpreadProp)*(1-isoA)*exposedPeriod*Ea -> lIa + cumIa", # Development of disease among non-isolated exposed low spreaders
    "Ea -> (1-hospExpA)*isoA*exposedPeriod*Ea -> pIa + cumIa",        # Isolation of exposed (non-hospitalized)
    "Ea -> hospExpA*isoA*exposedPeriod*Ea -> Ha + cumIa",                       # New hospitalizations from exposed
    "hIa -> mSa*(1-hospPostA)*initInfectiousPeriod*hIa -> pIa",              # Transition from high spreader Infectious to Post-Infectious because of symptom monitoring
    "lIa -> mSa*(1-hospPostA)*initInfectiousPeriod*lIa -> pIa",              # Transition from low spreader Infectious to Post-Infectious because of symptom monitoring
    "hIa -> mSa*hospPostA*initInfectiousPeriod*hIa -> Ha",
    "lIa -> mSa*hospPostA*initInfectiousPeriod*lIa -> Ha",
    "hIa -> (1-mSa)*initInfectiousPeriod*hIa -> uhIa",              # Transition to high spread unknown-infectious if not identified by symptoms/testing
    "lIa -> (1-mSa)*initInfectiousPeriod*lIa -> ulIa",              # Transition to low spread unknown-infectious if not identified by symptoms/testing
    "pIa -> (1-nonHospDeathA)*postInfectiousPeriod*pIa -> Ra",          # Recovery among post-infectious
    "pIa -> nonHospDeathA*postInfectiousPeriod*pIa -> Ma",              # Deaths among non-hospitalized COVID patients (post-infectious)
    "uhIa -> (1-nonHospDeathA)*unknownInfectiousPeriod*uhIa -> Ra",     # Recovera among high spread unknown-infectious
    "ulIa -> (1-nonHospDeathA)*unknownInfectiousPeriod*ulIa -> Ra",     # Recovera among low spread unknown-infectious
    "uhIa -> nonHospDeathA*unknownInfectiousPeriod*uhIa -> Ma",         # Deaths among non-hospitalized COVID patients (high spread unknown-infectious)
    "ulIa -> nonHospDeathA*unknownInfectiousPeriod*ulIa -> Ma",         # Deaths among non-hospitalized COVID patients (low spread unknown-infectious)
    "Ha -> (1-hospDeathA)*hospPeriod*Ha -> Ra",                         # Recovera among hospitalized
    "Ha -> hospDeathA*Ha -> Ma",                                        # Deaths among hospitalized COVID patients
    "Ra -> (1-reSuscepRate)*tempImmPeriod*Ra -> Ima",                      # Development of immunity
    
    "Rs -> reSuscepRate*tempImmPeriod*Rs -> Ss",                           # Recovereds who become susceptible again
    paste0("Ss -> ",
           "(ksf*(betaI*(1/phi)*(spreadRatio*hIk+(1/spreadRatio)*lIk)+",
           "betaU*(1/phi)*(spreadRatio*uhIk+(1/spreadRatio)*ulIk)+",
           "betaP*pIk)+",
           "ysf*(betaI*(1/phi)*(spreadRatio*hIy+(1/spreadRatio)*lIy)+",
           "betaU*(1/phi)*(spreadRatio*uhIy+(1/spreadRatio)*ulIy)+",
           "betaP*pIy)+",
           "asf*(betaI*(1/phi)*(spreadRatio*hIa+(1/spreadRatio)*lIa)+",
           "betaU*(1/phi)*(spreadRatio*uhIa+(1/spreadRatio)*ulIa)+",
           "betaP*pIa)+",
           "ssf*(betaI*(1/phi)*(spreadRatio*hIs+(1/spreadRatio)*lIs)+",
           "betaU*(1/phi)*(spreadRatio*uhIs+(1/spreadRatio)*ulIs)+",
           "betaP*pIs))*",
           "season*betaRandomizer*Ss/",
           "(Sk+Sy+Sa+Ss+Ek+Ey+Ea+Es+hIk+hIy+hIa+hIs+lIk+lIy+lIa+lIs+uhIk+uhIy+uhIa+uhIs+ulIk+ulIy+ulIa+ulIs+Rk+Ry+Ra+Rs+Imk+Imy+Ima+Ims)",
           " -> Es"
    ),
    "Es -> highSpreadProp*(1-isoS)*exposedPeriod*Es -> hIs + cumIs",     # Development of disease among non-isolated exposed high spreaders
    "Es -> (1-highSpreadProp)*(1-isoS)*exposedPeriod*Es -> lIs + cumIs", # Development of disease among non-isolated exposed low spreaders
    "Es -> (1-hospExpS)*isoS*exposedPeriod*Es -> pIs + cumIs",        # Isolation of exposed (non-hospitalized)
    "Es -> hospExpS*isoS*exposedPeriod*Es -> Hs + cumIs",                       # New hospitalizations from exposed
    "hIs -> mSs*(1-hospPostS)*initInfectiousPeriod*hIs -> pIs",              # Transition from high spreader Infectious to Post-Infectious because of symptom monitoring
    "lIs -> mSs*(1-hospPostS)*initInfectiousPeriod*lIs -> pIs",              # Transition from low spreader Infectious to Post-Infectious because of symptom monitoring
    "hIs -> mSs*hospPostS*initInfectiousPeriod*hIs -> Hs",
    "lIs -> mSs*hospPostS*initInfectiousPeriod*lIs -> Hs",
    "hIs -> (1-mSs)*initInfectiousPeriod*hIs -> uhIs",              # Transition to high spread unknown-infectious if not identified by symptoms/testing
    "lIs -> (1-mSs)*initInfectiousPeriod*lIs -> ulIs",              # Transition to low spread unknown-infectious if not identified by symptoms/testing
    "pIs -> (1-nonHospDeathS)*postInfectiousPeriod*pIs -> Rs",          # Recovery among post-infectious
    "pIs -> nonHospDeathS*postInfectiousPeriod*pIs -> Ms",              # Deaths among non-hospitalized COVID patients (post-infectious)
    "uhIs -> (1-nonHospDeathS)*unknownInfectiousPeriod*uhIs -> Rs",     # Recovers among high spread unknown-infectious
    "ulIs -> (1-nonHospDeathS)*unknownInfectiousPeriod*ulIs -> Rs",     # Recovers among low spread unknown-infectious
    "uhIs -> nonHospDeathS*unknownInfectiousPeriod*uhIs -> Ms",         # Deaths among non-hospitalized COVID patients (high spread unknown-infectious)
    "ulIs -> nonHospDeathS*unknownInfectiousPeriod*ulIs -> Ms",         # Deaths among non-hospitalized COVID patients (low spread unknown-infectious)
    "Hs -> (1-hospDeathS)*hospPeriod*Hs -> Rs",                         # Recovers among hospitalized
    "Hs -> hospDeathS*Hs -> Ms",                                        # Deaths among hospitalized COVID patients
    "Rs -> (1-reSuscepRate)*tempImmPeriod*Rs -> Ims"                      # Development of immunity
  )
  
  return(transitions)
}










