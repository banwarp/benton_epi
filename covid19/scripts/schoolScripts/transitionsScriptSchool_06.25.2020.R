# transitionsScriptSchool_06.25.2020.R

transitionsGroupQuarantine <- function() {
  
  transitions <-  c(
    "Rs -> reSuscepRate*tempImmPeriod*Rs -> Ss",                     # Recovereds who become susceptible again: 
    paste0("Ss ->",
           "outOfSchool*",
           "(ssD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "stD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*Ss/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSymps + cumIs"
    ),                                                               # Infection: 
    "preSymps -> symptomaticProp*preSympPeriod*preSymps -> Symps",        # Development of symptomatic disease among exposed: stationary
    "preSymps -> (1-symptomaticProp)*preSympPeriod*preSymps -> aSymps",     # Development of asymptomatic disease among exposed: stationary
    "Symps -> symptomaticPeriod*Symps -> postSymps",                             # Transition from symptomatic disease to recovered 
    "postSymps -> postSympPeriod*postSymps -> Rs",                # Transition from post symptomatic to recovered
    "aSymps -> aSympPeriod*aSymps -> Rs",                            # Transition from asymptomatic disease to recovered: 
    "Isos -> isoPeriod*Isos -> Rs",                                # Transition from isolated to recovered
    "Rs -> (1-reSuscepRate)*tempImmPeriod*Rs -> Ims",                # Development of immunity: 
    
    "Rt -> reSuscepRate*tempImmPeriod*Rt -> St",                     # Recovereds who become susceptible again: 
    paste0("St ->",
           "outOfSchool*",
           "(stD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "ttD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*St/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSympt + cumIt"
    ),                                                               # Infection: 
    "preSympt -> symptomaticProp*preSympPeriod*preSympt -> Sympt",        # Development of symptomatic disease among exposed: stationary
    "preSympt -> (1-symptomaticProp)*preSympPeriod*preSympt -> aSympt",     # Development of asymptomatic disease among exposed: stationary
    "Sympt -> symptomaticPeriod*Sympt -> postSympt",                             # Transition from symptomatic disease to recovered 
    "postSympt -> postSympPeriod*postSympt -> Rt",                 # Transition from post symptomatic to recovered
    "aSympt -> symptomaticPeriod*aSympt -> Rt",                            # Transition from asymptomatic disease to recovered: 
    "Isot -> isoPeriod*Isot -> Rt",                                # Transition from isolated to recovered
    "Rt -> (1-reSuscepRate)*tempImmPeriod*Rt -> Imt"                # Development of immunity: 
  )
  
  return(transitions)
}

transitionsIndividualIsolation <- function() {
  transitions <-  c(
    "Rs -> reSuscepRate*tempImmPeriod*Rs -> Ss",                     # Recovereds who become susceptible again: 
    paste0("Ss ->",
           "outOfSchool*",
           "(ssD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "stD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*Ss/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSymps + cumIs"
    ),                                                               # Infection: 
    "preSymps -> preDetectSuccess*preSympPeriod*preSymps -> Isos",        # Identification and Isolation of presymptomatic
    "preSymps -> (1-preDetectSuccess)*symptomaticProp*preSympPeriod*preSymps -> Symps",        # Development of symptomatic disease among exposed: stationary
    "preSymps -> (1-symptomaticProp)*preSympPeriod*preSymps -> aSymps",     # Development of asymptomatic disease among exposed: stationary
    "Symps -> sympDetectSuccess*symptomaticPeriod*Symps -> Isos",                             # Identification and Isolation of symptomatic
    "Symps -> (1-sympDetectSuccess)*symptomaticPeriod*Symps -> postSymps",                             # Transition from symptomatic disease to recovered 
    "postSymps -> postDetectSuccess*postSympPeriod*postSymps -> Isos",                # Identification and Isolation of postsymptomatic
    "postSymps -> (1-postDetectSuccess)*postSympPeriod*postSymps -> Rs",                # Transition from post symptomatic to recovered
    "aSymps -> asympDetectSuccess*aSympPeriod*aSymps -> Isos",                            # Identification and Isolation of asymptomatic
    "aSymps -> (1-asympDetectSuccess)*aSympPeriod*aSymps -> Rs",                            # Transition from asymptomatic disease to recovered: 
    "Isos -> isoPeriod*Isos -> Rs",                                # Transition from isolated to recovered
    "Rs -> (1-reSuscepRate)*tempImmPeriod*Rs -> Ims",                # Development of immunity: 
    
    "Rt -> reSuscepRate*tempImmPeriod*Rt -> St",                     # Recovereds who become susceptible again: 
    paste0("St ->",
           "outOfSchool*",
           "(stD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "ttD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*St/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSympt + cumIt"
    ),                                                               # Infection: 
    "preSympt -> preDetectSuccess*preSympPeriod*preSympt -> Isot",        # Identification and Isolation of presymptomatic
    "preSympt -> (1-preDetectSuccess)*symptomaticProp*preSympPeriod*preSympt -> Sympt",        # Development of symptomatic disease among exposed: stationary
    "preSympt -> (1-symptomaticProp)*preSympPeriod*preSympt -> aSympt",     # Development of asymptomatic disease among exposed: stationary
    "Sympt -> sympDetectSuccess*symptomaticPeriod*Sympt -> Isot",                             # Identification and Isolation of symptomatic
    "Sympt -> (1-sympDetectSuccess)*symptomaticPeriod*Sympt -> postSympt",                             # Transition from symptomatic disease to recovered 
    "postSympt -> postDetectSuccess*postSympPeriod*postSympt -> Isot",                # Identification and Isolation of postsymptomatic
    "postSympt -> (1-postDetectSuccess)*postSympPeriod*postSympt -> Rt",                # Transition from post symptomatic to recovered
    "aSympt -> asympDetectSuccess*aSympPeriod*aSympt -> Isot",                            # Identification and Isolation of asymptomatic
    "aSympt -> (1-asympDetectSuccess)*aSympPeriod*aSympt -> Rt",                            # Transition from asymptomatic disease to recovered: 
    "Isot -> isoPeriod*Isot -> Rt",                                # Transition from isolated to recovered
    "Rt -> (1-reSuscepRate)*tempImmPeriod*Rt -> Imt"                # Development of immunity: 
  )
  
  return(transitions)
}