# transitionsScriptSchool_10.07.2020.R

# changes from transitionsScriptSchool_07.04.2020.R
# added cumulative isolation compartment for cohort testing

# changes from transitionsScriptSchool_06.25.2020.R
# Consolidated individual and group transitions into one set to match main code

transitionsScript <- function() {
  transitions <-  c(
    "Rs -> reSuscepRate*tempImmPeriod*Rs -> Ss",                     # Recovereds who become susceptible again: 
    paste0("Ss ->",
           "outOfSchool*",
           "R0ComplianceFrac*",
           "R0VentilationFrac*",
           "(ssD*(betaPre*preSymps+",
           "betaSymp*Symps+",
           "betaPost*postSymps+",
           "betaA*aSymps)+",
           "tsD*(betaPre*preSympt+",
           "betaSymp*Sympt+",
           "betaPost*postSympt+",
           "betaA*aSympt))*",
           "season*betaRandomizer*Ss/",
           "(Ss+St+preSymps+preSympt+Symps+Sympt+postSymps+postSympt+aSymps+aSympt+Rs+Rt+Ims+Imt)",
           " -> preSymps + cumIs"
    ),                                                               # Infection: 
    "preSymps -> preDetectSuccess*preSymps -> Isos + cumIso",        # Identification and Isolation of presymptomatic
    "preSymps -> (1-preDetectSuccess)*symptomaticProps*preSympPeriod*preSymps -> Symps",        # Development of symptomatic disease among exposed: stationary
    "preSymps -> (1-preDetectSuccess)*(1-symptomaticProps)*preSympPeriod*preSymps -> aSymps",     # Development of asymptomatic disease among exposed: stationary
    "Symps -> sympDetectSuccess*Symps -> Isos + cumIso",                             # Identification and Isolation of symptomatic
    "Symps -> (1-sympDetectSuccess)*symptomaticPeriod*Symps -> postSymps",                             # Transition from symptomatic disease to recovered 
    "postSymps -> postDetectSuccess*postSymps -> Isos + cumIso",                # Identification and Isolation of postsymptomatic
    "postSymps -> (1-postDetectSuccess)*postSympPeriod*postSymps -> Rs",                # Transition from post symptomatic to recovered
    "aSymps -> asympDetectSuccess*aSymps -> Isos + cumIso",                            # Identification and Isolation of asymptomatic
    "aSymps -> (1-asympDetectSuccess)*aSympPeriod*aSymps -> Rs",                            # Transition from asymptomatic disease to recovered: 
    "Isos -> isoPeriod*Isos -> Rs",                                # Transition from isolated to recovered
    "Rs -> (1-reSuscepRate)*tempImmPeriod*Rs -> Ims",                # Development of immunity: 
    
    "Rt -> reSuscepRate*tempImmPeriod*Rt -> St",                     # Recovereds who become susceptible again: 
    paste0("St ->",
           "outOfSchool*",
           "R0ComplianceFrac*",
           "R0VentilationFrac*",
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
    "preSympt -> preDetectSuccess*preSympt -> Isot + cumIso",        # Identification and Isolation of presymptomatic
    "preSympt -> (1-preDetectSuccess)*symptomaticPropt*preSympPeriod*preSympt -> Sympt",        # Development of symptomatic disease among exposed: stationary
    "preSympt -> (1-preDetectSuccess)*(1-symptomaticPropt)*preSympPeriod*preSympt -> aSympt",     # Development of asymptomatic disease among exposed: stationary
    "Sympt -> sympDetectSuccess*Sympt -> Isot + cumIso",                             # Identification and Isolation of symptomatic
    "Sympt -> (1-sympDetectSuccess)*symptomaticPeriod*Sympt -> postSympt",                             # Transition from symptomatic disease to recovered 
    "postSympt -> postDetectSuccess*postSympt -> Isot + cumIso",                # Identification and Isolation of postsymptomatic
    "postSympt -> (1-postDetectSuccess)*postSympPeriod*postSympt -> Rt",                # Transition from post symptomatic to recovered
    "aSympt -> asympDetectSuccess*aSympt -> Isot + cumIso",                            # Identification and Isolation of asymptomatic
    "aSympt -> (1-asympDetectSuccess)*aSympPeriod*aSympt -> Rt",                            # Transition from asymptomatic disease to recovered: 
    "Isot -> isoPeriod*Isot -> Rt",                                # Transition from isolated to recovered
    "Rt -> (1-reSuscepRate)*tempImmPeriod*Rt -> Imt"                # Development of immunity: 
  )
  
  return(transitions)
}