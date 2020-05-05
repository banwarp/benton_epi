# studentEventFunction.R

studentEventFunction <- function(
                                  compart,                   # compartment index
                                  fracTable = studentPropTable,
                                  nodeList = studentNodes,
                                  pop = studentPop,
                                  nT = numTrials,
                                  enn = N,
                                  sRD = studentReturnDate,
                                  sRS = studentReturnSpread
                                ) {
  frac <- fracTable[compart, "frac"]
  numStud <- round(pop*frac,0)
  nodeTable <- data.frame(table(sample(nodeList,numStud,replace=TRUE)))
  nodeTable$Var1 <- as.numeric(as.character(nodeTable$Var1))
  
  if(nrow(nodeTable) > 0) {
    trialMat <- rep(0:(nT-1),each=nrow(nodeTable))*enn
    nodeAssignments <- data.frame(node=trialMat+nodeTable$Var1,number=rep(nodeTable$Freq,nT))
    
    eventDF <- data.frame(
      event = rep("enter",nrow(nodeAssignments)),
      time = rep(sRD+sample(-sRS:sRS,nrow(nodeTable),replace=TRUE),nT),
      node = nodeAssignments$node,
      dest = 0,
      n = nodeAssignments$number,
      proportion = 0,
      select = compart,
      shift = 0
    )
    
    return(eventDF)
  } else {return(NULL)}
  
 
}
