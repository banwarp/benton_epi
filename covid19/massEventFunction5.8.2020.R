# massEventFunction5.8.2020.R

# changes from studentEventFunction.R
# added function for super spread event

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

superFunction <- function(
                          ind,
                          IList=parmList$superInfections,
                          nodeList=sort(sample(1:N,min(N,parmList$superNodes))),
                          dayList=parmList$superDate,
                          spreadList=parmList$superSpread,
                          startDay=startofSimDay,
                          nT = numTrials,
                          enn = N
                          ) {
  I <- IList[[ind]]
  nodes <- nodeList[[ind]]
  day <- dayList[[ind]]
  spread <- spreadList[[ind]]
  if(is.character(day)) {
    day <- as.numeric(as.Date(day)) - as.numeric(as.Date("2020-01-01")) - startDay
  }
  nodeTable <- data.frame(table(sample(nodes,I,replace=TRUE)))
  nodeTable$Var1 <- as.numeric(as.character(nodeTable$Var1))
  if(nrow(nodeTable) > 0) {
    trialMat <- rep(0:(nT-1),each=nrow(nodeTable))*enn
    nodeAssignments <- data.frame(node=trialMat+nodeTable$Var1,number=rep(nodeTable$Freq,nT))  
    
    eventDF <- data.frame(
      event = rep("intTrans",nrow(nodeAssignments)),
      time = rep(day+sample(-spread:spread,nrow(nodeTable),replace=TRUE),nT),
      node = nodeAssignments$node,
      dest = 0,
      n = nodeAssignments$number,
      proportion = 0,
      select = 1,
      shift = 1
    )
    
    return(eventDF)
  } else {return(NULL)}
}