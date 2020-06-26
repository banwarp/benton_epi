# eventFunctionsAgeGroups_06.19.2020.R

# changes from eventFunctions5.16.2020.R
# Incorporated age groups into events

# changes from eventFunctions5.14.2020.R
# Added option to limit super-spreader/massEntry event contact to certain nodeGroups

# changes from massEventFunction5.8.2020.R
# added transferEventFunction to this script to clean up main script

# changes from massEntryEventFunction.R
# added function for super spread event


nodeSampleFun <- function(x,nL=nList,nTL=numTList,out=outLogic) {
  if(out & length(nL[[x]]) == 1) {
      return(rep(nL[[x]],nTL[[x]]))
  } else { return(sample(nL[[x]],nTL[[x]],replace=TRUE)) }
}

destSampleFun <- function(x,nL=nList,nTL=numTList,ntm=nTM,out=outLogic) {
  if(out) {
    destSet <- ntm[-which(ntm$node %in% nL[[x]] | ntm$trial != x),(node)]
    return(sample(destSet,nTL[[x]],replace=TRUE))
  } else { return(sample(nL[[x]],nTL[[x]],replace=TRUE)) }
}

# function to change node=dest transfers to different dest
destChangeFun <- function(n=node,t=trial,g=nodeGroup,ntm){
  trialNodes <- ntm[which(ntm$trial==t & ntm$nodeGroup==g),(node)]
  if(length(trialNodes)>2){
    return(sample(trialNodes[-which(trialNodes==n)],1))  
  } else if(length(trialNodes)==2) {return(trialNodes[-which(trialNodes==n)])}
  else {return(trialNodes)}
}

transferFunction <- function(
                              ind,
                              numTList = numInTransferList,
                              tS = tspan,
                              nT = numTrials,
                              nList = nodeAndGroupList,
                              minProp = inGroupTransferMinProp,
                              maxProp = inGroupTransferMaxProp,
                              selectCol = ncol(E),
                              nTM = nodeTrialMat,
                              u0S = u0$S,
                              outLogic = FALSE
                            ) {
  
  numTList <- numTList[[ind]]
  nList <- nList[[ind]]
  minProp <- minProp[[ind]]
  maxProp = maxProp[[ind]]
  
  if(max(unlist(numTList))>0 & (length(nList[[1]])>1 | outLogic)) {
    # Transfer data frame
    transferEvents <- data.frame(
      event = "extTrans",
      time = unlist(lapply(numTList, function(x) sample(tS,x,replace=TRUE))),
      node = unlist(lapply(c(1:nT),nodeSampleFun,nL=nList,nTL=numTList,out=outLogic)),
      dest = unlist(lapply(c(1:nT),destSampleFun,nL=nList,nTL=numTList,ntm=nTM,out=outLogic)),
      n = 0,
      proportion = unlist(lapply(numTList,function(x) round(runif(x,min=minProp,max=maxProp),3))),
      select = selectCol,
      shift = 0
    )
  
    # changing node=dest transfers to different dest, dropping transfers if there is only 1 node per group
    sameIndices <- which(transferEvents$node==transferEvents$dest)
    transferSame <- transferEvents[sameIndices,]
    transferSame <- as.data.table(transferSame)[,rn:=as.numeric(rownames(transferSame))]
    transferSame <- merge(transferSame,nTM,by="node")
    transferSame <- transferSame[order(rn)]
    newDest <- transferSame[,destChangeFun(node,trial,nodeGroup,ntm=nTM),by=seq_len(nrow(transferSame))]
    transferEvents[sameIndices,"dest"] <- newDest[,V1]
    sameIndices <- which(transferEvents$node==transferEvents$dest)
    if(length(sameIndices)>0) {transferEvents <- transferEvents[-sameIndices,]}
    
    return(transferEvents)
  } else {return(NULL)}
}


massEntryEventFunction <- function(
                                  compart,                   # compartment index
                                  fracTable = massEntryPropTable,
                                  nodeList = massEntryNodes,
                                  pop = massEntryPop,
                                  nT = numTrials,
                                  enn = N,
                                  sRD = massEntryReturnDate,
                                  sRS = massEntryReturnSpread
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
                          IList=superInfections,
                          ageGroup = superAgeGroup,
                          nodeList=sort(sample(1:N,min(N,superNodes))),
                          dayList=superDate,
                          spreadList=superSpread,
                          startDay=startofSimDay,
                          nT = numTrials,
                          enn = N
                          ) {
  I <- IList[[ind]]
  aG <- superAgeGroup[[ind]]
  nodes <- nodeList[[ind]]
  day <- dayList[[ind]]
  spread <- spreadList[[ind]]
  if(is.character(day)) {
    day <- as.numeric(as.Date(day)) - as.numeric(as.Date("2020-01-01")) - startDay
  }
  if(length(nodes)>1) {
    nodeTable <- data.frame(table(sample(nodes,I,replace=TRUE)))
  } else {
    nodeTable <- data.frame(table(rep(nodes,I)))
  }
  
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
      select = 0+aG,
      shift = 1
    )
    
    return(eventDF)
  } else {return(NULL)}
}