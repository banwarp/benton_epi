# eventFunctionsSchool_06.20.2020.R

# changes from eventFunctions5.16.2020.R
# Changed the code to have stationary transfers
# Changed the subroutines to handle "cross" transfers
# Added code for transitory transfers

# changes from eventFunctions5.14.2020.R
# Added option to limit infection-spreader/massEntry event contact to certain nodeGroups

# changes from massEventFunction5.8.2020.R
# added transferEventFunction to this script to clean up main script

# changes from massEntryEventFunction.R
# added function for infection spread event

nodeSampleFun <- function(x,nL=nList,nTL=numTList,out=outLogic) {
  if(out %in% c("out","cross") & length(nL[[x]]) == 1) {
    return(rep(nL[[x]],nTL[[x]]))
  } else {
    return(sample(nL[[x]],nTL[[x]],replace=TRUE))
  }
}

destSampleFun <- function(x,nL=nList,nTL=numTList,ntm=nTM,out=outLogic) {
  if(out == "out") {
    destSet <- ntm[-which(ntm$node %in% nL[[x]] | ntm$trial != x),(node)]
    return(sample(destSet,nTL[[x]],replace=TRUE))
  } else if(out == "cross") {
    destSet <- ntm[-which(ntm$trial != x),(node)]
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
                              popType = i,
                              numTList = numInTransferList,
                              tS = tspan,
                              nT = numTrials,
                              nList = nodeAndGroupList,
                              minProp = iGTMinP,
                              maxProp = iGTMaxP,
                              selectCol = length(compartments)+ceiling(i/2),
                              nTM = nodeTrialMat,
                              u0S = u0$S,
                              outLogic = "in"
                            ) {
  
  numTList <- numTList[[ind]]
  nList <- nList[[ind]]
  minProp <- minProp[[ind]]
  maxProp = maxProp[[ind]]
  enn <- nrow(nTM)/length(unique(nTM$trial))
  
  maxt <- max(tS)/24
  
  if(ind == 1 | ind == 3) {
    tS <- unlist(lapply(0:(maxt-1),function(x) c(8:16)+(24*x)))
  } else {
    tS <- unlist(lapply(0:(maxt-1),function(x) c(1:7,17:24)+(24*x)))
  }
  
  if(max(unlist(numTList)) > 0 & (length(nList[[1]])>1 | outLogic != "in")) {
  
  # Transfer data frame
  # Stationary transfers  
  if(popType < 3){
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
  } else {  #### Transitory transfers
    time <- lapply(numTList, function(x) sample(tS,x,replace=TRUE))
    node <- lapply(time, function(x) rep(1:enn,length(x)))
    time <- lapply(time, function(x) sort.int(rep(x,enn),method="quick"))
    node <- lapply(1:nT, function(x) node[[x]] + enn*(x-1))
    if(outLogic == "in") {
      dest <- lapply(node, function(x) data.frame(cbind("node" = x, "nodeGroup" = nTM$nodeGroup[1:enn])))
      dest <- lapply(dest, function(x) split(x$node,x$nodeGroup))
      dest <- lapply(dest,function(x) x[which(lapply(x,function(y) length(unique(y)))>1)])
      time <- unlist(time)[which(unlist(node) %in% unlist(dest))]
      node <- unlist(dest)
      dest <- lapply(dest, function(x) lapply(x, function(y) y + c(rep(1,length(unique(y))-1),-(length(unique(y))-1))))
      dest <- unlist(dest)
    } else if(outLogic == "out") {
      allDest <- lapply(node, function(x) data.frame(cbind("node" = x, "nodeGroup" = nTM$nodeGroup[1:enn])))
      allDest <- lapply(allDest, function(x) split(x$node,x$nodeGroup))
      destSet <- lapply(allDest, function(x) lapply(x, function(y) nTM[-which(nTM$node %in% y | nTM$trial != nTM[nTM$node == min(y), (trial)]),node]))
      dest <- lapply(1:length(destSet), function(x) lapply(1:length(destSet[[x]]), function(y) 
        rep(destSet[[x]][[y]],ceiling(length(allDest[[x]][[y]])/length(destSet[[x]][[y]])))[1:length(allDest[[x]][[y]])]))
      dest <- lapply(dest, function(x) lapply(x, function(y) sample(y,length(y))))
      time <- unlist(time)
      node <- unlist(allDest)
      dest <- unlist(dest)
    } else {
      dest <- unlist(lapply(node, function(x) x + c(rep(1,enn-1),-(enn-1))))
      node <- unlist(node)
      time <- unlist(time)
    }
    
    transferEvents <- data.frame(
      event = "extTrans",
      time = time,
      node = node,
      dest = dest,
      n = 0,
      proportion = round(runif(length(time),min=minProp,max=maxProp),3),
      select = selectCol,
      shift = 0
    )
  }
  
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
  
  if(popType < 3) {  ### Only retransfer for stationary 
    reTransferEvents <- transferEvents
    reTransferEvents$node <- transferEvents$dest
    reTransferEvents$dest <- transferEvents$node
    transferEvents <- rbind(transferEvents,reTransferEvents)
  }
  
  return(transferEvents)
  } else {return(NULL)}
}

infectionFunction <- function(
                          ind,
                          IList=infectionInfections,
                          nodeList=sort(sample(1:N,min(N,infectionNodes))),
                          dayList=infectionDate,
                          spreadList=infectionSpread,
                          cohortList = infectionCohort,
                          startDay=startofSimDay,
                          nT = numTrials,
                          enn = N
                          ) {
  I <- IList[[ind]]
  nodes <- nodeList[[ind]]
  day <- dayList[[ind]]
  spread <- spreadList[[ind]]
  cohort <- cohortList[[ind]]
  if(is.character(day)) {
    day <- as.numeric(as.Date(day)) - as.numeric(as.Date("2020-01-01")) - startDay
  }
  day <- day*24+8 # infection starts at 8am
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
      select = cohort,
      shift = 1
    )
    
    cumIDF <- eventDF
    cumIDF[,"event"] <- "enter"
    cumIDF[,"select"] <- cohort+14
    cumIDF[,"shift"] <- 0
    
    eventDF <- rbind(eventDF,cumIDF)
    
    return(eventDF)
  } else {return(NULL)}
}
