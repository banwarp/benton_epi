# simInfPlottingFunction_03.05.2023.R

# changes from simInfPlottingFunction_02.16.2023.R
# added functionality to plot central tendency of each kmeans cluster

# changes from simInfPlottingFunction_06.20.2020.R
# added in functionality for kmeans clusterings of trajectories
# added in functionality for finding a central tendency trajectory

# changes from simInfPlottingFunction5.30.2020.R
# Option for user to choose title for data, e.g. "Infections" instead of "I"
# Changed newI logic to accommodate multiple cumI compartments

# changes from simInfPlottingFunction5.19.2020.R
# Added option to plot all trajectories and highlight some random trajectories

# changes from simInfPlottingFunction5.14.2020.R
# Fixed confidence interval calculation

# changes from simInfPlottingFunction5.6.2020.R
# added functionality to plot sums of compartments

# changes from simInfPlottingFunction5.5.2020.R
# added option to plot cumulative infections and daily infections

# changes from simInfPlottingFunction.R
# instead of making a data table of the full trajectory, extract just the compartments/variables that are being plotted

library(ggplot2)
library(data.table)
library(viridis)
library(zoo)
library(FCPS)

simInfPlottingFunction <- function(
    result,                    # trajectory of the result
    table = "U",               # which table to get data from
    compts = "I",              # compartments/variables that will be plotted
    newI = FALSE,              # Logical to plot new infections
    groups = NULL,             # List of groups to aggregate and plot
    sumGroups = TRUE,          # Automatically sums over all groups
    uNames = names(u0),        # list of compartments
    vNames = NULL,             # list of variables
    rollM = 1,                 # number of days for rolling average
    allTraj = FALSE,           # logical to plot all simulation trajectories or just median and spread
    kmeansClustering = FALSE,  # logical for kmeans clustering of trajectories
    kmeansCompts= "cumI",      # compartments for kmeans clustering and central tendency
    kmeansCenters=1,           # number of kmeans centers
    includeCentralTendency = FALSE, # whether to include a single trajectory that is the central tendency according to the cumulative norm
    includeCentralTendencyCluster = FALSE, # whether to include a central tendency for each cluster according to the cumulative norm
    separateClusters = FALSE, # whether to graph clusters separately and return as a list of plots
    plotRandomTrajs = 0,         # If allTraj = true, can plot random trajectories, select the number desired
    percentile = .5,           # percentile of simulations to plot, defaults to median
    includeMedian = FALSE,      # Whether to include median
    confIntv = .90,            # two-sided confidence interval for plotting spread
    nTM = nodeTrialMat,        # node-Trial matrix
    tS = tspan,                # length of simulation
    enn = N,                   # number of nodes per trial
    nT = numTrials,            # number of trials in simulation
    startDate = startofSimDay, # start date of simulation
    byHour=FALSE,              # time span is by hours
    dateBreaks = "1 month",    # plot parameter: Date axis format
    dataTitle = NULL,          # plot parameter: Title of data
    titleString = NULL,        # plot parameter: Title of plot
    xString = NULL,            # plot parameter: Title of x axis
    yString = NULL,            # plot parameter: Title of y axis
    lString = "Median",
    clString = "Median",       # plot parameter: Title of color legend
    slString = "",             # plot parameter: Title of size legend
    cString = ""              # plot parameter: Plot caption
) {
  
  # turns two-sided confidence interval into one-sided
  confIntv <- .5+confIntv/2
  
  plotRandomTrajs <- min(plotRandomTrajs,nT)
  
  if(is.null(dataTitle)) { dataTitle <- compts}
  
  # if groups are being plotted separately, only plot median and spread, not all trajectories
  if(!sumGroups) {allTraj <- FALSE}
  
  # plotting all trajectories requires only one compartment
  if(allTraj) {compts <- compts[1]}
  
  # if table is U or V
  if(table == "U") {
    if(newI) {
      comptList <- unlist(strsplit(compts,"_"))
      uIndex <- which(uNames %in% comptList)
      uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),uIndex,FUN = "+"))
      names(uIndices) <- comptList
      dt <- data.table(apply(uIndices,2,function(x) result@U[x,]))
      dt[,compts:=rowSums(.SD)][,(comptList):=NULL][,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
        traj[,newI:=c(0,diff(compts)), by = trial]
      } else {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
        traj[,newI:=c(0,diff(compts)), by = .(trial,nodeGroup)]
      }
      traj[,compts:=NULL]
      setcolorder(traj,c("time","trial","newI","nodeGroup"))
      compts <- "newI"
    } else if(grepl("_",compts)){
      comptList <- unlist(strsplit(compts,"_"))
      uIndex <- which(uNames %in% comptList)
      uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),uIndex,FUN = "+"))
      names(uIndices) <- comptList
      dt <- data.table(apply(uIndices,2,function(x) result@U[x,]))
      dt[,(compts):=rowSums(.SD)][,(comptList):=NULL][,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
      } else {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
      }
    } else {
      uIndex <- which(uNames %in% compts)
      uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),uIndex,FUN = "+"))
      names(uIndices) <- compts
      dt <- data.table(apply(uIndices,2,function(x) result@U[x,]))
      dt[,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
      } else {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
      }
    }
  } else if (table == "V") {
    groups <- NULL
    sumGroups <- TRUE
    vIndex <- which(vNames %in% compts)
    vIndices <- data.table(outer((0:(enn*nT-1))*length(vNames),vIndex,FUN="+"))
    names(vIndices) <- compts
    dt <- data.table(apply(vIndices,2,function(x) result@V[x,]))/enn
    dt[,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
    traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
  }
  
  # median and low/high
  if(sumGroups) {
    trajCT <- traj[,lapply(.SD,quantile,probs = percentile),by=.(time)]
    medTraj <- trajCT[,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
    medTrajMelt <- melt.data.table(medTraj,measure.vars=compts)
    names(medTrajMelt) <- c("medVar","medVal")
    medTrajMelt$medVar = dataTitle
    
    if(!allTraj){
      trajSpread <- split(traj[,lapply(.SD,quantile,probs = c((1-confIntv),confIntv)),by=.(time)],1:2)
      lowTraj <- trajSpread[[1]][,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
      highTraj <- trajSpread[[2]][,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
      lowTrajMelt <- melt.data.table(lowTraj,measure.vars=compts)
      names(lowTrajMelt) <- c("lowVar","lowVal")
      highTrajMelt <- melt.data.table(highTraj,measure.vars=compts)
      names(highTrajMelt) <- c("highVar","highVal")
    }
  } else {
    trajCT <- traj[,lapply(.SD,quantile,probs = percentile),by=.(time,nodeGroup)]
    medTraj <- trajCT[,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
    medTrajMelt <- melt.data.table(medTraj,measure.vars=compts)
    names(medTrajMelt) <- c("nodeGroup","medVar","medVal")
    medTrajMelt$nodeGroup <- as.factor(medTrajMelt$nodeGroup)
    
    if(!allTraj) {
      trajSpread <- split(traj[,lapply(.SD,quantile,probs = c((1-confIntv),confIntv)),by=.(time,nodeGroup)],1:2)
      lowTraj <- trajSpread[[1]][,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
      highTraj <- trajSpread[[2]][,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
      lowTrajMelt <- melt.data.table(lowTraj,measure.vars=compts)
      names(lowTrajMelt) <- c("nodeGroup","lowVar","lowVal")
      highTrajMelt <- melt.data.table(highTraj,measure.vars=compts)
      names(highTrajMelt) <- c("nodeGroup","highVar","highVal")
    }
  }
  
  # kmeans clustering
  if(kmeansClustering) {
    if(grepl("_",kmeansCompts)){
      k.comptList <- unlist(strsplit(kmeansCompts,"_"))
      k.uIndex <- which(uNames %in% k.comptList)
      k.uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),k.uIndex,FUN = "+"))
      names(k.uIndices) <- k.comptList
      dt <- data.table(apply(k.uIndices,2,function(x) result@U[x,]))
      dt[,(kmeansCompts):=rowSums(.SD)][,(k.comptList):=NULL][,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        k.traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
      } else {
        k.traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
      }
    } else {
      k.uIndex <- which(uNames %in% kmeansCompts)
      k.uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),k.uIndex,FUN = "+"))
      names(k.uIndices) <- kmeansCompts
      dt <- data.table(apply(k.uIndices,2,function(x) result@U[x,]))
      dt[,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        k.traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
      } else {
        k.traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
      }
    }
    k.traj <- k.traj[,c(2:3)]
    traj.mat <- matrix(data.frame(k.traj)[,2],nrow=nT)
    traj.kmeans <- kmeans(traj.mat,centers=kmeansCenters)
    traj.dist <- as.matrix(dist(traj.mat,diag=TRUE,upper=TRUE))
    centralTend <- which.min(apply(traj.dist,2,function(x) sum(x*x)))
    
  } else {
    traj.kmeans <- data.frame("cluster" = rep(1,nT))
    includeCentralTendencyCluster <- FALSE
    separateClusters <- FALSE
  }
  
    # creating date list
  if(byHour) {
    dateList <- seq(ISOdate(2020,1,1)+startDate*(24*60*60)-(11*60*60),ISOdate(2020,1,1)+startDate*(24*60*60)-(12*60*60)+max(tS)*60*60,by="hours")
  } else {
    dateList <- rep(seq(as.Date(startDate+rollM-1,origin="2020-01-01"),by="day",length.out = max(tS) - rollM+1),length(compts))
  }
  
  # creating data frame for ggplot
  if(!allTraj) {
    trajDF <- data.frame(cbind(medTrajMelt,lowTrajMelt,highTrajMelt,dateList))
  }  else {
    trajDF <- data.frame(cbind(medTrajMelt,dateList))
  }
  
  # creating a data frame with all trajectories for graphing all trajectories and/or for finding the central tendency
  allTrajDF <- data.frame(cbind(traj[,c(2:4)],rep(dateList,each=nT)),rep(traj.kmeans$cluster,length(tS)))
  names(allTrajDF) <- c("trial","compartment","nodeGroup","dateList","cluster")
  allTrajDF$trial <- factor(allTrajDF$trial)
  allTrajDF$cluster <- factor(allTrajDF$cluster)
  
  
  # plotting
  p <- ggplot()
  
  if(allTraj){
    p <- p +
      geom_line(data=allTrajDF, aes(x=dateList,y=compartment,group=trial,color= cluster), size=.1, alpha = 0.3) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    if(plotRandomTrajs > 0) {
      randomTrials <- sample(nT,plotRandomTrajs)
      randomTrajDF <- allTrajDF[allTrajDF$trial %in% randomTrials,]
      colorVec <- rep(viridis(plotRandomTrajs),each=max(tS))
      p <- p + geom_line(data=randomTrajDF, aes(x=randomTrajDF$dateList,y=randomTrajDF$compartment, group=randomTrajDF$trial), color=colorVec, size=.75)
    }
  } else {
    p <- p +
      geom_ribbon(data = trajDF,
                  aes(x = dateList, ymin=lowVal, ymax=highVal, color = medVar),linetype=2, alpha = 0.1, size = 1.03125) +
      theme_bw()
  }
  
  if(includeMedian) {
    
    if(sumGroups) {
      p <- p + geom_line(data = trajDF, aes(x=dateList,y=medVal,color=medVar), size = 1.125)
    } else {
      p <- p + geom_line(data = trajDF, aes(x=dateList,y=medVal,color=nodeGroup), size = 1.125)
    }
  }
  
  if(includeCentralTendency) {
    centralTendDF <- allTrajDF[allTrajDF$trial==centralTend,][,c("trial","compartment","dateList")]
    centralTendDF$trial <- "Central Tendency"
    p <- p + geom_line(data = centralTendDF, aes(x=dateList,y=compartment,color=trial), size = 1.125)
  }
  
  if(includeCentralTendencyCluster) {
    # tendency of each cluster
    cTfunction <- function(clus) {
      dt <- traj.dist[which(traj.kmeans$cluster==clus),which(traj.kmeans$cluster==clus)]
      clustCT <- which.min(apply(dt,2,function(x) sum(x*x)))
    }
    cTList <- as.numeric(names(unlist(lapply(sort(unique(traj.kmeans$cluster)),cTfunction))))
    cTDF <- data.frame(trial=cTList,Cluster=1:length(cTList),Relative.Likelihood =as.numeric(table(allTrajDF$cluster)*kmeansCenters/(nT*length(tS))))
    
    centralTendClust <- merge(allTrajDF[allTrajDF$trial %in% cTList,][,c("trial","compartment","dateList")],cTDF,by="trial")
    centralTendClust$Cluster <- factor(centralTendClust$Cluster)
    p <- p + geom_line(data = centralTendClust, aes(x=dateList,y=compartment,color=Cluster,size=Cluster)) +
      scale_size_manual(values = cTDF$Relative.Likelihood)
  }
  
  p <- p +
    labs(title = titleString, x = xString, y = yString, color=clString,size=slString, caption = cString)+
    theme(title=element_text(size=16),
          axis.text=element_text(size=14),
          axis.text.x=element_text(angle=90),
          legend.title = element_text(size=14),
          legend.text= element_text(size=14))
  
  if(!byHour) {
    p <- p + scale_x_date(date_labels="%b %d", date_breaks=dateBreaks)
  }
  
  if(!separateClusters) {
    return(p)  
  } else {
    p1 <- ggplot()
    pbuild <- ggplot_build(p)
    pcolors <- unique(pbuild$data[[2]][["colour"]])
    
    separateClusterFunction <- function(clust) {
      pl <- p1 +
        geom_line(data=allTrajDF[allTrajDF$cluster ==clust,], aes(x=dateList,y=compartment,group=trial),color= pcolors[clust], size=.1, alpha = 0.3) +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_line(data = centralTendClust[centralTendClust$Cluster==clust,], aes(x=dateList,y=compartment),color=pcolors[clust],size=1.25) +
        labs(title = paste0(titleString,"; Cluster ",clust), x = xString, y = yString, color=lString, caption = cString)+
        theme(title=element_text(size=16),
              axis.text=element_text(size=14),
              axis.text.x=element_text(angle=90),
              legend.title = element_text(size=14),
              legend.text= element_text(size=14))
    }
    
    pList <- lapply(c(1:kmeansCenters),separateClusterFunction)
    return(list(p,pList))
    }
  
}