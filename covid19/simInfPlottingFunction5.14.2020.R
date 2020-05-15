# simInfPlottingFunction5.8.2020.R

# changes from simInfPlottingFunction5.6.2020.R
# added functionality to plot sums of compartments

# changes from simInfPlottingFunction5.5.2020.R
# added option to plot cumulative infections and daily infections

# changes from simInfPlottingFunction.R
# instead of making a data table of the full trajectory, extract just the compartments/variables that are being plotted

library(ggplot2)
library(data.table)
library(zoo)

simInfPlottingFunction <- function(
                                   result,                    # trajectory of the result
                                   table = "U",               # which table to get data from
                                   compts = "I",              # compartments/variables that will be plotted
                                   groups = NULL,             # List of groups to aggregate and plot
                                   sumGroups = TRUE,          # Automatically sums over all groups
                                   uNames = names(u0),        # list of compartments
                                   vNames = NULL,             # list of variables
                                   rollM = 1,                 # number of days for rolling average
                                   confIntv = .95,            # confidence interval for plotting spread
                                   nTM = nodeTrialMat,        # node-Trial matrix
                                   tS = tspan,                # length of simulation
                                   enn = N,                   # number of nodes per trial
                                   nT = numTrials,            # number of trials in simulation
                                   startDate = startofSimDay, # start date of simulation
                                   dateBreaks = "1 month",    # plot parameter: Date axis format
                                   titleString = NULL,        # plot parameter: Title of plot
                                   xString = NULL,            # plot parameter: Title of x axis
                                   yString = NULL,            # plot parameter: Title of y axis
                                   lString = "Compartment"    # plot parameter: Title of legend
                                   ) {
  
  # if table is U or V
  if(table == "U") {
    if(compts == "newI") {
      uIndex <- which(uNames == "cumI")
      uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),uIndex,FUN = "+"))
      names(uIndices) <- "cumI"
      dt <- data.table(apply(uIndices,2,function(x) result@U[x,]))
      dt[,trial:=rep.int(nTM$trial,max(tS))][,nodeGroup:=rep.int(nTM$nodeGroup,max(tS))][,time:=sort(rep.int(tS,enn*nT),method="quick")]
      if(!is.null(groups)) {dt <- dt[nodeGroup %in% groups]}
      if(sumGroups) {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
        traj[,newI:=c(0,diff(cumI)), by = trial]
      } else {
        traj <- dt[,lapply(.SD,sum),by=.(time,trial,nodeGroup)]
        traj[,newI:=c(0,diff(cumI)), by = .(trial,nodeGroup)]
      }
      traj[,cumI:=NULL]
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
  
  # central tendency and low/high
  if(sumGroups) {
    trajCT <- traj[,lapply(.SD,median),by=.(time)]
    medTraj <- trajCT[,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
    trajSpread <- split(traj[,lapply(.SD,quantile,probs = c((1-confIntv),confIntv)),by=.(time)],1:2)
    lowTraj <- trajSpread[[1]][,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
    highTraj <- trajSpread[[2]][,lapply(.SD,rollmean,k=rollM),.SDcols=compts]
    
    medTrajMelt <- melt.data.table(medTraj,measure.vars=compts)
    names(medTrajMelt) <- c("medVar","medVal")
    lowTrajMelt <- melt.data.table(lowTraj,measure.vars=compts)
    names(lowTrajMelt) <- c("lowVar","lowVal")
    highTrajMelt <- melt.data.table(highTraj,measure.vars=compts)
    names(highTrajMelt) <- c("highVar","highVal")
    
  } else {
    trajCT <- traj[,lapply(.SD,median),by=.(time,nodeGroup)]
    medTraj <- trajCT[,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
    trajSpread <- split(traj[,lapply(.SD,quantile,probs = c((1-confIntv),confIntv)),by=.(time,nodeGroup)],1:2)
    lowTraj <- trajSpread[[1]][,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
    highTraj <- trajSpread[[2]][,lapply(.SD,rollmean,k=rollM),by=.(nodeGroup),.SDcols=compts]
    
    medTrajMelt <- melt.data.table(medTraj,measure.vars=compts)
    names(medTrajMelt) <- c("nodeGroup","medVar","medVal")
    medTrajMelt$nodeGroup <- as.factor(medTrajMelt$nodeGroup)
    lowTrajMelt <- melt.data.table(lowTraj,measure.vars=compts)
    names(lowTrajMelt) <- c("nodeGroup","lowVar","lowVal")
    highTrajMelt <- melt.data.table(highTraj,measure.vars=compts)
    names(highTrajMelt) <- c("nodeGroup","highVar","highVal")
  }
  
  # creating date list
  dateList <- rep(seq(as.Date(startDate+rollM-1,origin="2020-01-01"),by="day",length.out = max(tS) - rollM+1),length(compts))
  
  # creating data frame for ggplot
  trajDF <- data.frame(cbind(medTrajMelt,lowTrajMelt,highTrajMelt,dateList))
  
  if(sumGroups) {
    p <- ggplot(data = trajDF, aes(x=dateList,y=medVal,color=medVar))
  } else {
    p <- ggplot(data = trajDF, aes(x=dateList,y=medVal,color=nodeGroup))
  }
  
  p <- p +
    geom_line(size=1.125) +
    geom_ribbon(aes(ymin=trajDF$lowVal, ymax=trajDF$highVal),linetype=2, alpha = 0.1, size = 1.03125) +
    labs(title = titleString, x = xString, y = yString, color=lString)+
    scale_x_date(date_labels="%b %d", date_breaks=dateBreaks)+
    theme_bw()+
    theme(title=element_text(size=20),
          axis.text=element_text(size=14),
          axis.text.x=element_text(angle=90),
          legend.title = element_text(size=14),
          legend.text= element_text(size=14))
  
  return(p)
}