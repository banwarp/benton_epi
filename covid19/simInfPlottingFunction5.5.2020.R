# simInfPlottingFunction5.5.2020.R

# changes from simInfPlottingFunction.R
# instead of making a data table of the full trajectory, extract just the compartments/variables that are being plotted

library(ggplot2)
library(data.table)
library(zoo)

simInfPlottingFunction <- function(
                                   result,                    # trajectory of the result
                                   table = "U",               # which table to get data from
                                   compts = "I",              # compartments/variables that will be plotted
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
    uIndex <- which(uNames %in% compts)
    uIndices <- data.table(outer((0:(enn*nT-1))*length(uNames),uIndex,FUN = "+"))
    names(uIndices) <- compts
    dt <- data.table(apply(uIndices,2,function(x) result@U[x,]))
    dt[,trial:=rep.int(nTM$trial,max(tS))]
    dt[,time:=sort(rep.int(tspan,enn*nT),method="quick")]
    traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
    
  } else if (table == "V") {
    vIndex <- which(vNames %in% compts)
    vIndices <- data.table(outer((0:(enn*nT-1))*length(vNames),vIndex,FUN="+"))
    names(vIndices) <- compts
    dt <- data.table(apply(vIndices,2,function(x) result@V[x,]))/enn
    dt[,trial:=rep.int(nTM$trial,max(tS))]
    dt[,time:=sort(rep.int(tspan,enn*nT),method="quick")]
    traj <- dt[,lapply(.SD,sum),by=.(time,trial)]
  }
  
  # central tendency and low/high
  trajCT <- traj[,lapply(.SD,median),by=.(time)]
  trajSpread <- split(traj[,lapply(.SD,quantile,probs = c((1-confIntv),confIntv)),by=.(time)],1:2)
 
  # getting element of trajList whose trial corresponds to the minSSETrial for low, central, and high
  lowTraj <- as.data.table(apply(trajSpread[[1]][,..compts],2,rollmean,k=rollM))
  lowTrajrollM <- melt.data.table(lowTraj,measure.vars=compts)
  names(lowTrajrollM) <- c("lowVar","lowVal")
  
  highTraj <- as.data.table(apply(trajSpread[[2]][,..compts],2,rollmean,k=rollM))
  highTrajrollM <- melt.data.table(highTraj,measure.vars=compts)
  names(highTrajrollM) <- c("highVar","highVal")
  
  medTraj <- as.data.table(apply(trajCT[,..compts],2,rollmean,k=rollM))
  medTrajrollM <- melt.data.table(medTraj,measure.vars=compts)
  names(medTrajrollM) <- c("medVar","medVal")
  
  # creating date list
  dateList <- rep(seq(as.Date(91+rollM-1,origin="2020-01-01"),by="day",length.out = max(tS) - rollM+1),length(compts))
  
  # creating data frame for ggplot
  trajDF <- data.frame(cbind(medTrajrollM,lowTrajrollM,highTrajrollM,dateList))
  
  p <- ggplot(data = trajDF, aes(x=dateList,y=medVal,color=medVar)) +
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