# simInfPlottingFunction.R

library(ggplot2)
library(data.table)
library(zoo)

simInfPlottingFunction <- function(
                                   dt,                        # trajectory of the result
                                   compts,                    # compartments that will be plotted
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
  
  # function to compute SSE
  sseFunction <- function(trial, tendList=trajSpread) {
    lapply(tendList, function(x) norm(matrix(x$I-trial$I),type="f"))
  }
  
  # Sums all nodes for each trial
  traj <- merge(x=dt,y=nTM,by="node")
  traj <- traj[,lapply(.SD,sum),by=.(time,trial)]
  
  if(compts %in% c("prevalence","phi")) {
    traj[,prevalence:=prevalence/enn]
    traj[,phi:=phi/enn]
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
