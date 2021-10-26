# https://github.com/braydengerrard/COVIDVaccineForecast/blob/main/Covidvaccines.R

# libraries
library(sf)
library(tidyverse)
library(tidycensus)
library(ggplot2)
library(data.table)
library(ggmap)
library(RColorBrewer)
library(lubridate)
library(eeptools)
library(kableExtra)
library(dplyr)
library(scales)
library(ggplot2)
library(survival)
library(mgcv)

# setting date of birth cutoff
dobCutoff <- "2009-10-22"

# reading in vaccine data
setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton")
vax <- readxl::read_xlsx("alertDataBenton.xlsx",
                         col_types = c("numeric","text","text","text","date",rep("skip",26),"date",rep("skip",6),"text","text",rep("skip",5)))

# # temporary dataset for development
# vsave <- vax
# vax <- vax[c(1:100),]

# separate incomplete from complete and recombine into wider format
vaxInc <- vax %>% filter(MFR %in% c("Pfizer","Moderna") & DOSE == 1)
vaxComp <- vax %>% filter(MFR == "Janssen" | DOSE >= 2)

names(vaxInc) <- c("CLIENT_ID","LAST_NAME","FIRST_NAME","MIDDLE_NAME","BIRTHDATE","FIRSTDOSE_DATE","MFR","FIRSTDOSE")
names(vaxComp) <- c("CLIENT_ID","LN","FN","MN","BD","COMPLETEDOSE_DATE","Manu","COMPLETEDOSE")

vaxWide <- data.table(merge(vaxInc,vaxComp, by="CLIENT_ID",all=TRUE))
vaxWide[is.na(vaxWide$LAST_NAME),c("LAST_NAME","FIRST_NAME","MIDDLE_NAME","BIRTHDATE","MFR")] <- vaxWide[is.na(vaxWide$LAST_NAME),c("LN","FN","MN","BD","Manu")]
vaxWide[is.na(vaxWide$COMPLETEDOSE),"COMPLETEDOSE"] <- 0

vax <- vaxWide %>% select(c(CLIENT_ID,LAST_NAME,FIRST_NAME,MIDDLE_NAME,BIRTHDATE,FIRSTDOSE_DATE,MFR,COMPLETEDOSE_DATE,COMPLETEDOSE))

# paste name and DOB
vax$nameDOB <- paste0(vax$LAST_NAME,vax$BIRTHDATE)

# read in student data
student <- readxl::read_excel("csdRE10042021.xlsx",
                              col_types = c("text","text",rep("skip",5),"date","text","text","skip",rep("text",8),rep("skip",6),rep("text",4)))
student$fullAddress <- paste0(student$`Physical Address`,", ",student$`Physical City`,", OR ",student$`Physical Zip`)

# restrict to only Benton County
countyResidence <- data.table(table(student$`Resident County`))
names(countyResidence) <- c("County of residence","Number of students")
student <- student[student$`Resident County` == "Benton",]

# restrict to students 12 years and older
student <- student[student$DOB < dobCutoff,]

student$`Legal Last Name` <- toupper(student$`Legal Last Name`)
student$nameDOB <- paste0(student$`Legal Last Name`,student$DOB)

# flag vaccinated students
student$vax <- 0
student[student$nameDOB %in% vax$nameDOB,"vax"] <- 1

# merge vaccination data
student <- merge(student,vax,by = "nameDOB",all.x=TRUE)

# transformation function
f <- function (x, a, b) log((x - a) / (b - x))
# inverse transformation function
finv <- function (x, a, b) (b * exp(x) + a) / (exp(x) + 1)

#############
#############
#############

vaxFunction <- function(school, schoolLabel, raceEthnicity, cmvrLista, vaxdt = student) {
  
  if(cmvrLista[1]<1) {
  
    vaxdt <- subset(vaxdt, `School Name` %in% school)
    
    if(raceEthnicity == "Latinx") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "Hispanic")
    } else if(raceEthnicity == "White") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "White")
    } else if(raceEthnicity == "Black") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "Black/African American")
    } else if(raceEthnicity == "AIAN") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "American Indian/Alaska Native")
    } else if(raceEthnicity == "Asian") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "Asian")
    } else if(raceEthnicity == "NHPI") {
      vaxdt <- subset(vaxdt,`FAFSA Race` == "Native Hawaiian/Pacific Islander")
    }
    
    vaxBD <- data.table(dplyr::count(as_tibble(vaxdt[vaxdt$vax==1,]), FIRSTDOSE_DATE))
    vaxBD$cumulative <- cumsum(vaxBD$n)
    vaxBD$percent <- vaxBD$cumulative/nrow(vaxdt)
    if(length(which(is.na(vaxBD$FIRSTDOSE_DATE))) > 0) {vaxBD <- vaxBD[-which(is.na(vaxBD$FIRSTDOSE_DATE)),]}
    
    intervalCutoff <- vaxBD$FIRSTDOSE_DATE[which(vaxBD$percent==max(vaxBD$percent))]
    
    vaxBD$Day <- as.Date(vaxBD$FIRSTDOSE_DATE)-as.Date(vaxBD$FIRSTDOSE_DATE[1])
    
    # counterfactual maximum vaccination rate
    cmvrL <- cmvrLista[1]
    cmvrU <- cmvrLista[2]
    
    # lower estimate
    
    vaxBD$percentT <- f(vaxBD$percent,0,cmvrL)
    
    train <- subset(vaxBD,FIRSTDOSE_DATE < "2021-09-01")
    
    gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = train)
    
    p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
    
    vaxPrediction <- finv(p$fit,0,cmvrL)
    
    # Create 95% confidence intervals
    lwr <- finv(p$fit - (1.96 * p$se.fit),0,cmvrL)
    
    # create data table with estimates and actuals
    predDT <- data.table(data.frame(Day = seq(0:500),
                                    FIRSTDOSE_DATE = seq(min(as.Date(vaxBD$FIRSTDOSE_DATE),na.rm=TRUE), min(as.Date(vaxBD$FIRSTDOSE_DATE),na.rm=TRUE)+500, by="days"),
                                    prediction = vaxPrediction,lwr=lwr))
    predDT <- merge(predDT,vaxBD[,c("Day","percent")],by="Day",all.x=TRUE)
    predDT[predDT$lwr < max(vaxBD$percent) & (predDT$FIRSTDOSE_DATE > intervalCutoff),"lwr"] <- max(vaxBD$percent)
    predDT[predDT$FIRSTDOSE_DATE <= intervalCutoff,"lwr"] <- predDT$prediction[predDT$FIRSTDOSE_DATE <= intervalCutoff]
    
    # upper estimate
    vaxBD$percentT <- f(vaxBD$percent,0,cmvrU)
    
    gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = vaxBD)
    
    p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
    
    vaxPrediction <- finv(p$fit,0,cmvrU)
    
    # Create 95% confidence intervals
    upr <- finv(p$fit + (1.96 * p$se.fit),0,cmvrU)
    predDT$upr <- upr
    predDT[predDT$FIRSTDOSE_DATE <= intervalCutoff,"upr"] <- predDT$prediction[predDT$FIRSTDOSE_DATE <= intervalCutoff]
    
    # plot results
    p <- ggplot(predDT, aes(x=FIRSTDOSE_DATE, y=prediction)) +
      geom_line(color="red")+
      geom_point(aes(x=FIRSTDOSE_DATE, y=percent), size=1) +
      geom_ribbon(aes(ymin=lwr, ymax=upr, x=FIRSTDOSE_DATE, fill="band"), alpha=0.3, fill="deepskyblue") +
      scale_y_continuous(label=scales::percent, limits=c(0,1)) +
      labs(y = "Vaccinated",
           title = paste0("Forecast of CSD Vaccination rate - age 12+; ",raceEthnicity),
           subtitle = ifelse(length(school)<4,paste(school,collapse=", "),schoolLabel)) +
      theme_bw()+
      theme(axis.title.x = element_blank())
    
    setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/csd")
    ggsave(paste0("CSDForecast",raceEthnicity,schoolLabel,Sys.Date(),".png"), plot = p, width = 9, height = 6, units = "in")

    return(p)
  } else {return(NULL)}
}

schoolList <- list(unique(student$`School Name`)[c(1,2,5)],
                   unique(student$`School Name`)[c(3,4,6)],
                   unique(student$`School Name`)[c(1:7)],
                   unique(student$`School Name`)[1],
                   unique(student$`School Name`)[2],
                   unique(student$`School Name`)[3],
                   unique(student$`School Name`)[4],
                   unique(student$`School Name`)[5],
                   unique(student$`School Name`)[6])
schoolLL <- as.list(rep(c("HS","MS","District","CVHS","CHS","CMS","FS","CP","LP"),9))
cmvrList <- list(c(.71,.8),
                 c(.69,.5),
                 c(.7,.8),
                 c(.74,.82),
                 c(.71,.8),
                 c(.66,.75),
                 c(.78,.88),
                 c(.78,.88),
                 c(.7,.8),
                 c(.52,.62),
                 c(.43,.55),
                 c(.46,.57),
                 c(.55,.65),
                 c(.7,.8),
                 c(1,1), # observation set too small for fitting model
                 c(1,1), # observation set too small for fitting model
                 c(1,1), # observation set too small for fitting model
                 c(.44,.54),
                 c(.76,.85),
                 c(.73,.8),
                 c(.74,.82),
                 c(.76,.85),
                 c(.78,.85),
                 c(.69,.8),
                 c(.8,.9),
                 c(.6,.8),
                 c(.8,.9)
                 )
reList <- as.list(rep(c("Total","Latinx","White"),each=9)) 

vaxForecasts <- mapply(vaxFunction,schoolList,schoolLL,reList,cmvrList)

# https://github.com/braydengerrard/COVIDVaccineForecast/blob/main/Covidvaccines.R