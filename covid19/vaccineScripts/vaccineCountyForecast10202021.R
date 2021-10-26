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


# transformation function
f <- function (x, a, b) log((x - a) / (b - x))
# inverse transformation function
finv <- function (x, a, b) (b * exp(x) + a) / (exp(x) + 1)

### read in shapefiles
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
shp <- st_read("tl_2020_41_bg.shp")

shp <- shp[shp$COUNTYFP %in% c("003","041","043"),]
shp <- shp[shp$GEOID != "410419901000",]

# getting population files
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
bgs <- data.table(read.table("lblPopData2020.txt",sep=",",header=TRUE))

# filter to county
bgs <- bgs[bgs$COUNTY==3,]

bgs <- merge(shp,bgs,by="GEOID")

# reading in vaccine data
# reading in vaccine data
setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton")
vax <- readxl::read_xlsx("alertDataBenton.xlsx")
vax <- vax %>% select(CLIENT_ID,BIRTHDATE,ETHNICITY,WHITE,BLACK,INDIAN,ASIAN,HAWAIIAN,OTHER,LATITUDE,LONGITUDE,VACCINATION_DATE,MFR,DOSE)

vaxInc <- vax %>% filter(MFR %in% c("Pfizer","Moderna") & DOSE == 1)
vaxComp <- vax %>% filter(MFR == "Janssen" | DOSE >= 2)

vaxWide <- data.table(merge(vaxInc,vaxComp, by=names(vaxInc)[c(1:11)],all=TRUE))
names(vaxWide)[c(12:17)] <- c("DATE1","MFR1","DOSE1","DATEC","MFRC","DOSEC")

vaxWide[is.na(vaxWide$MFR1),"MFR1"] <- vaxWide[is.na(vaxWide$MFR1),"MFRC"]
vaxWide[is.na(vaxWide$DOSEC),"DOSEC"] <- 0

vax <- vaxWide %>% select(CLIENT_ID,BIRTHDATE,DATE1,MFR1,DATEC,DOSEC,ETHNICITY,WHITE,BLACK,INDIAN,ASIAN,HAWAIIAN,OTHER,LATITUDE,LONGITUDE)

#### County wide analysis ####
#### County wide analysis ####
#### County wide analysis ####
#### County wide analysis ####

#### Race Ethnicity ####
#### Race Ethnicity ####
bgsREpop <- data.table(bgs[,c("COUNTYFP","Total","Latinx","White","Black","AIAN","Asian","NHPI")])[,-9][,lapply(.SD,sum),by=.(COUNTYFP)]

reFunction <- function(raceEthnicity, cmvrLista, vaxdt = vax, popdt = bgsREpop) {
  if(raceEthnicity == "Latinx") {
    vaxdt <- vaxdt[vaxdt$ETHNICITY %in% c("2135-2","Hispanic"),]
    pop <- popdt$Latinx
  } else if(raceEthnicity == "White") {
    vaxdt <- vaxdt[vaxdt$WHITE == "Y",]
    pop <- popdt$White
  } else if(raceEthnicity == "Black") {
    vaxdt <- vaxdt[vaxdt$BLACK == "Y",]
    pop <- popdt$Black
  } else if(raceEthnicity == "AIAN") {
    vaxdt <- vaxdt[vaxdt$INDIAN == "Y",]
    pop <- popdt$AIAN
  } else if(raceEthnicity == "Asian") {
    vaxdt <- vaxdt[vaxdt$ASIAN == "Y",]
    pop <- popdt$Asian
  } else if(raceEthnicity == "NHPI") {
    vaxdt <- vaxdt[vaxdt$HAWAIIAN == "Y",]
    pop <- popdt$NHPI
  } else {
    pop <- popdt$Total
  }
  
  # daily vaccination totals
  
  vaxBD <- data.table(dplyr::count(as_tibble(vaxdt), DATE1))
  vaxBD <- vaxBD[-which(is.na(vaxBD$DATE1)),]
  
  vaxBD$cumulative <- cumsum(vaxBD$n)
  vaxBD$percent <- vaxBD$cumulative/pop
  
  vaxBD$Day <- as.Date(vaxBD$DATE1)-min(as.Date(vaxBD$DATE1))
  
  # counterfactual maximum vaccination rate
  cmvrL <- cmvrLista[1]
  cmvrU <- cmvrLista[2]
  
  # lower estimate
  
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrL)
  
  train <- subset(vaxBD,DATE1 < "2021-09-01")
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = train)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrL)
  
  # Create 95% confidence intervals
  lwr <- finv(p$fit - (1.96 * p$se.fit),0,cmvrL)
  
  # create data table with estimates and actuals
  predDT <- data.table(data.frame(Day = seq(0:500),
                                  DATE1 = seq(min(as.Date(vaxBD$DATE1),na.rm=TRUE), min(as.Date(vaxBD$DATE1),na.rm=TRUE)+500, by="days"),
                                  prediction = vaxPrediction,lwr=lwr))
  predDT <- merge(predDT,vaxBD[,c("Day","percent")],by="Day",all.x=TRUE)
  
  # upper estimate
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrU)
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = vaxBD)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrU)
  
  # Create 95% confidence intervals
  upr <- finv(p$fit + (1.96 * p$se.fit),0,cmvrU)
  predDT$upr <- upr
  
  # plot results
  p <- ggplot(predDT, aes(x=DATE1, y=prediction)) +
    geom_line(color="red")+
    geom_point(aes(x=DATE1, y=percent), size=1) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=DATE1, fill="band"), alpha=0.3, fill="deepskyblue") +
    scale_y_continuous(label=scales::percent, limits=c(0,1)) +
    labs(y = "Vaccinated",
         title = paste0("Forecast of Benton County Vaccination rate - ",raceEthnicity,"; all ages"),
         subtitle = "Generalized Additive Model") +
    theme_bw()+
    theme(axis.title.x = element_blank())
  
  setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/forecasts")
  ggsave(paste0("BentonCountyForecast",raceEthnicity,Sys.Date(),".png"), plot = p, width = 9, height = 6, units = "in")
  
  return(p)
}

cmvrList <- list(c(.670,.75),c(.42,.5),c(.5,.6),c(.5,.6),c(.45,.55),c(.65,.75),c(.62,.7))
reList <- list("Total","Latinx","AIAN","Asian","Black","NHPI","White")

vaxRaceEthnicity <- mapply(reFunction,reList, cmvrList)



#### Age groups ####
#### Age groups ####
vax$Age <- age_calc(as.Date(vax$BIRTHDATE), enddate = Sys.Date(), units = "years", precise = FALSE)

ageFunction <- function(ageGroup, ageGroupLabel, cmvrLista, pop, vaxdt = vax) {
  vaxdt <- vaxdt[vaxdt$Age >= ageGroup[1],]
  vaxdt <- vaxdt[vaxdt$Age <= ageGroup[2],]

  # daily vaccination totals
  vaxBD <- data.table(dplyr::count(as_tibble(vaxdt), DATE1))
  vaxBD <- vaxBD[-which(is.na(vaxBD$DATE1)),]
  
  vaxBD$cumulative <- cumsum(vaxBD$n)
  vaxBD$percent <- vaxBD$cumulative/pop
  
  vaxBD$Day <- as.Date(vaxBD$DATE1)-min(as.Date(vaxBD$DATE1))
  
  # counterfactual maximum vaccination rate
  cmvrL <- cmvrLista[1]
  cmvrU <- cmvrLista[2]
  
  # lower estimate
  
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrL)

  train <- subset(vaxBD,DATE1 < "2021-09-01")
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = train)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrL)
  
  # Create 95% confidence intervals
  lwr <- finv(p$fit - (1.96 * p$se.fit),0,cmvrL)
  
  # create data table with estimates and actuals
  predDT <- data.table(data.frame(Day = seq(0:500),
                                  DATE1 = seq(min(as.Date(vaxBD$DATE1),na.rm=TRUE), min(as.Date(vaxBD$DATE1),na.rm=TRUE)+500, by="days"),
                                  prediction = vaxPrediction,lwr=lwr))
  predDT <- merge(predDT,vaxBD[,c("Day","percent")],by="Day",all.x=TRUE)
  
  # upper estimate
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrU)
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = vaxBD)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrU)
  
  # Create 95% confidence intervals
  upr <- finv(p$fit + (1.96 * p$se.fit),0,cmvrU)
  predDT$upr <- upr
  
  # plot results
  p <- ggplot(predDT, aes(x=DATE1, y=prediction)) +
    geom_line(color="red")+
    geom_point(aes(x=DATE1, y=percent), size=1) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=DATE1, fill="band"), alpha=0.3, fill="deepskyblue") +
    scale_y_continuous(label=scales::percent, limits=c(0,1)) +
    labs(y = "Vaccinated",
         title = paste0("Forecast of Benton County Vaccination rate - ",ageGroupLabel),
         subtitle = "Generalized Additive Model") +
    theme_bw()+
    theme(axis.title.x = element_blank())
  
  setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/forecasts")
  ggsave(paste0("BentonCountyForecast",ageGroupLabel,Sys.Date(),".png"), plot = p, width = 9, height = 6, units = "in")
  
  return(p)
}

ageGroupL <- list(c(12,17),c(18,29),c(45,64),c(65,74))
ageGroupLL <- list("12 to 17","18 to 29","45 to 64","65 to 74")
cmvrList <- list(c(.63,.7),c(.62,.7),c(.76,.85),c(1,1))
popList <- list(5985,27610,20096,9535)

# ageGroupL <- list(c(12,17),c(18,29),c(30,44),c(45,64),c(65,74),c(75,120))
# ageGroupLL <- list("12 to 17","18 to 29","30 to 44","45 to 64","65 to 74","75 and older")
# cmvrList <- list(c(.63,.7),c(.62,.7),c(1,1),c(.76,.85),c(1,1),c(1,1))
# popList <- list(5985,27610,11142,20096,9535,5813)

vaxAge <- mapply(ageFunction,ageGroupL,ageGroupLL,cmvrList,popList)


#### Communities ####
#### Communities ####

# geocode data
vaxCoords <- vax[-which(is.na(vax$LATITUDE)),]
vaxSF <- st_as_sf(vaxCoords,coords=c("LONGITUDE","LATITUDE"),crs="NAD83")

# mapping vaccine residence to bgs
vaxBGs <- st_join(vaxSF,bgs,join=st_within)

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton")
bgComm <- read_excel("bgCommunity.xlsx")
bgComm$GEOID <- as.character(bgComm$GEOID)
vaxBGs <- merge(vaxBGs,bgComm,by="GEOID")

bgsPop <- merge(bgs,bgComm,by="GEOID")
popList <- data.table(bgsPop[,c("community","Total")])[,-3][,lapply(.SD,sum),by=.(community)]

commFunction <- function(community, pop, cmvrLista, vaxdt = vaxBGs) {
  vaxdt <- vaxdt[vaxdt$community == community,]
  
  # daily vaccination totals
  
  vaxBD <- data.table(dplyr::count(as_tibble(vaxdt), DATE1))
  vaxBD <- vaxBD[-which(is.na(vaxBD$DATE1)),]
  
  vaxBD$cumulative <- cumsum(vaxBD$n)
  vaxBD$percent <- vaxBD$cumulative/pop
  
  ###
  print(max(vaxBD$percent))
  ###
  
  vaxBD$Day <- as.Date(vaxBD$DATE1)-min(as.Date(vaxBD$DATE1))
  
  # counterfactual maximum vaccination rate
  cmvrL <- cmvrLista[1]
  cmvrU <- cmvrLista[2]
  
  # lower estimate
  
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrL)
  
  train <- subset(vaxBD,DATE1 < "2021-09-01")
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = train)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrL)
  
  # Create 95% confidence intervals
  lwr <- finv(p$fit - (1.96 * p$se.fit),0,cmvrL)
  
  # create data table with estimates and actuals
  predDT <- data.table(data.frame(Day = seq(0:500),
                                  DATE1 = seq(min(as.Date(vaxBD$DATE1),na.rm=TRUE), min(as.Date(vaxBD$DATE1),na.rm=TRUE)+500, by="days"),
                                  prediction = vaxPrediction,lwr=lwr))
  predDT <- merge(predDT,vaxBD[,c("Day","percent")],by="Day",all.x=TRUE)
  
  # upper estimate
  vaxBD$percentT <- f(vaxBD$percent,0,cmvrU)
  
  gam_m <- gam(percentT ~ s(as.numeric(Day)),method = "REML",data = vaxBD)
  
  p <- predict(gam_m, data.frame(Day = seq(0:500)), type = "link", se.fit = TRUE)
  
  vaxPrediction <- finv(p$fit,0,cmvrU)
  
  # Create 95% confidence intervals
  upr <- finv(p$fit + (1.96 * p$se.fit),0,cmvrU)
  predDT$upr <- upr
  
  # plot results
  p <- ggplot(predDT, aes(x=DATE1, y=prediction)) +
    geom_line(color="red")+
    geom_point(aes(x=DATE1, y=percent), size=1) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=DATE1, fill="band"), alpha=0.3, fill="deepskyblue") +
    scale_y_continuous(label=scales::percent, limits=c(0,1)) +
    labs(y = "Vaccinated",
         title = paste0("Forecast of Benton County Vaccination rate - ",community,"; all ages"),
         subtitle = "Generalized Additive Model") +
    theme_bw()+
    theme(axis.title.x = element_blank())
  
  setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/forecasts")
  ggsave(paste0("BentonCountyForecast",community,Sys.Date(),".png"), plot = p, width = 9, height = 6, units = "in")

  return(p)
}

cmvrList <- list(c(.67,.78),c(.62,.7),c(.62,.7),c(.65,.7),c(.51,.6),c(.57,.65),c(.63,.7))

vaxCommunity <- mapply(commFunction,popList$community,popList$Total,cmvrList)


# https://github.com/braydengerrard/COVIDVaccineForecast/blob/main/Covidvaccines.R