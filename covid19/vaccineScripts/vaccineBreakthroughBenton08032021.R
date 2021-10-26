# vaccineBreakthroughBenton08032021.R

caseFile <- "bcases02012020_10242021.xlsx"

# libraries
library(sf)
library(tidyverse)
library(tidycensus)
library(data.table)

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

# reading in case data
setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/benton")
case <- readxl::read_xlsx(caseFile,
                          col_types = c("numeric",rep("skip",22),"date",rep("skip",5),"date","text","date"))

# # temp dataset for development
# csave <- case
# case <- case[c(1:100),]

# keep only useful data
names(case) <- c("CaseID","ReportDateLHD","TrueCaseDate","LastFirstMid","DOB")
case$LastFirstMid <- toupper(case$LastFirstMid)
case <- case %>% separate(LastFirstMid, c("Last","FirstMid"), sep=", ")
case <- case %>% separate(FirstMid,c("First","Mid"))
case$nameDOB <- paste0(case$Last,case$DOB)

# keep only vax and case records that match a case
vaxcase <- vax[which(vax$nameDOB %in% case$nameDOB),]
vaxcase <- data.table(merge(vaxcase,case,by="nameDOB"))

# flag case before vaccination; case after vaccine complete
vaxcase$prePostVax <- 0
vaxcase[which(vaxcase$FIRSTDOSE_DATE < (as.Date(vaxcase$TrueCaseDate) - 14)),"prePostVax"] <- 1
vaxcase[which(vaxcase$COMPLETEDOSE_DATE < (as.Date(vaxcase$TrueCaseDate) - 14)),"prePostVax"] <- 2

# merge back to cases
caseVax <- data.table(merge(case,vaxcase[,c("CaseID","prePostVax","COMPLETEDOSE_DATE","MFR")],all=TRUE))
caseVax[is.na(caseVax$prePostVax),"prePostVax"] <- -1
caseVax$daysAfterComplete <- caseVax$TrueCaseDate-caseVax$COMPLETEDOSE_DATE-14
caseVax[caseVax$prePostVax<2,"daysAfterComplete"] <- NA

caseVax$AgeAtInfection <- floor(as.numeric((caseVax$TrueCaseDate-caseVax$DOB)/365))

caseVax$prePostVax <- plyr::mapvalues(caseVax$prePostVax, from = c(-1,0,1,2),to=c("Unvaccinated","Unvaccinated at time of infection","1st dose only + 14 days","Complete vaccine + 14 days"))

caseVax <- caseVax[-which(duplicated(caseVax$CaseID)),]

caseVaxPBI <- caseVax %>% select(c(CaseID,TrueCaseDate,COMPLETEDOSE_DATE,prePostVax))


# plot(caseVax$COMPLETEDOSE_DATE,caseVax$daysAfterComplete)


setwd("../alertData/benton")
# 
# write.table(caseVax,file=paste0("breakthroughCases",Sys.Date(),".csv"),col.names=TRUE,sep=",",row.names=FALSE)
write.table(caseVaxPBI,file="caseVaxPBI.csv",col.names=TRUE,sep=",",row.names=FALSE)

# time series of cases with vaccine status
