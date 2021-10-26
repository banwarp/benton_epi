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
library(scales)
library(jsonlite)
library(stringi)
library(rvest)
library(rgdal)

# read in school districts
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020/schoolDistricts")
ssd <- st_read("tl_2020_41_unsd.shp")

# filter to counties
ssd <- ssd[ssd$NAME %in% c("Lincoln County School District",
                           "Monroe School District 1J",
                           "Alsea School District 7J",
                           "Corvallis School District 509J",
                           "Philomath School District 17J",
                           "Greater Albany School District 8J",
                           "Central Linn School District 552",
                           "Harrisburg School District 7J",
                           "Lebanon Community School District 9",
                           "Santiam Canyon School District 129J",
                           "Scio School District 95",
                           "Sweet Home School District 55"),]

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
student <- student[student$DOB<"2009-10-01",]

student$`Legal Last Name` <- toupper(student$`Legal Last Name`)
student$nameDOB <- paste0(student$`Legal Last Name`,student$DOB)

# flag vaccinated students
student$vax <- 0
student[student$nameDOB %in% vax$nameDOB,"vax"] <- 1

# read in geocoded addresses
setwd("L:/Health/Epidemiology/Banwarth_Epi/GIS/GIS_files/address_matching")
studentGeocode <- st_read("csdMatchedAddresses10042021.shp")
studentGeocode <- studentGeocode %>% select("Student_DB")
names(studentGeocode)[1] <- "Student DBN"

# merge to data
student <- merge(studentGeocode,student)


#################################
#################################
############ Blocks #############
#################################
#################################

# read in block shapefile
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
shp <- st_read("tl_2020_41_tabblock20.shp")

shp <- shp[shp$COUNTYFP %in% c("003","041","043"),]

# intersect
student <- st_transform(student,crs = st_crs(shp))
studentBlocks <- st_join(student,shp,join=st_within)

unvaxCount <- dplyr::count(as_tibble(studentBlocks[studentBlocks$vax==0,]), GEOID20)
unvaxBlock <- merge(shp,unvaxCount, by="GEOID20",all.x=TRUE)
unvaxBlock[is.na(unvaxBlock$n),"n"] <- 0

## joining neighboring blocks
unvaxCentroid <- st_centroid(unvaxBlock)
htree <- hclust(dist(data.frame(rownames=unvaxCentroid$GEOID20, x=as.numeric(unvaxCentroid$INTPTLON20),
                       y=as.numeric(unvaxCentroid$INTPTLAT20))), method="complete")
unvaxBlock$cluster <- cutree(htree,h=8)

unvaxBlockCluster <- group_by(unvaxBlock,cluster) %>% summarize(n = sum(n),do_union = TRUE)

## only mapping blocks with 15 or more unvaccinated
unvaxBlockCluster$n15 <- ifelse(unvaxBlockCluster$n > 14,unvaxBlockCluster$n,NA)

# mapping Benton County block clusters
crds <- c(left = -123.35, bottom = 44.5, right = -123.22843, top = 44.62)

basemap <- get_stamenmap(crds, zoom = 12, maptype = "toner-background")
lines <- get_stamenmap(crds, zoom = 12, maptype = "toner-lines")
labels <- get_stamenmap(crds, zoom = 14, maptype = "toner-labels")

m1 <- ggmap(basemap) +
  geom_sf(data = unvaxBlockCluster, aes(fill = n15, alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
  geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
  scale_fill_viridis_b(option="plasma", name="Number of unvaccinated students")+
  labs(title="Count of unvaccinated students by block cluster - Corvallis",
       subtitle =paste("Restricted to block groups with 15+ unvaccinated students;", as.Date(max(vax$FIRSTDOSE_DATE,na.rm=TRUE)),sep=" "))+
  theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
  inset_ggmap(lines) +
  inset_ggmap(labels)

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/csd")

ggsave(paste0("CorvallisUnvaxCountBlockCluster",Sys.Date(),".jpg"), plot = m1, width = 6, height = 6, units = "in")


#################################
#################################
######### Block Groups ##########
#################################
#################################

# read in block group shapefile
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
shp <- st_read("tl_2020_41_bg.shp")

shp <- shp[shp$COUNTYFP %in% c("003","041","043"),]
shp <- shp[shp$GEOID != "410419901000",]

# intersect
student <- st_transform(student,crs = st_crs(shp))
studentBGs <- st_join(student,shp,join=st_within)

unvaxCount <- dplyr::count(as_tibble(studentBGs[studentBGs$vax==0,]), GEOID)
unvaxBG <- merge(shp,unvaxCount, by="GEOID",all.x=TRUE)

unvaxBG$n15 <- ifelse(unvaxBG$n > 14,unvaxBG$n,NA)

# mapping Corvallis School District block groups

crdsB <- c(left = -123.4748, bottom = 44.38875, right = -123.1623, top = 44.72024)
basemapB <- get_stamenmap(crdsB, zoom = 11, maptype = "toner-background")
linesB <- get_stamenmap(crdsB, zoom = 11, maptype = "toner-lines")
labelsB <- get_stamenmap(crdsB, zoom = 11, maptype = "toner-labels")

m2 <- ggmap(basemapB) +
  geom_sf(data = unvaxBG, aes(fill = n15, alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
  geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
  scale_fill_viridis_b(option="plasma", name="Number of unvaccinated students")+
  labs(title="Count of unvaccinated students by block group - Corvallis School District",
       subtitle =paste("Restricted to block groups with 15+ unvaccinated students;", as.Date(max(vax$FIRSTDOSE_DATE,na.rm=TRUE)),sep=" "))+
  theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
  inset_ggmap(linesB) +
  inset_ggmap(labelsB)

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/csd")

ggsave(paste0("CSDUnvaxCountBlockGroup",Sys.Date(),".jpg"), plot = m2, width = 6, height = 6, units = "in")
