# vaccine08182021function.R
# libraries
library(sf)
library(tidyverse)
library(tidycensus)
library(ggplot2)
library(data.table)
library(ggmap)
library(RColorBrewer)
library(lubridate)
library(readxl)

# 
# # reading in shapefiles
# setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
# 
# header <- read.delim("orgeo2020.pl", header=FALSE, colClasses="character", sep="|")
# part1  <- read.delim("or000012020.pl",  header=FALSE, colClasses="character", sep="|")
# part2  <- read.delim("or000022020.pl",  header=FALSE, colClasses="character", sep="|")
# part3  <- read.delim("or000032020.pl",  header=FALSE, colClasses="character", sep="|")
# 
# 
# # -----------------------------
# colnames(header) <- c("FILEID", "STUSAB", "SUMLEV", "GEOVAR", "GEOCOMP", "CHARITER", "CIFSN", "LOGRECNO", "GEOID",
#                       "GEOCODE", "REGION", "DIVISION", "STATE", "STATENS", "COUNTY", "COUNTYCC", "COUNTYNS", "COUSUB",
#                       "COUSUBCC", "COUSUBNS", "SUBMCD", "SUBMCDCC", "SUBMCDNS", "ESTATE", "ESTATECC", "ESTATENS",
#                       "CONCIT", "CONCITCC", "CONCITNS", "PLACE", "PLACECC", "PLACENS", "TRACT", "BLKGRP", "BLOCK",
#                       "AIANHH", "AIHHTLI", "AIANHHFP", "AIANHHCC", "AIANHHNS", "AITS", "AITSFP", "AITSCC", "AITSNS",
#                       "TTRACT", "TBLKGRP", "ANRC", "ANRCCC", "ANRCNS", "CBSA", "MEMI", "CSA", "METDIV", "NECTA",
#                       "NMEMI", "CNECTA", "NECTADIV", "CBSAPCI", "NECTAPCI", "UA", "UATYPE", "UR", "CD116", "CD118",
#                       "CD119", "CD120", "CD121", "SLDU18", "SLDU22", "SLDU24", "SLDU26", "SLDU28", "SLDL18", "SLDL22",
#                       "SLDL24", "SLDL26", "SLDL28", "VTD", "VTDI", "ZCTA", "SDELM", "SDSEC", "SDUNI", "PUMA", "AREALAND",
#                       "AREAWATR", "BASENAME", "NAME", "FUNCSTAT", "GCUNI", "POP100", "HU100", "INTPTLAT", "INTPTLON",
#                       "LSADC", "PARTFLAG", "UGA")
# colnames(part1) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
#                      paste0("P00", c(10001:10071, 20001:20073)))
# colnames(part2) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
#                      paste0("P00", c(30001:30071, 40001:40073)),
#                      paste0("H00", 10001:10003))
# colnames(part3) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
#                      paste0("P00", 50001:50010))
# 
# # -----------------------------
# # Merge the data
# # -----------------------------
# combine <- Reduce(function(x,y) {merge(x, y, by=c("LOGRECNO", "STUSAB", "FILEID", "CHARITER"))}, list(header[,-7], part1[,-4], part2[,-4], part3))
# 
# # -----------------------------
# # Order the data
# # -----------------------------
# combine <- combine[order(combine$LOGRECNO), c("FILEID", "STUSAB", "SUMLEV", "GEOVAR", "GEOCOMP", "CHARITER", "CIFSN", "LOGRECNO", "GEOID",
#                                               "GEOCODE", "REGION", "DIVISION", "STATE", "STATENS", "COUNTY", "COUNTYCC", "COUNTYNS", "COUSUB",
#                                               "COUSUBCC", "COUSUBNS", "SUBMCD", "SUBMCDCC", "SUBMCDNS", "ESTATE", "ESTATECC", "ESTATENS",
#                                               "CONCIT", "CONCITCC", "CONCITNS", "PLACE", "PLACECC", "PLACENS", "TRACT", "BLKGRP", "BLOCK",
#                                               "AIANHH", "AIHHTLI", "AIANHHFP", "AIANHHCC", "AIANHHNS", "AITS", "AITSFP", "AITSCC", "AITSNS",
#                                               "TTRACT", "TBLKGRP", "ANRC", "ANRCCC", "ANRCNS", "CBSA", "MEMI", "CSA", "METDIV", "NECTA",
#                                               "NMEMI", "CNECTA", "NECTADIV", "CBSAPCI", "NECTAPCI", "UA", "UATYPE", "UR", "CD116", "CD118",
#                                               "CD119", "CD120", "CD121", "SLDU18", "SLDU22", "SLDU24", "SLDU26", "SLDU28", "SLDL18", "SLDL22",
#                                               "SLDL24", "SLDL26", "SLDL28", "VTD", "VTDI", "ZCTA", "SDELM", "SDSEC", "SDUNI", "PUMA", "AREALAND",
#                                               "AREAWATR", "BASENAME", "NAME", "FUNCSTAT", "GCUNI", "POP100", "HU100", "INTPTLAT", "INTPTLON",
#                                               "LSADC", "PARTFLAG", "UGA", paste0("P00", c(10001:10071, 20001:20073)), paste0("P00", c(30001:30071, 40001:40073)),
#                                               paste0("H00", 10001:10003), paste0("P00", 50001:50010))]
# rownames(combine) <- 1:nrow(combine)
# combine[,c(98:ncol(combine))] <- apply(combine[,c(98:ncol(combine))],2,as.numeric)
# 
# lbl <- combine[combine$COUNTY %in% c("003","041","043"),]
# 
# 
# lbl <- as.data.table(lbl)[,c("GEOCODE","COUNTY","SUMLEV",paste0("P00", c(10001:10071, 20002)))]
# names(lbl)[c(1,4,ncol(lbl))] <- c("GEOID","Total","Latinx")
# lbl$White <- rowSums(lbl[,paste0("P00",c(10003,10011:10015,10027:10036,10048:10057,10064:10068,10071))])
# lbl$Black <- rowSums(lbl[,paste0("P00",c(10004,10011,10016:10019,10027:10030,10037:10042,10048:10053,10058:10061,10064:10067,10069,10071))])
# lbl$AIAN <- rowSums(lbl[,paste0("P00",c(10005,10012,10016,10020:10022,10027,10031:10033,10037:10039,10043,10045,10048:10050,10054:10056,10058:10060,10062,10064:10066,10068,10069,10071))])
# lbl$Asian <- rowSums(lbl[,paste0("P00",c(10006,10013,10017,10020,10023,10024,10028,10031,10034,10035,10037,10040,10041,10043,10044,10046,10048,10051,10052,10054,10055,10057:10059,10061,10062,10064,10065,10067:10069,10071))])
# lbl$NHPI <- rowSums(lbl[,paste0("P00",c(10007,10014,10018,10021,10023,10025,10029,10032,10034,10036,10038,10040,10042,10043,10045,10046,10049,10051,10053,10054,10056:10058,10060,10061,10062,10064,10066:10069,10071))])
# 
# 
# bgs <- lbl[lbl$SUMLEV=="150",c("GEOID","COUNTY","Total","Latinx","White","Black","AIAN","Asian","NHPI")]
# write.table(bgs,"lblPopData2020.txt",sep=",",row.names=FALSE)


### read in shapefiles
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
shp <- st_read("tl_2020_41_bg.shp")

shp <- shp[shp$COUNTYFP %in% c("003","041","043"),]
shp <- shp[shp$GEOID != "410419901000",]

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

####################################
####################################
##### Mapping Function
####################################
####################################

mapFunction <- function(pL,crds,zL,noL=FALSE,dt,rE,date,mt) {
  
  # mapping Benton County block groups
  basemap <- get_stamenmap(crds, zoom = zL[1], maptype = mt)
  lines <- get_stamenmap(crds, zoom = zL[2], maptype = "toner-lines")
  if(!noL) {labels <- get_stamenmap(crds, zoom = zL[3], maptype = "toner-labels")}
  
  m1 <- ggmap(basemap) +
    geom_sf(data = dt, aes(fill = factor(vaxCat), alpha = .1),inherit.aes=FALSE)+
    geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
    coord_sf(datum=NA) +
    scale_fill_manual(values = c(rev(brewer.pal(4,"YlOrRd")))) +
    labs(title=paste(pL,rE,date,"Vaccination Rate",sep = " "))+
    theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
    inset_ggmap(lines) 
  
  if(!noL){m1 <- m1 + inset_ggmap(labels)}
  
  m2 <- ggmap(basemap) +
    geom_sf(data = dt, aes(fill = factor(needCat), alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
    geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
    scale_fill_manual(values = c(brewer.pal(4,"YlOrRd"))) +
    labs(title=paste(pL,rE,date,"Unvaccinated Count",sep = " "))+
    theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
    inset_ggmap(lines)
  
  if(!noL){m2 <- m2 + inset_ggmap(labels)}
  
  m3 <- ggmap(basemap) +
    geom_sf(data = dt, aes(fill = highValue, alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
    geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
    scale_fill_viridis_b(option="plasma")+
    labs(title=paste(pL,rE,date,"High value for outreach",sep = " "))+
    theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
    inset_ggmap(lines)
  
  if(!noL){m3 <- m3 + inset_ggmap(labels)}
  
  m4 <- ggmap(basemap) +
    geom_sf(data = dt, aes(fill = increase, alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
    geom_sf(data=ssd,aes(),color="blue",fill=NA,inherit.aes = FALSE,size=1)+
    scale_fill_viridis_b(option="plasma")+
    labs(title=paste(pL,rE,date,"Increase from past week",sep = " "))+
    theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
    inset_ggmap(lines)
  
  if(!noL){m4 <- m4 + inset_ggmap(labels)}
  
  ggsave(paste0(pL,rE,date,"VaxRate.jpg"), plot = m1, width = 6, height = 6, units = "in")
  ggsave(paste0(pL,rE,date,"UnvaxCount.jpg"), plot = m2, width = 6, height = 6, units = "in")
  ggsave(paste0(pL,rE,date,"Value.jpg"), plot = m3, width = 6, height = 6, units = "in")
  ggsave(paste0(pL,rE,date,"IncreaseFromPastWeek.jpg"), plot = m4, width = 6, height = 6, units = "in")
  
  return(list(m1,m2,m3,m4))
}


####################################
####################################
########## School Function
####################################
####################################

schoolFunction <- function(raceEthnicity,dt,district=ssd) {
  # filter to raceEthnicity
  dt <- dt[dt$raceEthnicity == raceEthnicity,]
  
  # merge to district shapefile
  dt <- st_intersection(dt,district)
  
  districtVax <- dt %>% data.table() %>% select(c(pop,n,stillNeeds,increase,NAME))
  districtVax <- districtVax[,lapply(.SD,sum),by=.(NAME)]
  names(districtVax) <- c("School District","Population","Number vaxxed","Unmet need","Increase from previous week")
  districtVax$"Vax rate" <- districtVax$"Number vaxxed"/districtVax$Population
  districtVax$raceEthnicity <- raceEthnicity
  
  return(districtVax)
}

####################################
####################################
########## Main Function
####################################
####################################

vaxFunction <- function(raceEthnicity,wd,alertCounty,shpCounty,placeList,coordsList,zoomList,shpFile=shp, noLabels=FALSE,mapType="toner-background",bgC = FALSE) {
  
  # reading in vaccine data
  setwd(wd)
  vax <- readxl::read_xlsx(paste0("alertData",alertCounty,".xlsx"))
  
  # temporary dataset for development
  # vsave <- vax
  # vax <- vax[c(1:100),]
  
  # remove duplicates - whether complete or in progress is fine
  vax <- vax[!duplicated(vax$CLIENT_ID),]
  
  if(raceEthnicity == "Latinx") {
    vax <- vax[vax$ETHNICITY %in% c("2135-2","Hispanic"),] 
  } else if(raceEthnicity == "White") {
    vax <- vax[vax$WHITE == "Y",]
  } else if(raceEthnicity == "Black") {
    vax <- vax[vax$BLACK == "Y",]
  } else if(raceEthnicity == "AIAN") {
    vax <- vax[vax$INDIAN == "Y",]
  } else if(raceEthnicity == "Asian") {
    vax <- vax[vax$ASIAN == "Y",]
  } else if(raceEthnicity == "NHPI") {
    vax <- vax[vax$HAWAIIAN == "Y",]
  }
  
  # geocode data
  vaxCoords <- vax[,c("CLIENT_ID","VACCINATION_DATE","LATITUDE","LONGITUDE")]
  vaxCoords <- vaxCoords[-which(is.na(vaxCoords$LATITUDE)),]
  vaxSF <- st_as_sf(vaxCoords,coords=c("LONGITUDE","LATITUDE"),crs="NAD83")
  
  # getting shape files
  setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
  bgs <- data.table(read.table("lblPopData2020.txt",sep=",",header=TRUE))
  
  # filter to county
  bgs <- bgs[bgs$COUNTY==shpCounty,]
  
  bgs <- merge(shp,bgs,by="GEOID")
  
  # mapping vaccine residence to bgs
  vaxBGs <- st_join(vaxSF,bgs,join=st_within)
  
  vaxCount <- dplyr::count(as_tibble(vaxBGs), GEOID)
  vaxCountLW <- dplyr::count(as_tibble(vaxBGs[vaxBGs$VACCINATION_DATE < (max(vaxBGs$VACCINATION_DATE) - as.difftime(6, unit="days")),]), GEOID)
   
  vaxCount <- merge(vaxCount,vaxCountLW,by="GEOID")
  names(vaxCount) <- c("GEOID","n","nLW")
  
  # merging back to spatial data
  vaxBGm <- merge(bgs,vaxCount, by="GEOID",all.x=TRUE)
  vaxBGm[is.na(vaxBGm$n),"n"] <- 0
   
  # calculating percent vaccination
  if(raceEthnicity =="Latinx") {
    vaxBGm$pop <- vaxBGm$Latinx
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$Latinx
  } else if(raceEthnicity == "White") {
    vaxBGm$pop <- vaxBGm$White
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$White
  } else if(raceEthnicity == "Black") {
    vaxBGm$pop <- vaxBGm$Black
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$Black
  } else if(raceEthnicity == "AIAN") {
    vaxBGm$pop <- vaxBGm$AIAN
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$AIAN
  } else if(raceEthnicity == "Asian") {
    vaxBGm$pop <- vaxBGm$Asian
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$Asian
  } else if(raceEthnicity == "NHPI") {
    vaxBGm$pop <- vaxBGm$NHPI
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$NHPI
  } else {
    vaxBGm$pop <- vaxBGm$Total
    # vaxBGm$percent <- vaxBGm$n/vaxBGm$Total
  }
  
  vaxBGm$percent <- vaxBGm$n/vaxBGm$pop
  vaxBGm$raceEthnicity <- raceEthnicity
  
  # capping percent at 100%
  vaxBGm[which(vaxBGm$percent>1),"percent"] <- 1
  
  # categorical vax rate; 25; 50; 75; 100%
  vaxBGm$vaxCat <- .25
  vaxBGm[which(vaxBGm$percent > .25), "vaxCat"] <- .5
  vaxBGm[which(vaxBGm$percent > .5), "vaxCat"] <- .75
  vaxBGm[which(vaxBGm$percent > .75), "vaxCat"] <- 1

  vaxBGm$stillNeeds <- pmax(0,vaxBGm$pop - vaxBGm$n)
    
  vaxBGm$needCat <- quantile(vaxBGm$stillNeeds,.5)
  vaxBGm[which(vaxBGm$stillNeeds > quantile(vaxBGm$stillNeeds,.5)), "needCat"] <- quantile(vaxBGm$stillNeeds,.75)
  vaxBGm[which(vaxBGm$stillNeeds > quantile(vaxBGm$stillNeeds,.75)), "needCat"] <- quantile(vaxBGm$stillNeeds,.9)
  vaxBGm[which(vaxBGm$stillNeeds > quantile(vaxBGm$stillNeeds,.9)), "needCat"] <- max(vaxBGm$stillNeeds)
  
  vaxBGm$highValue <- (vaxBGm$percent/max(vaxBGm$percent))*vaxBGm$stillNeeds/max(vaxBGm$stillNeeds)
  
  vaxBGm$increase <- vaxBGm$n-vaxBGm$nLW
  
  # mapping
  setwd(paste0(wd,"./maps"))
  # using most recent vaccination date
  maps <- mapply(mapFunction,placeList,coordsList,zoomList,noL=noLabels, dt = list(vaxBGm), rE=raceEthnicity,date=max(vaxBGs$VACCINATION_DATE),mt=mapType)

  if(bgC == TRUE){
    setwd(wd)
    bgComm <- read_excel("bgCommunity.xlsx")
    bgComm$GEOID <- as.character(bgComm$GEOID)
    vaxBGm <- merge(vaxBGm,bgComm,by="GEOID")
    
    vaxComm <- vaxBGm %>% data.table() %>% select(c(pop,n,increase,community))
    vaxComm <- vaxComm[,lapply(.SD,sum),by=.(community)]
    names(vaxComm) <- c("community","pop","n","increase")
    vaxComm$percent <- vaxComm$n/vaxComm$pop
    vaxComm$need <- vaxComm$pop - vaxComm$n
    vaxComm$raceEthnicity <- raceEthnicity
    
    setwd(paste0(wd,"./maps"))
    write.csv(vaxComm,paste0("vaxByCommunity",raceEthnicity,max(vaxBGs$VACCINATION_DATE),".csv"))
    
    return(list(vaxBGm,vaxComm))
  }
  
  return(vaxBGm)
  
}

### Benton County maps

bentonVax <- lapply(list("Total","Latinx","AIAN","Asian","Black","NHPI","White"), vaxFunction,
       wd = "L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton",
       alertCounty = "Benton",
       shpCounty = 3,
       placeList = list("Benton","Corvallis"),
       coordsList = list(c(left = -123.82135, bottom = 44.28859, right = -123.07, top = 44.71759),
                         c(left = -123.35, bottom = 44.5, right = -123.22843, top = 44.62)),
       zoomList = list(rep(10,3),c(12,12,14)),
       bgC = TRUE)

### Lincoln County maps

lincolnVax <- lapply(list("Total","Latinx","AIAN","Asian","Black","NHPI","White"), vaxFunction,
       wd = "L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/lincoln",
       alertCounty = "Lincoln",
       shpCounty = 41,
       placeList = list("Lincoln","Newport","LCity"),
       coordsList = list(c(left = -124.16130, bottom = 44.25869, right = -123.57903, top = 45.07131),
                         c(left = -124.09110, bottom = 44.54155, right = -123.9, top = 44.70307),
                         c(left = -124.05831, bottom = 44.8, right = -123.77, top = 45.075)),
       zoomList = list(rep(10,3),c(12,13,12),c(12,13,12)),
       mapType="watercolor")

### Linn County maps
linnVax <- lapply(list("Total","Latinx","AIAN","Asian","Black","NHPI","White"), vaxFunction,
       wd = "L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/linn",
       alertCounty = "Linn",
       shpCounty = 43,
       placeList = list("Linn","Albany","Lebanon","AlbanyDetail"),
       coordsList = list(c(left = -123.30696, bottom = 44.18669, right = -121.77986, top = 44.80771),
                         c(left = -123.2, bottom = 44.55, right = -122.99048, top = 44.67),
                         c(left = -123.02, bottom = 44.46767, right = -122.80319, top = 44.605),
                         c(left = -123.13505, bottom = 44.59459, right = -123.05509, top = 44.64104)),
       zoomList = list(rep(10,3),c(12,13,12),c(12,13,12),c(15,15)),
       noLabels = list(FALSE,FALSE,FALSE,TRUE),
       bgC = TRUE)

####################################
####################################
########## Communities
####################################
####################################

# Benton

bV <- unlist(bentonVax,recursive=FALSE)
bV <- split(bV, rep(1:2, length = length(bV)))
bV2 <- bV$"2"
BentonCommVax <- do.call(rbind.data.frame, bV2)

names(BentonCommVax) <- c("Community","Population","Number Vaxxed","Increase from previous week","Vax rate","Unmet need","Race Ethnicity")

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton/maps")
write_csv(BentonCommVax,paste0("bentonCommVax",Sys.Date(),".csv"))

# Linn

lnV <- unlist(linnVax,recursive=FALSE)
lnV <- split(lnV, rep(1:2, length = length(lnV)))
lnV2 <- lnV$"2"
LinnCommVax <- do.call(rbind.data.frame, lnV2)
names(LinnCommVax) <- c("Community","Population","Number Vaxxed","Increase from previous week","Vax rate","Unmet need","Race Ethnicity")

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/linn/maps")
write_csv(LinnCommVax,paste0("linnCommVax",Sys.Date(),".csv"))


####################################
####################################
########## Schools
####################################
####################################

bentonVaxTable <- data.table(do.call(rbind.data.frame,bV$"1"))
linnVaxTable <- data.table(do.call(rbind.data.frame,lnV$"1"))

lblVaxTable <- rbind(bentonVaxTable,linnVaxTable) %>% select(GEOID,raceEthnicity,pop,n,percent,stillNeeds,increase)

lblVaxShp <- merge(shp,lblVaxTable,by="GEOID")

schoolVax <- lapply(list("Total","Latinx","AIAN","Asian","Black","NHPI","White"),schoolFunction,dt=lblVaxShp,district =)

schoolVax <- data.table(do.call(rbind.data.frame,schoolVax))

setwd("L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton")

write.csv(schoolVax,paste0("schoolVax",Sys.Date(),".csv"))

# ## hard code for testing
# wd <- "L:/Health/CD/Coronavirus - Covid19 - SARS CoV 2/Epi Data/Data - Orpheus Exports/banwarth/raw/alertData/benton"
# alertCounty <- "Benton"
# raceEthnicity <- "Total"
# shpCounty <- 3
# placeList <- list("Benton","Corvallis")
# coordsList <- list(c(left = -123.82135, bottom = 44.28859, right = -123.07, top = 44.71759),
#                    c(left = -123.34464, bottom = 44.5, right = -123.22843, top = 44.62))
# zoomList <- list(rep(10,3),c(12,12,14))
# noLabels <- FALSE
# mapType <- "toner-background"


# https://r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
# https://mattherman.info/blog/point-in-poly/
