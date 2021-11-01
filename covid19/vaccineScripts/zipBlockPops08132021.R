# updated in BentonOpera 08162021

library(sf)
library(tidyverse)
library(tidycensus)
library(ggplot2)
library(data.table)
library(ggmap)
library(RColorBrewer)
library(lwgeom)

# reading in shapefiles
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")

header <- read.delim("orgeo2020.pl", header=FALSE, colClasses="character", sep="|")
part1  <- read.delim("or000012020.pl",  header=FALSE, colClasses="character", sep="|")
part2  <- read.delim("or000022020.pl",  header=FALSE, colClasses="character", sep="|")
part3  <- read.delim("or000032020.pl",  header=FALSE, colClasses="character", sep="|")


# -----------------------------
colnames(header) <- c("FILEID", "STUSAB", "SUMLEV", "GEOVAR", "GEOCOMP", "CHARITER", "CIFSN", "LOGRECNO", "GEOID",
                      "GEOCODE", "REGION", "DIVISION", "STATE", "STATENS", "COUNTY", "COUNTYCC", "COUNTYNS", "COUSUB",
                      "COUSUBCC", "COUSUBNS", "SUBMCD", "SUBMCDCC", "SUBMCDNS", "ESTATE", "ESTATECC", "ESTATENS",
                      "CONCIT", "CONCITCC", "CONCITNS", "PLACE", "PLACECC", "PLACENS", "TRACT", "BLKGRP", "BLOCK",
                      "AIANHH", "AIHHTLI", "AIANHHFP", "AIANHHCC", "AIANHHNS", "AITS", "AITSFP", "AITSCC", "AITSNS",
                      "TTRACT", "TBLKGRP", "ANRC", "ANRCCC", "ANRCNS", "CBSA", "MEMI", "CSA", "METDIV", "NECTA",
                      "NMEMI", "CNECTA", "NECTADIV", "CBSAPCI", "NECTAPCI", "UA", "UATYPE", "UR", "CD116", "CD118",
                      "CD119", "CD120", "CD121", "SLDU18", "SLDU22", "SLDU24", "SLDU26", "SLDU28", "SLDL18", "SLDL22",
                      "SLDL24", "SLDL26", "SLDL28", "VTD", "VTDI", "ZCTA", "SDELM", "SDSEC", "SDUNI", "PUMA", "AREALAND",
                      "AREAWATR", "BASENAME", "NAME", "FUNCSTAT", "GCUNI", "POP100", "HU100", "INTPTLAT", "INTPTLON",
                      "LSADC", "PARTFLAG", "UGA")
colnames(part1) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
                     paste0("P00", c(10001:10071, 20001:20073)))
colnames(part2) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
                     paste0("P00", c(30001:30071, 40001:40073)),
                     paste0("H00", 10001:10003))
colnames(part3) <- c("FILEID", "STUSAB", "CHARITER", "CIFSN", "LOGRECNO",
                     paste0("P00", 50001:50010))

# -----------------------------
# Merge the data
# -----------------------------
combine <- Reduce(function(x,y) {merge(x, y, by=c("LOGRECNO", "STUSAB", "FILEID", "CHARITER"))}, list(header[,-7], part1[,-4], part2[,-4], part3))

# -----------------------------
# Order the data
# -----------------------------
combine <- combine[order(combine$LOGRECNO), c("FILEID", "STUSAB", "SUMLEV", "GEOVAR", "GEOCOMP", "CHARITER", "CIFSN", "LOGRECNO", "GEOID",
                                              "GEOCODE", "REGION", "DIVISION", "STATE", "STATENS", "COUNTY", "COUNTYCC", "COUNTYNS", "COUSUB",
                                              "COUSUBCC", "COUSUBNS", "SUBMCD", "SUBMCDCC", "SUBMCDNS", "ESTATE", "ESTATECC", "ESTATENS",
                                              "CONCIT", "CONCITCC", "CONCITNS", "PLACE", "PLACECC", "PLACENS", "TRACT", "BLKGRP", "BLOCK",
                                              "AIANHH", "AIHHTLI", "AIANHHFP", "AIANHHCC", "AIANHHNS", "AITS", "AITSFP", "AITSCC", "AITSNS",
                                              "TTRACT", "TBLKGRP", "ANRC", "ANRCCC", "ANRCNS", "CBSA", "MEMI", "CSA", "METDIV", "NECTA",
                                              "NMEMI", "CNECTA", "NECTADIV", "CBSAPCI", "NECTAPCI", "UA", "UATYPE", "UR", "CD116", "CD118",
                                              "CD119", "CD120", "CD121", "SLDU18", "SLDU22", "SLDU24", "SLDU26", "SLDU28", "SLDL18", "SLDL22",
                                              "SLDL24", "SLDL26", "SLDL28", "VTD", "VTDI", "ZCTA", "SDELM", "SDSEC", "SDUNI", "PUMA", "AREALAND",
                                              "AREAWATR", "BASENAME", "NAME", "FUNCSTAT", "GCUNI", "POP100", "HU100", "INTPTLAT", "INTPTLON",
                                              "LSADC", "PARTFLAG", "UGA", paste0("P00", c(10001:10071, 20001:20073)), paste0("P00", c(30001:30071, 40001:40073)),
                                              paste0("H00", 10001:10003), paste0("P00", 50001:50010))]
rownames(combine) <- 1:nrow(combine)
combine[,c(98:ncol(combine))] <- apply(combine[,c(98:ncol(combine))],2,as.numeric)

lbl <- combine[combine$COUNTY %in% c("003","041","043"),]

# extract county and desired population fields

# Benton County
benton <- lbl[lbl$COUNTY == "003",c("GEOCODE",paste0("P00", c(10001:10071, 20002)))]
names(benton)[1] <- "GEOID20"

bshp <- st_read("tl_2020_41_tabblock20.shp")

bentonBlocks <- merge(bshp,benton,by="GEOID20")

zips <- get_acs(geography="zcta",
                variables="B01001_001",
                year=2019,
                state="OR",
                geometry="True")
zips <- zips[which(zips$GEOID %in% c("97330","97333","97370","97321","97324","97326","97331","97456")),]

blockZips <- st_join(bentonBlocks,zips,largest=TRUE)

# sum over zip code
zipPop <- as.data.table(blockZips)[,c("GEOID",paste0("P00", c(10001:10071, 20002)))][,lapply(.SD,sum),by=.(GEOID)]

# create columns for race/ethnicity
names(zipPop)[c(2,ncol(zipPop))] <- c("Total","Latinx")
zipPop$White <- rowSums(zipPop[,paste0("P00",c(10003,10011:10015,10027:10036,10048:10057,10064:10068,10071))])
zipPop$Black <- rowSums(zipPop[,paste0("P00",c(10004,10011,10016:10019,10027:10030,10037:10042,10048:10053,10058:10061,10064:10067,10069,10071))])
zipPop$AIAN <- rowSums(zipPop[,paste0("P00",c(10005,10012,10016,10020:10022,10027,10031:10033,10037:10039,10043,10045,10048:10050,10054:10056,10058:10060,10062,10064:10066,10068,10069,10071))])
zipPop$Asian <- rowSums(zipPop[,paste0("P00",c(10006,10013,10017,10020,10023,10024,10028,10031,10034,10035,10037,10040,10041,10043,10044,10046,10048,10051,10052,10054,10055,10057:10059,10061,10062,10064,10065,10067:10069,10071))])
zipPop$NHPI <- rowSums(zipPop[,paste0("P00",c(10007,10014,10018,10021,10023,10025,10029,10032,10034,10036,10038,10040,10042,10043,10045,10046,10049,10051,10053,10054,10056:10058,10060,10061,10062,10064,10066:10069,10071))])

write.csv(zipPop[,c("GEOID","Total","Latinx","White","Black","AIAN","Asian","NHPI")],file="zipPopBenton.csv", row.names=FALSE)



# Lincoln County
lincoln <- lbl[lbl$COUNTY == "041",c("GEOCODE",paste0("P00", c(10001:10071, 20002)))]
names(lincoln)[1] <- "GEOID20"

lshp <- st_read("tl_2020_41_tabblock20.shp")

lincolnBlocks <- merge(lshp,lincoln,by="GEOID20")

zips <- get_acs(geography="zcta",
                variables="B01001_001",
                year=2019,
                state="OR",
                geometry="True")
zips <- zips[which(zips$GEOID %in% c("97341",
                                     "97343",
                                     "97357",
                                     "97364",
                                     "97365",
                                     "97366",
                                     "97367",
                                     "97368",
                                     "97369",
                                     "97376",
                                     "97380",
                                     "97388",
                                     "97391",
                                     "97394",
                                     "97498")),]

blockZips <- st_join(lincolnBlocks,zips,largest=TRUE)

# sum over zip code
zipPop <- as.data.table(blockZips)[,c("GEOID",paste0("P00", c(10001:10071, 20002)))][,lapply(.SD,sum),by=.(GEOID)]

# create columns for race/ethnicity
names(zipPop)[c(2,ncol(zipPop))] <- c("Total","Latinx")
zipPop$White <- rowSums(zipPop[,paste0("P00",c(10003,10011:10015,10027:10036,10048:10057,10064:10068,10071))])
zipPop$Black <- rowSums(zipPop[,paste0("P00",c(10004,10011,10016:10019,10027:10030,10037:10042,10048:10053,10058:10061,10064:10067,10069,10071))])
zipPop$AIAN <- rowSums(zipPop[,paste0("P00",c(10005,10012,10016,10020:10022,10027,10031:10033,10037:10039,10043,10045,10048:10050,10054:10056,10058:10060,10062,10064:10066,10068,10069,10071))])
zipPop$Asian <- rowSums(zipPop[,paste0("P00",c(10006,10013,10017,10020,10023,10024,10028,10031,10034,10035,10037,10040,10041,10043,10044,10046,10048,10051,10052,10054,10055,10057:10059,10061,10062,10064,10065,10067:10069,10071))])
zipPop$NHPI <- rowSums(zipPop[,paste0("P00",c(10007,10014,10018,10021,10023,10025,10029,10032,10034,10036,10038,10040,10042,10043,10045,10046,10049,10051,10053,10054,10056:10058,10060,10061,10062,10064,10066:10069,10071))])

write.csv(zipPop[,c("GEOID","Total","Latinx","White","Black","AIAN","Asian","NHPI")],file="zipPopLincoln.csv", row.names=FALSE)

# assign blocks w/0 zips to zip codes based on proximity
