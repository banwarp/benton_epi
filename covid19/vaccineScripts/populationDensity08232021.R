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

bshp <- st_read("tl_2020_41_tabblock20.shp")

mapFunction <- function(rE,pL,crds,zL,noL=FALSE,dt,mt,wd) {
  
  setwd(paste0(wd,"/popDensityMaps"))
  
  dt <- dt[,which(names(dt)==rE)]
  names(dt)[1] <- "density"
  
  quants <- quantile(dt$density,c(.3,.5,.7,.85,1),na.rm=TRUE)
  dt$DensityPercentile  <- with(dt, factor(ifelse(density < quants[1], 0.3, 
                                        ifelse(density < quants[2], 0.5,
                                         ifelse(density < quants[3], 0.7,
                                          ifelse(density < quants[4], 0.85,
                                           ifelse(density < quants[5], 1, NA)))))))
  
  # mapping County block groups
  basemap <- get_stamenmap(crds, zoom = zL[1], maptype = mt)
  lines <- get_stamenmap(crds, zoom = zL[2], maptype = "toner-lines")
  if(!noL) {labels <- get_stamenmap(crds, zoom = zL[3], maptype = "toner-labels")}
  
  m1 <- ggmap(basemap) +
    geom_sf(data = dt, aes(fill = DensityPercentile, alpha = .1),inherit.aes=FALSE)+coord_sf(datum=NA) +
    scale_fill_manual(values = c(brewer.pal(5,"Spectral"))) +
    labs(title=paste(pL,rE,"Population density",sep = " "))+
    theme(title=element_text(size=8), legend.text=element_text(size=6),axis.text=element_text(size=6))+
    inset_ggmap(lines)
  
  if(!noL){m1 <- m1 + inset_ggmap(labels)}
  
  ggsave(paste0(pL,rE,"PopDens.jpg"), plot = m1, width = 6, height = 6, units = "in")
  
}


# extract county and desired population fields
densityFunction <- function(wd,countyName,countyCode,placeList,coordsList,zoomList,shpFile=bshp, noLabels=FALSE,mapType="toner-background",repEach) {
  
  cty <- lbl[lbl$COUNTY == countyCode,c("GEOCODE","AREALAND",paste0("P00", c(10001:10071, 20002)))]
  names(cty)[1] <- "GEOID20"
  
  ctyBlocks <- merge(shpFile,cty,by="GEOID20")
  
  # create columns for race/ethnicity
  ctyBlocks$White <- rowSums(data.table(ctyBlocks)[,paste0("P00",c(10003,10011:10015,10027:10036,10048:10057,10064:10068,10071))])
  ctyBlocks$Black <- rowSums(data.table(ctyBlocks)[,paste0("P00",c(10004,10011,10016:10019,10027:10030,10037:10042,10048:10053,10058:10061,10064:10067,10069,10071))])
  ctyBlocks$AIAN <- rowSums(data.table(ctyBlocks)[,paste0("P00",c(10005,10012,10016,10020:10022,10027,10031:10033,10037:10039,10043,10045,10048:10050,10054:10056,10058:10060,10062,10064:10066,10068,10069,10071))])
  ctyBlocks$Asian <- rowSums(data.table(ctyBlocks)[,paste0("P00",c(10006,10013,10017,10020,10023,10024,10028,10031,10034,10035,10037,10040,10041,10043,10044,10046,10048,10051,10052,10054,10055,10057:10059,10061,10062,10064,10065,10067:10069,10071))])
  ctyBlocks$NHPI <- rowSums(data.table(ctyBlocks)[,paste0("P00",c(10007,10014,10018,10021,10023,10025,10029,10032,10034,10036,10038,10040,10042,10043,10045,10046,10049,10051,10053,10054,10056:10058,10060,10061,10062,10064,10066:10069,10071))])
  
  ctyBlocks <- ctyBlocks[,c("GEOID20","AREALAND","P0010001","P0020002","White","Black","AIAN","Asian","NHPI")]
  ctyBlocks$AREALAND <- as.numeric(ctyBlocks$AREALAND)
  names(ctyBlocks)[c(1:4)] <- c("GEOID","Area","Total","Latinx")
  
  ctyDensity <- data.table(apply(data.table(ctyBlocks)[,c(3:9)],2,function(x) as.numeric(x)/ctyBlocks$Area)[,c(1:7)])
  ctyDensity$GEOID <- ctyBlocks$GEOID
  
  ctyBlocks <- merge(ctyBlocks,ctyDensity,by="GEOID")
  ctyBlocks[ctyBlocks$Total.x == 0,"Total.y"] <- NA
  ctyBlocks[ctyBlocks$Latinx.x == 0,"Latinx.y"] <- NA
  ctyBlocks[ctyBlocks$Black.x == 0,"Black.y"] <- NA
  ctyBlocks[ctyBlocks$AIAN.x == 0,"AIAN.y"] <- NA
  ctyBlocks[ctyBlocks$Asian.x == 0,"Asian.y"] <- NA
  ctyBlocks[ctyBlocks$NHPI.x == 0,"NHPI.y"] <- NA
  ctyBlocks[ctyBlocks$White.x == 0,"White.y"] <- NA
  
  
  mapply(mapFunction,rep(c("Total.y","Latinx.y","AIAN.y","Black.y","NHPI.y","White.y"),each=repEach),
          placeList,coordsList,zoomList,noL=TRUE,dt = list(ctyBlocks),mt=mapType,wd)
  
  
}

# Benton County
densityFunction("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020",
                "Benton",
                "003",
                list("Benton","Corvallis"),
                list(c(left = -123.82135, bottom = 44.28859, right = -123.07, top = 44.71759),
                     c(left = -123.35, bottom = 44.5, right = -123.22843, top = 44.62)),
                list(rep(10,3),c(12,12,14)),
                shpFile=bshp,
                noLabels=TRUE,
                mapType="toner-background",
                repEach=2)

# Lincoln County
densityFunction("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020",
                "Lincoln",
                "041",
                list("Lincoln","Newport","LCity"),
                list(c(left = -124.16130, bottom = 44.25869, right = -123.57903, top = 45.07131),
                     c(left = -124.09110, bottom = 44.54155, right = -123.9, top = 44.70307),
                     c(left = -124.05831, bottom = 44.8, right = -123.77, top = 45.075)),
                list(rep(10,3),c(12,13,12),c(12,13,12)),
                shpFile=bshp,
                noLabels=TRUE,
                mapType="toner-background",
                repEach=3)

# Linn County
densityFunction("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020",
                "Linn",
                "043",
                list("Linn","Albany","Lebanon","AlbanyDetail"),
                list(c(left = -123.30696, bottom = 44.18669, right = -121.77986, top = 44.80771),
                     c(left = -123.2, bottom = 44.55, right = -122.99048, top = 44.67),
                     c(left = -123.02, bottom = 44.46767, right = -122.80319, top = 44.605),
                     c(left = -123.13505, bottom = 44.59459, right = -123.05509, top = 44.64104)),
                list(rep(10,3),c(12,13,12),c(12,13,12),c(15,15)),
                shpFile=bshp,
                noLabels=TRUE,
                mapType="toner-background",
                repEach=4)

  
  # # for testing
  # wd <- setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/census2020")
  # countyName <- "Benton"
  # countyCode <- "003"
  # placeList <- list("Corvallis")
  # coordsList <- list(c(left = -123.34464, bottom = 44.5, right = -123.22843, top = 44.62))
  # zoomList <- list(c(12,12,14))
  # mapType <- "toner-background"
  # bshp <- st_read("tl_2020_41_tabblock20.shp"
  


