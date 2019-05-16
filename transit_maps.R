# transit_maps.R
# Access to public transit


# ArcMap Procedure
# 
# 1. Select desired transit stops and create new layer from selection.
# 2. Export selection as shape file.
# 3. Geoprocessing - buffer: Choose layer, name buffer, set distance, Dissolve = All.
# 4. Add Census blocks to map.
# 5. Clip Census blocks with all-transit buffer layer
# 6. Union clipped blocks with all blocks, save as new shapefile and add layer
# 7. Arc Toolbox:Data Management:Features:Multipart to Singlepart to separate non-contiguous blocks.
# 8. Delete extraneous fields
# 9. Open data table of clipped/unioned census blocks
# 10. Add fields: ClipArea, Latitude, Longitude as doubles
# 11. Calculate geometry: ClipArea = Area in square meters
# 12. Calculate geometry: Coordinates of centroids
# 13. Clip unioned layer with desired buffer, save as new shapefile; repeae for all different route types



# R Procedure
# 1. Read in shapefile with all transit stops
# 2. Calculate PctArea from the two ALANDs and ClipArea
# 3. Merge block groups
# 4. Sum Pop10 over block group
# 5. Divide Pop10 by bgPop10 to get percent of population in each block group
# 6. Compute BufferPop17 etc. by (Pop10/bgPop10)*(Pop17)*PctArea
## 7. Compute ExcludePop17 ect. by (Pop10/bgPop10)*(Pop17)*(1-PctArea)
## 8. Run bap subroutine - different commands for all transit stops compared to other types
# 9. Drop all rows where PctArea == 0
# 10. Compute percent of population within 1/4 mile of any transit stop
# 11. Repeat steps c(1:9)-[6] for other types of transit stops.

# power bi tips
# No relationships between bap_transit and pop_type_table
# FilterMeasure = 
# VAR ranks = max(pop_type_table[pop_rank])
# RETURN(IF(COUNTROWS(FILTER(bap_transit,bap_transit[bap_rank]<=ranks))>0,1,0))

setwd("L:/Health/Epidemiology/Banwarth_Epi/GIS/GIS_files/lbl_transit")

library(maptools)
library(rgdal)
library(data.table)
library(reshape2)
library(plyr)

# read in common datasets

# read in most recent ACS blockgroup populations
bg_transit <- read.table('blockgroup_lbl_transit_2017.txt',sep=',',header=TRUE)
bg_transit$COUNTYFP10 <- as.numeric(substr(as.character(bg_transit$blockgroup),start=3,stop=5))

# all blockgroups population in 2010
bg_pop10 <- read.csv('lbl_bg_pop10.csv',header=TRUE,sep=',')
names(bg_pop10) <- c('blockgroup','bgPop10')

# city blocks
city_blocks <- read.csv('city_blocks.csv',header=TRUE,sep=",")
city_blocks$urban <- 1
city_blocks$blockgroup <- round(as.numeric(as.character(city_blocks$GEOID10))*1e-3,0)
city_blockgroups <- as.data.table(city_blocks[,c('POP10','blockgroup')])[,lapply(.SD,sum),by=.(blockgroup)]
names(city_blockgroups) <- c('blockgroup','urbanPop10')

# data set with variable names and ranks
variable_df <- read.csv('variable_list.csv',sep=',',header=TRUE)

# Extracting lists of names
bg_colnames <- names(bg_transit)[-c(14,16)]
block_colnames <- paste('block',bg_colnames,sep="")
buffer_colnames <- paste('Buffer',bg_colnames,sep="")



# identifying city blockgroups
bg_transit <- merge(bg_transit,city_blockgroups,by='blockgroup',all.x=TRUE)
bg_transit[is.na(bg_transit$urbanPop10),'urbanPop10'] <- 0
bg_transit <- merge(bg_transit,bg_pop10,by='blockgroup')

# calculating percent of population in urban part of blockgroups
bg_urban_colnames <- paste('urban',bg_colnames,sep="")
bg_transit[,bg_urban_colnames] <- round(bg_transit[,bg_colnames]*bg_transit$urbanPop10/bg_transit$bgPop10,0)

# total populations across regions for different population group types
total_pops <- colSums(bg_transit[,bg_colnames])
pop_pcts <- total_pops/total_pops[1]
variable_df$pop_pct <- pop_pcts
variable_df$pop_rank <- rank(pop_pcts)
variable_df <- variable_df[order(variable_df$pop_rank,decreasing=TRUE),]
write.table(variable_df,'variable_table.txt',row.names=FALSE,sep=',')

# Define coordinate system object
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# number of rows in bap table
baprows <- 20000

# functions

# transit_readin
transit_readin <- function(shapefn,transit_type,typerank){
  transit_df <- data.frame(rgdal::readOGR(shapefn))
  transit_df$FID_lbl_bl <- as.numeric(as.character(transit_df$FID_lbl_bl))
  
  if(typerank == 1){
    transit_df$ALAND10 <- pmax(as.numeric(as.character(transit_df$ALAND10)),as.numeric(as.character(transit_df$ALAND10_1)),na.rm=TRUE)
    transit_df$PctArea <- transit_df$ClipArea/transit_df$ALAND10
    transit_df$PctArea[which(transit_df$PctArea==Inf)] <- 0
    transit_df$GEOID10 <- pmax(as.numeric(as.character(transit_df$GEOID10)),as.numeric(as.character(transit_df$GEOID10_1)),na.rm = TRUE)
    transit_df$POP10 <- pmax(as.numeric(as.character(transit_df$POP10)),as.numeric(as.character(transit_df$POP10_1)),na.rm=TRUE)
    transit_df$COUNTYFP10 <- pmax(as.numeric(as.character(transit_df$COUNTYFP10)),as.numeric(as.character(transit_df$COUNTYFP_1)),na.rm=TRUE)
    transit_df[,transit_type] <- 0
    transit_df[transit_df$FID_lbl_bl!=-1,transit_type] <- 1
    transit_df <- transit_df[,c('FID_lbl_bl','GEOID10','ALAND10','POP10','ClipArea','Latitude','Longitude','PctArea',transit_type)]
    return(transit_df)
  } else {
    transit_df <- data.frame(transit_df[,'FID_lbl_bl'],rep(1,nrow(transit_df)))
    names(transit_df) <- c('FID_lbl_bl',transit_type)
    transit_df <- merge(tcomb,transit_df,by='FID_lbl_bl',all.x=TRUE)
    transit_df[is.na(transit_df)] <- 0
    return(transit_df)
  }
}


# Function to add/subtract population to buffers.
addpop_function <- function(x,tdf,bgdf) {
  x1 <- tcomb[,x]
  x2 <- bg_transit[,sub('Buffer','',x)]
  tadd <- sum(x2)-sum(x1)
  xsamp <- sample(length(x1[x1>0]),abs(tadd))
  if(tadd > 0) {
    x1[x1>0][xsamp] <- x1[x1>0][xsamp]+1
  } else if(tadd < 0) {
    x1[x1>0][xsamp] <- x1[x1>0][xsamp]-1
  }
  return(x1)
}


# function for flagging which population type
poptype_fun <- function(sbpname,bapstart,bapend) {
  flag <- rownames(bapstart[rep(row.names(bapstart),bapstart[,sbpname]),])
  exp_df <- data.frame(flag,rep(1,length(flag)))
  names(exp_df) <- c('RowNames',paste(sbpname,'Flag',sep=""))
  bapend <- merge(bapend,exp_df,by='RowNames',all.x=TRUE)
  return(bapend)
}





### Step 1: Read in transit shapefiles and combine

routelist <- c('AllRoutes','Local','Intercity','Hourly','Intermittent','Infrequent')

# transit_readin appends each tcomb to make a total tcomb=transit_combined
tcomb <- transit_readin('lbl_blocks_union.shp','AllRoutes',1)
tcomb <- transit_readin('lbl_blocks_loca.shp','Local',2)
tcomb <- transit_readin('lbl_blocks_intercity.shp','Intercity',2)
tcomb <- transit_readin('lbl_blocks_hourly.shp','Hourly',2)
tcomb <- transit_readin('lbl_blocks_intermittent.shp','Intermittent',2)
tcomb <- transit_readin('lbl_blocks_infrequent.shp','Infrequent',2)


# Step 2: calculate block and buffer populations

# identifying blockgroups
tcomb$blockgroup <- round(as.numeric(as.character(tcomb$GEOID10))*1e-3,0)

# merge data sets
tcomb <- merge(tcomb,bg_transit,by='blockgroup')
# tcomb <- merge(tcomb,bg_pop10,by='blockgroup')

# divide Pop10 by blockgroup pop: blockPop10Pct to get each block's share of the blockgroup population
tcomb$blockPop10Pct <- tcomb$POP10/tcomb$bgPop10
tcomb$blockPop10Pct[is.na(tcomb$blockPop10Pct)] <- 0

# multiply Pop17 etc. by blockPop10Pct= blockPop17 etc. to get each block's Pop17 etc.
tcomb[,block_colnames] <- tcomb[,bg_colnames]*tcomb$blockPop10Pct
tcomb[which(tcomb$bgPop10==0),block_colnames] <- 0

# multiply Pop17block etc. by PctArea = BufferPop17 etc.
tcomb[,buffer_colnames] <- round(tcomb[,block_colnames]*tcomb$PctArea,0)

# randomly assigning extra people to get to total population
tcomb[,buffer_colnames] <- as.data.frame(lapply(buffer_colnames,function(x) addpop_function(x)))




# Step 3: Compute blocks as points
bap <- tcomb[,c('FID_lbl_bl','Latitude','Longitude','COUNTYFP10',routelist,buffer_colnames)]

# one row for each person
bap_exp <- bap[rep(row.names(bap),bap$BufferPop17),]
bap_exp$RowNames <- rownames(bap_exp)
bap_exp$BufferPop17Flag <- 1

# indicating population type
subpop_names <- buffer_colnames[-1]
for(i in 1:length(subpop_names)) {
  bap_exp <- poptype_fun(subpop_names[i],bap,bap_exp)
}
bap_exp[is.na(bap_exp)] <- 0

# randomizing location
bap_exp$Latitude <- bap_exp$Latitude+runif(nrow(bap_exp),0,.01)
bap_exp$Longitude <- bap_exp$Longitude+runif(nrow(bap_exp),0,.01)

# reordering/dropping columns
bap_exp <- bap_exp[,c('FID_lbl_bl','Latitude','Longitude','COUNTYFP10',routelist,paste(buffer_colnames,'Flag',sep="")),]
names(bap_exp) <- c('FID_lbl_bl','Latitude','Longitude','COUNTYFP10',routelist,bg_colnames)

bap_exp$geography <- mapvalues(bap_exp$COUNTYFP10,from=c(3,41,43),to=c('Benton County','Lincoln County','Linn County'))

# generating bap_rank column
bap_exp$bap_rank <- 0

for(i in 1:length(bg_colnames)) {
  bapcolname <- sub('Buffer','',variable_df[i,'variable'])
  bap_exp[bap_exp[,bapcolname]==1,'bap_rank'] <- length(bg_colnames)+1-i
}

# sample to reduce size
bap_exp <-  bap_exp[sample(nrow(bap_exp),baprows),]

# replicating bap_exp; one for each route type
bap_exp <- bap_exp[rep(row.names(bap_exp),length(routelist)),]
bap_exp$Access <- 0
bap_exp$'Route type' <- rep(c('All routes','Local routes','Intercity routes',
                              'Hourly routes','Intermittent routes','Infrequent routes'),each=baprows)
bap_exp$route_rank <- rep(c(1,2,3,4,5,6),each=baprows)
for(i in 1:length(routelist)) {
  j <- (i-1)*baprows+1
  k <- i*baprows
  bap_exp$Access[c(j:k)] <- bap_exp[c(j:k),routelist[i]]
}

# keeping only useful columns
bap_exp <- bap_exp[,c('geography','Latitude','Longitude','bap_rank','Access','Route type','route_rank')]

write.table(variable_df,'pop_type_table.txt',row.names=FALSE,sep=',')
write.table(bap_exp,'bap_transit.txt',row.names = FALSE,sep=',')

rm(bap_exp)
rm(bap)
rm(bapcolname)


# Step 4: compute percent of population w/in 1/4 mile: BufferPop17Pct

# need to use bg_transit in denominator because tcomb only has blocks that overlap the buffer

# Dropping all blocks that are outside of the buffer to match the "All transit" data frame with the other routes data frames.
tcomb <- tcomb[which(tcomb$AllRoutes == 1),]

bufferpct_function <- function(transitCol,transitName,transitRank) {
  # keeping only desired transit_type
  tdf <- tcomb[tcomb[,transitCol]==1,]
  
  # three county region
  regionBufferPops <- colSums(tdf[,buffer_colnames])
  regionTotalPops <- colSums(bg_transit[,bg_colnames])
  regionBufferPcts <- regionBufferPops/regionTotalPops
  
  # each county
  countyBufferPops <- as.data.table(tdf[,c('COUNTYFP10',buffer_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
  countyTotalPops <- as.data.table(bg_transit[,c('COUNTYFP10',bg_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
  countyTotalPops <- countyTotalPops[which(countyTotalPops$COUNTYFP10 %in% unique(countyBufferPops$COUNTYFP10)),]
  countyBufferPcts <- countyBufferPops[,-1]/countyTotalPops[,-1]
  
  # three county urban/rural
  # city blocks
  tdf <- merge(tdf,city_blocks[,c('GEOID10','urban')],by='GEOID10',all.x=TRUE)
  tdf[is.na(tdf$urban),'urban'] <- 0
  
  # whole region urban/rural
  regionUrbanBufferPops <- as.data.table(tdf[,c('urban',buffer_colnames)])[,lapply(.SD,sum),by=.(urban)]
  regionUrbanTotalPops <- colSums(bg_transit[,bg_urban_colnames])
  regionUrbanTotalPops <- rbind(regionUrbanTotalPops,regionTotalPops-regionUrbanTotalPops)
  regionUrbanBufferPcts <- regionUrbanBufferPops[,-1]/regionUrbanTotalPops
  regionUrbanTotalPops <- data.frame(regionUrbanTotalPops)
  
  # each county urban/rural
  countyUrbanBufferPops <- as.data.table(tdf[,c('COUNTYFP10','urban',buffer_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10,urban)]
  countyUrbanBufferPops <- countyUrbanBufferPops[order(countyUrbanBufferPops$COUNTYFP10,decreasing=FALSE),]
  countyUrbanBufferPops <- countyUrbanBufferPops[order(countyUrbanBufferPops$urban,decreasing=TRUE),]
  countyUrbanTotalPops <- as.data.table(bg_transit[,c('COUNTYFP10',bg_urban_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
  countyUrbanTotalPops <- countyUrbanTotalPops[which(countyUrbanTotalPops$COUNTYFP10 %in% unique(countyUrbanBufferPops$COUNTYFP10)),]
  countyRuralTotalPops <- cbind(countyUrbanTotalPops$COUNTYFP10,countyTotalPops[,-1]-countyUrbanTotalPops[,-1])
  countyUrbanTotalPops <- rbind(countyUrbanTotalPops,countyRuralTotalPops,use.names=FALSE)
  countyUrbanBufferPcts <- countyUrbanBufferPops[,-c(1,2)]/countyUrbanTotalPops[,-1]
  
  # making one datatable
  region_df <- data.frame(cbind(regionBufferPops,regionTotalPops,regionBufferPcts))
  colnames(region_df) <- c('BufferPop','TotalPop','BufferPct')
  region_df$variable <- rownames(region_df)
  countyBufferPcts$COUNTYFP10 <- as.factor(unique(countyBufferPops$COUNTYFP10))
  county_df <- data.frame(cbind(data.table::melt(countyBufferPops,id.vars='COUNTYFP10'),
                     data.table::melt(countyTotalPops,id.vars='COUNTYFP10'),
                     data.table::melt(countyBufferPcts,id.vars='COUNTYFP10')))
  names(county_df) <- c('geography','variable','BufferPop','a','b','TotalPop','c','d','BufferPct')
  county_df$geography <- mapvalues(county_df$geography,from=c(3,41,43),to=c('Benton County','Lincoln County','Linn County'))
  region_df$geography <- as.factor('Region')
  region_df$urban <- 2
  county_df$urban <- 2
  
  regionUrbanTotalPops$urban <- as.factor(c(1,0))
  regionUrbanBufferPcts$urban <- as.factor(c(1,0))
  region_udf <- data.frame(cbind(data.table::melt(regionUrbanBufferPops,id.vars='urban'),
                           data.table::melt(regionUrbanTotalPops,id.vars='urban'),
                           data.table::melt(regionUrbanBufferPcts,id.vars='urban')))
  names(region_udf) <- c('urban','variable','BufferPop','a','b','TotalPop','c','d','BufferPct')
  region_udf$geography <- 'Region'
  
  names(countyUrbanTotalPops)[1] <- 'COUNTYFP10'
  countyUrbanTotalPops$urban <- as.factor(rep(c(1,0),each=length(unique(countyUrbanBufferPops$COUNTYFP10))))
  countyUrbanBufferPcts$COUNTYFP10 <- as.factor(rep(unique(countyUrbanBufferPops$COUNTYFP10),2))
  countyUrbanBufferPcts$urban <- as.factor(rep(c(1,0),each=length(unique(countyUrbanBufferPops$COUNTYFP10))))
  county_udf <- data.frame(cbind(melt(countyUrbanBufferPops,id.vars=c('COUNTYFP10','urban')),
                                 melt(countyUrbanTotalPops,id.vars=c('COUNTYFP10','urban')),
                                 melt(countyUrbanBufferPcts,id.vars=c('COUNTYFP10','urban'))))
  names(county_udf) <- c('geography','urban','variable','BufferPop','a','b','c','TotalPop','d','e','f','BufferPct')
  county_udf$geography <- mapvalues(county_udf$geography,from=c(3,41,43),to=c('Benton County','Lincoln County','Linn County'))
  
  df <- data.frame(rbind(region_df[,c('geography','urban','variable','BufferPop','TotalPop','BufferPct')],
                         county_df[,c('geography','urban','variable','BufferPop','TotalPop','BufferPct')],
                         region_udf[,c('geography','urban','variable','BufferPop','TotalPop','BufferPct')],
                         county_udf[,c('geography','urban','variable','BufferPop','TotalPop','BufferPct')]))
  
  df$urban <- mapvalues(df$urban,from=c(0,1,2),to=c('Rural','City','All communities'))
  df$type <- transitName
  df$typerank <- transitRank
  
  df <- merge(df,variable_df,by='variable')
  
  return(df)
}

alldf <- bufferpct_function('AllRoutes','All routes',1)
localdf <- bufferpct_function('Local','Local routes',2)
intercitydf <- bufferpct_function('Intercity','Intercity routes',3)
hourlydf <- bufferpct_function('Hourly','Hourly routes',4)
intermittentdf <- bufferpct_function('Intermittent','Intermittent routes',5)
infrequentdf <- bufferpct_function('Infrequent','Infrequent routes',6)

df <- data.frame(rbind(alldf,localdf,intercitydf,hourlydf,intermittentdf,infrequentdf))
write.csv(df,'transit_access.csv',row.names=FALSE)

