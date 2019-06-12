# parks_maps.R
# Access to public parks


# ArcMap Procedure
# 
# 1. Select desired parks and create new layer from selection.
# 2. Export selection as shape file.
# 3. Make a "parks_with_roads" buffer at 50 feet to exclude any roads running through the middle of parks from being included when mapping points
# 4. Clip all blocks using the parks_with_roads buffer. Union the excluded blocks back to the other blocks, but zero out the population in the partial blocks within the parks areas.
# 5. Geoprocessing - buffer: Choose layer, name buffer, set distance, Dissolve = All.
# 6. Clip Census blocks with all-parks buffer layer
# 7. Union clipped blocks with all blocks, save as new shapefile and add layer
# 8. Arc Toolbox:Data Management:Features:Multipart to Singlepart to separate non-contiguous blocks.
# 9. Delete extraneous fields
# 10. Open data table of clipped/unioned census blocks; add unique identifier
# 11. Add fields: ClipArea, Latitude, Longitude as doubles
# 12. Calculate geometry: ClipArea = Area in square meters
# 13. Calculate geometry: Coordinates of centroids
# 14. Clip unioned layer with desired buffer, save as new shapefile; repeat for all different park types



# R Procedure
# 0. Uses lookup_table_2017_transit.csv for blockgroup subpopulations
# 1. Read in shapefile with all parks stops: parks_all_blocks_mp.shp
# 2. Calculate PctArea from the two ALANDs and ClipArea
# 3. Merge block groups
# 4. Sum Pop10 over block group
# 5. Divide Pop10 by bgPop10 to get percent of population in each block group
# 6. Compute BufferPop17 etc. by (Pop10/bgPop10)*(Pop17)*PctArea
## 7. Compute ExcludePop17 ect. by (Pop10/bgPop10)*(Pop17)*(1-PctArea)
## 8. Run bap subroutine - different commands for all parks compared to other types
# 9. Drop all rows where PctArea == 0
# 10. Compute percent of population within 1/4 mile of any park
# 11. Repeat steps c(1:9)-[6] for other types of parks.

# power bi tips
# No relationships between bap_parks and pop_type_table
# FilterMeasure = 
# VAR ranks = max(pop_type_table[pop_rank])
# RETURN(IF(COUNTROWS(FILTER(bap_parks,bap_parks[bap_rank]<=ranks))>0,1,0))

setwd("L:/Health/Epidemiology/Banwarth_Epi/GIS/GIS_files/benton_parks")

library(maptools)
library(rgdal)
library(data.table)
library(reshape2)
library(plyr)
library(sf)

# read in common datasets

# read in most recent ACS blockgroup populations
bg_parks <- read.csv('blockgroup_benton_parks_2017.csv',sep=',',header=TRUE)
bg_parks$COUNTYFP10 <- as.numeric(substr(as.character(bg_parks$blockgroup),start=3,stop=5))

# all blockgroups population in 2010
bg_pop10 <- read.csv('benton_bg_pop10.csv',header=TRUE,sep=',')
names(bg_pop10) <- c('blockgroup','bgPop10')

# # city blocks
# city_blocks <- read.csv('city_blocks.csv',header=TRUE,sep=",")
# city_blocks$urban <- 1
# city_blocks$blockgroup <- round(as.numeric(as.character(city_blocks$GEOID10))*1e-3,0)
# city_blockgroups <- as.data.table(city_blocks[,c('POP10','blockgroup')])[,lapply(.SD,sum),by=.(blockgroup)]
# names(city_blockgroups) <- c('blockgroup','urbanPop10')

# data set with variable names and ranks
variable_df <- read.csv('variable_list.csv',sep=',',header=TRUE)

# Extracting lists of names
bg_colnames <- names(bg_parks)[c(1:(ncol(bg_parks)-2))]
block_colnames <- paste('block',bg_colnames,sep="")
buffer_colnames <- paste('Buffer',bg_colnames,sep="")



# identifying city blockgroups
# bg_parks <- merge(bg_parks,city_blockgroups,by='blockgroup',all.x=TRUE)
# bg_parks[is.na(bg_parks$urbanPop10),'urbanPop10'] <- 0
# bg_parks <- merge(bg_parks,bg_pop10,by='blockgroup')
# 
# # calculating percent of population in urban part of blockgroups
# bg_urban_colnames <- paste('urban',bg_colnames,sep="")
# bg_parks[,bg_urban_colnames] <- round(bg_parks[,bg_colnames]*bg_parks$urbanPop10/bg_parks$bgPop10,0)

# total populations across regions for different population group types
total_pops <- colSums(bg_parks[,bg_colnames])
pop_pcts <- total_pops/total_pops[1]
variable_df$pop_pct <- pop_pcts
variable_df$pop_rank <- rank(pop_pcts)
variable_df <- variable_df[order(variable_df$pop_rank,decreasing=TRUE),]
write.table(variable_df,'variable_table.txt',row.names=FALSE,sep=',')

# Define coordinate system object
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# number of rows in bap table
baprows <- 10000

# functions

# parks_readin
parks_readin <- function(shapefn,parks_type,typerank){
  
  parks_shp <- st_read(shapefn)
  parks_shp <- st_transform(parks_shp,crswgs84)
  parks_shp <- as(parks_shp,'Spatial')

  parks_df <- data.frame(parks_shp)
  parks_df$FID <- as.numeric(as.character(parks_df$UniqueID))
  
  if(typerank == 1){
    parks_df$FID_parks_ <- as.numeric(as.character(parks_df$FID_parks_))
    parks_df$ALAND10 <- parks_df$Area_1
    parks_df$PctArea <- parks_df$ClipArea/parks_df$ALAND10
    parks_df$PctArea[which(parks_df$PctArea==Inf)] <- 0
    parks_df$GEOID10 <- parks_df$GEOID_1
    parks_df$POP10 <- parks_df$Populati_1
    parks_df$COUNTYFP10 <- parks_df$County_1
    parks_df[,parks_type] <- 0
    parks_df[parks_df$FID_parks_!=-1,parks_type] <- 1
    parks_df$polyID <- as.numeric(rownames(parks_df))
    parks_df <- parks_df[,c('FID','polyID','FID_parks_','GEOID10','ALAND10','POP10','ClipArea','Latitude','Longitude','PctArea',parks_type)]
    parks_shp <<- parks_shp
    return(parks_df)
  } else {
    parks_df <- data.frame(parks_df[,'FID'],rep(1,nrow(parks_df)))
    names(parks_df) <- c('FID',parks_type)
    parks_df <- merge(pcomb,parks_df,by='FID',all.x=TRUE)
    parks_df[is.na(parks_df)] <- 0
    parks_df[parks_df$AllParks==0,parks_type] <- 0
    return(parks_df)
  }
}


# Function to add/subtract population to buffers.
addpop_function <- function(x,tdf,bgdf) {
  x1 <- pcomb[,x]
  x2 <- bg_parks[,sub('Buffer','',x)]
  tadd <- sum(x2)-sum(x1)
  xsamp <- sample(length(x1[x1>0]),abs(tadd))
  if(tadd > 0) {
    x1[x1>0][xsamp] <- x1[x1>0][xsamp]+1
  } else if(tadd < 0) {
    x1[x1>0][xsamp] <- x1[x1>0][xsamp]-1
  }
  return(x1)
}

# function to make one row for each person in bap_exp
bap_fun <- function(x) {
  nx <- pcomb[pcomb$polyID==x,'BufferPop17']
  # print(c(x,nx))
  if(nx>0){
    rand_points <- data.frame(spsample(parks_shp@polygons[[x]]@Polygons[[1]],n=nx,'random'),x)
    return(rand_points)
  } else {
    return(NULL)
  }
}

# function for flagging which population type
poptype_fun <- function(sbpname) {
  cbuffer <- paste('Buffer',sbpname,sep='')
  bap_exp[,sbpname] <- 0
  indices <- unlist(lapply(unique(bap_exp$polyID),function(x) sample(which(bap_exp$polyID==x),min(bap_exp[bap_exp$polyID==x,cbuffer],bap_exp[bap_exp$polyID==x,'BufferPop17']))))
  bap_exp[indices,sbpname] <- 1
  return(bap_exp)
}


### Step 1: Read in parks shapefiles and combine

parklist <- c('AllParks','Natural','Landscaped','Hiking','Playground','SportsFields')

# parks_readin appends each pcomb to make a total pcomb=parks_combined
pcomb <- parks_readin('parks_blocks_mp.shp','AllParks',1)
pcomb <- parks_readin('parks_blocks_natural.shp','Natural',2)
pcomb <- parks_readin('parks_blocks_landscaping.shp','Landscaped',2)
pcomb <- parks_readin('parks_blocks_hiking.shp','Hiking',2)
pcomb <- parks_readin('parks_blocks_playground.shp','Playground',2)
pcomb <- parks_readin('parks_blocks_sportsfields.shp','SportsFields',2)


# Step 2: calculate block and buffer populations

# identifying blockgroups
pcomb$blockgroup <- round(as.numeric(as.character(pcomb$GEOID10))*1e-3,0)

# merge data sets
pcomb <- merge(pcomb,bg_parks,by='blockgroup')
pcomb <- merge(pcomb,bg_pop10,by='blockgroup')

# divide Pop10 by blockgroup pop: blockPop10Pct to get each block's share of the blockgroup population
pcomb$blockPop10Pct <- pcomb$POP10/pcomb$bgPop10
pcomb$blockPop10Pct[is.na(pcomb$blockPop10Pct)] <- 0

# multiply Pop17 etc. by blockPop10Pct= blockPop17 etc. to get each block's Pop17 etc.
pcomb[,block_colnames] <- pcomb[,bg_colnames]*pcomb$blockPop10Pct
pcomb[which(pcomb$bgPop10==0),block_colnames] <- 0

# multiply Pop17block etc. by PctArea = BufferPop17 etc.
pcomb[,buffer_colnames] <- round(pcomb[,block_colnames]*pcomb$PctArea,0)

# randomly assigning extra people to get to total population
pcomb[,buffer_colnames] <- as.data.frame(lapply(buffer_colnames,function(x) addpop_function(x)))


# Step 3: Compute blocks as points

# one row for every person
bap_coords <- data.frame(do.call('rbind',lapply(pcomb$polyID,function(x) bap_fun(x))))
names(bap_coords) <- c('xCoord','yCoord','polyID')

# merge to main table
bap_exp <- merge(bap_coords,pcomb,by='polyID',all.x=TRUE)

# keeping only necessary columns
bap_exp <- bap_exp[,c('polyID','xCoord','yCoord','FID_parks_','COUNTYFP10',parklist,buffer_colnames)]
# Each row is a person
bap_exp$Pop17 <- 1

# indicating population type
for(i in 2:length(bg_colnames)) {
  bap_exp <- poptype_fun(bg_colnames[i])
}

# reordering/dropping columns
bap_exp <- bap_exp[,c('polyID','xCoord','yCoord','COUNTYFP10',parklist,'Pop17',bg_colnames),]
bap_exp$geography <- mapvalues(bap_exp$COUNTYFP10,from=c(3,41,43),to=c('Benton County','Lincoln County','Linn County'))

# generating bap_rank column
bap_exp$bap_rank <- 0
for(i in 1:length(bg_colnames)) {
  # getting names from variable_df since that table is ordered by ranks
  bapcolname <- sub('Buffer','',variable_df[i,'variable'])
  bap_exp[bap_exp[,bapcolname]==1,'bap_rank'] <- length(bg_colnames)+1-i
}

# sample to reduce size
bap_exp <-  bap_exp[sample(nrow(bap_exp),baprows),]

# replicating bap_exp; one for each park type
bap_exp <- bap_exp[rep(row.names(bap_exp),length(parklist)),]
bap_exp$Access <- 0
bap_exp$'Park type' <- rep(c('All parks','Natural areas','Landscaped areas',
                              'Parks with hiking trails','Parks with playgrounds','Parks with sports fields'),each=baprows)
bap_exp$park_rank <- rep(c(1,2,3,4,5,6),each=baprows)
for(i in 1:length(parklist)) {
  j <- (i-1)*baprows+1
  k <- i*baprows
  bap_exp$Access[c(j:k)] <- bap_exp[c(j:k),parklist[i]]
}

# keeping only useful columns
bap_exp <- bap_exp[,c('geography','xCoord','yCoord','bap_rank','Access','Park type','park_rank')]

write.table(variable_df,'pop_type_table.txt',row.names=FALSE,sep=',')
write.table(bap_exp,'bap_parks.txt',row.names = FALSE,sep=',')

# removing tables now that the processing is done
rm(bap_coords)
rm(bap_exp)
rm(bapcolname)
rm(parks_shp)


# Step 4: compute percent of population w/in 1/4 mile: BufferPop17Pct

# need to use bg_parks in denominator because pcomb only has blocks that overlap the buffer

# Dropping all blocks that are outside of the buffer to match the "All parks" data frame with the other parks data frames.
pcomb <- pcomb[which(pcomb$AllParks == 1),]

bufferpct_function <- function(parksCol,parksName,parksRank) {
  # keeping only desired parks_type
  tdf <- pcomb[pcomb[,parksCol]==1,]
  
  # Benton County
  BufferPops <- colSums(tdf[,buffer_colnames])
  TotalPops <- colSums(bg_parks[,bg_colnames])
  BufferPcts <- BufferPops/TotalPops
  
  # reshaping the datatable
  df <- data.frame(cbind(BufferPops,TotalPops,BufferPcts))
  colnames(df) <- c('BufferPop','TotalPop','BufferPct')
  df$variable <- rownames(df)
  df$geography <- as.factor('Benton County')
  
  df$type <- parksName
  df$typerank <- parksRank
  
  df <- merge(df,variable_df,by='variable')
  
  return(df)
}

alldf <- bufferpct_function('AllParks','All parks',1)
localdf <- bufferpct_function('Natural','Natural areas',2)
intercitydf <- bufferpct_function('Landscaped','Landscaped areas',3)
hourlydf <- bufferpct_function('Hiking','Parks with hiking trails',4)
intermittentdf <- bufferpct_function('Playground','Parks with playgrounds',5)
infrequentdf <- bufferpct_function('SportsFields','Parks with sports fields',6)

df <- data.frame(rbind(alldf,localdf,intercitydf,hourlydf,intermittentdf,infrequentdf))
write.csv(df,'parks_access.csv',row.names=FALSE)



