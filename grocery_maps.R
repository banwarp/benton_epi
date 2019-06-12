# grocery_maps.R
# Access to public grocery


# ArcMap Procedure
# 
# 1. Select desired grocery stops and create new layer from selection.
# 2. Export selection as shape file.
# 3. Geoprocessing - buffer: Choose layer, name buffer, set distance, Dissolve = All.
# 4. Add Census blocks to map.
# 5. Clip Census blocks with all-grocery buffer layer
# 6. Union clipped blocks with all blocks, save as new shapefile and add layer
# 7. Arc Toolbox:Data Management:Features:Multipart to Singlepart to separate non-contiguous blocks.
# 8. Delete extraneous fields
# 9. Open data table of clipped/unioned census blocks
# 10. Add fields: ClipArea
# 11. Calculate geometry: ClipArea = Area in square meters
# 12. Clip unioned layer with desired buffer, save as new shapefile; repeat for all different route types



# R Procedure
# 1. Read in shapefile with all grocery stops
# 2. Calculate PctArea from the two ALANDs and ClipArea
# 3. Merge block groups
# 4. Sum Pop10 over block group
# 5. Divide Pop10 by bgPop10 to get percent of population in each block group
# 6. Compute BufferPop17 etc. by (Pop10/bgPop10)*(Pop17)*PctArea
## 7. Compute ExcludePop17 ect. by (Pop10/bgPop10)*(Pop17)*(1-PctArea)
## 8. Run bap subroutine - different commands for all grocery stops compared to other types
# 9. Drop all rows where PctArea == 0
# 10. Compute percent of population within 1/4 mile of any grocery stop
# 11. Repeat steps c(1:9)-[6] for other types of grocery stops.

# power bi tips
# No relationships between bap_grocery and pop_type_table
# FilterMeasure = 
# VAR ranks = max(pop_type_table[pop_rank])
# RETURN(IF(COUNTROWS(FILTER(bap_grocery,bap_grocery[bap_rank]<=ranks))>0,1,0))

setwd("L:/Health/Epidemiology/Banwarth_Epi/GIS/GIS_files/lbl_grocery")

library(maptools)
library(rgdal)
library(data.table)
library(reshape2)
library(plyr)
library(sf)

# read in common datasets

# read in most recent ACS blockgroup populations
bg_grocery <- read.table('blockgroup_lbl_grocery_2017.txt',sep=',',header=TRUE)
bg_grocery <- bg_grocery[,c('Pop17','Children','Seniors','White','Hispanic','Nonwhite','Poverty','SNAP','Disability','Renter','HousingCost','OHP','blockgroup')]
bg_grocery$COUNTYFP10 <- as.numeric(substr(as.character(bg_grocery$blockgroup),start=3,stop=5))

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
bg_colnames <- names(bg_grocery)[-c(13,14)]
block_colnames <- paste('block',bg_colnames,sep="")
buffer_colnames <- paste('Buffer',bg_colnames,sep="")



# identifying city blockgroups
bg_grocery <- merge(bg_grocery,city_blockgroups,by='blockgroup',all.x=TRUE)
bg_grocery[is.na(bg_grocery$urbanPop10),'urbanPop10'] <- 0
bg_grocery <- merge(bg_grocery,bg_pop10,by='blockgroup')

# calculating percent of population in urban part of blockgroups
bg_urban_colnames <- paste('urban',bg_colnames,sep="")
bg_grocery[,bg_urban_colnames] <- round(bg_grocery[,bg_colnames]*bg_grocery$urbanPop10/bg_grocery$bgPop10,0)

# total populations across regions for different population group types
total_pops <- colSums(bg_grocery[,bg_colnames])
pop_pcts <- total_pops/total_pops[1]
variable_df$pop_pct <- pop_pcts
variable_df$pop_rank <- rank(pop_pcts)
variable_df <- variable_df[order(variable_df$pop_rank,decreasing=TRUE),]
write.table(variable_df,'variable_table.txt',row.names=FALSE,sep=',')

# Define coordinate system object
crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# number of rows in bap table
baprows <- 20000

# grocery_readin
grocery_readin <- function(shapefn,grocery_type,typerank){
  
  grocery_shp <- st_read(shapefn)
  grocery_shp <- st_transform(grocery_shp,crswgs84)
  grocery_shp <- as(grocery_shp,'Spatial')
  
  grocery_df <- data.frame(grocery_shp)
  grocery_df$FID <- grocery_df$UniqueID
  
  if(typerank == 1){
    grocery_df$ALAND10 <- grocery_df$Area_1
    grocery_df$PctArea <- grocery_df$ClipArea/grocery_df$ALAND10
    grocery_df$PctArea[which(grocery_df$PctArea==Inf)] <- 0
    grocery_df$GEOID10 <- pmax(as.numeric(as.character(grocery_df$GEOID10)),as.numeric(as.character(grocery_df$GEOID10_1)),na.rm = TRUE)
    grocery_df$POP10 <- pmax(as.numeric(as.character(grocery_df$POP10)),as.numeric(as.character(grocery_df$POP10_1)),na.rm=TRUE)
    grocery_df$COUNTYFP10 <- pmax(as.numeric(as.character(grocery_df$COUNTYFP10)),as.numeric(as.character(grocery_df$COUNTYFP_1)),na.rm=TRUE)
    grocery_df[,grocery_type] <- 0
    grocery_df[grocery_df$FID_grocer!=-1,grocery_type] <- 1
    grocery_df$polyID <- as.numeric(rownames(grocery_df))
    grocery_df <- grocery_df[,c('FID','FID_grocer','polyID','GEOID10','ALAND10','POP10','ClipArea','PctArea',grocery_type)]
    grocery_shp <<- grocery_shp
    return(grocery_df)
  } else {
    grocery_df <- data.frame(grocery_df[,'FID'],rep(1,nrow(grocery_df)))
    names(grocery_df) <- c('FID',grocery_type)
    grocery_df <- merge(gcomb,grocery_df,by='FID',all.x=TRUE)
    grocery_df[is.na(grocery_df)] <- 0
    grocery_df[grocery_df$AllRoutes==0,grocery_type] <- 0
    return(grocery_df)
  }
}


# Function to add/subtract population to buffers.
addpop_function <- function(x,tdf,bgdf) {
  x1 <- gcomb[,x]
  x2 <- bg_grocery[,sub('Buffer','',x)]
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
  nx <- gcomb[gcomb$polyID==x,'BufferPop17']
  # print(c(x,nx))
  if(nx>0){
    rand_points <- data.frame(spsample(grocery_shp@polygons[[x]]@Polygons[[1]],n=nx,'random'),x)
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



### Step 1: Read in grocery shapefiles and combine

grocerylist <- c('AllStores','Groceries','WIC','Pantries','Convenience','Tobacco')

# grocery_readin appends each gcomb to make a total gcomb=grocery_combined
gcomb <- grocery_readin('grocery_blocks_mp.shp','AllStores',1)
gcomb <- grocery_readin('grocery_blocks_grocery.shp','Groceries',2)
gcomb <- grocery_readin('grocery_blocks_WIC.shp','WIC',2)
gcomb <- grocery_readin('grocery_blocks_pantry.shp','Pantries',2)
gcomb <- grocery_readin('grocery_blocks_convenience.shp','Convenience',2)
gcomb <- grocery_readin('grocery_blocks_tobacco.shp','Tobacco',2)

# Step 2: calculate block and buffer populations

# identifying blockgroups
gcomb$blockgroup <- round(as.numeric(as.character(gcomb$GEOID10))*1e-3,0)

# merge data sets
gcomb <- merge(gcomb,bg_grocery,by='blockgroup')
# gcomb <- merge(gcomb,bg_pop10,by='blockgroup')

# divide Pop10 by blockgroup pop: blockPop10Pct to get each block's share of the blockgroup population
gcomb$blockPop10Pct <- gcomb$POP10/gcomb$bgPop10
gcomb$blockPop10Pct[is.na(gcomb$blockPop10Pct)] <- 0

# multiply Pop17 etc. by blockPop10Pct= blockPop17 etc. to get each block's Pop17 etc.
gcomb[,block_colnames] <- gcomb[,bg_colnames]*gcomb$blockPop10Pct
gcomb[which(gcomb$bgPop10==0),block_colnames] <- 0

# multiply Pop17block etc. by PctArea = BufferPop17 etc.
gcomb[,buffer_colnames] <- round(gcomb[,block_colnames]*gcomb$PctArea,0)

# randomly assigning extra people to get to total population
gcomb[,buffer_colnames] <- as.data.frame(lapply(buffer_colnames,function(x) addpop_function(x)))



# Step 3: Compute blocks as points

# one row for every person
bap_coords <- data.frame(do.call('rbind',lapply(gcomb$polyID,function(x) bap_fun(x))))
names(bap_coords) <- c('xCoord','yCoord','polyID')

# merge to main table
bap_exp <- merge(bap_coords,gcomb,by='polyID',all.x=TRUE)

# keeping only necessary columns
bap_exp <- bap_exp[,c('polyID','xCoord','yCoord','FID_grocer','COUNTYFP10',grocerylist,buffer_colnames)]
# Each row is a person
bap_exp$Pop17 <- 1

# indicating population type
for(i in 2:length(bg_colnames)) {
  bap_exp <- poptype_fun(bg_colnames[i])
}

# reordering/dropping columns
bap_exp <- bap_exp[,c('polyID','xCoord','yCoord','COUNTYFP10',grocerylist,bg_colnames),]
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

# replicating bap_exp; one for each grocery type
bap_exp <- bap_exp[rep(row.names(bap_exp),length(grocerylist)),]
bap_exp$Access <- 0

bap_exp$'Store type' <- rep(c('All stores','Large or small groceries','WIC authorized stores',
                              'Food pantries','Convenience stores','Tobacco retailers'),each=baprows)

bap_exp$grocery_rank <- rep(c(1,2,3,4,5,6),each=baprows)
for(i in 1:length(grocerylist)) {
  j <- (i-1)*baprows+1
  k <- i*baprows
  bap_exp$Access[c(j:k)] <- bap_exp[c(j:k),grocerylist[i]]
}

# keeping only useful columns
bap_exp <- bap_exp[,c('geography','xCoord','yCoord','bap_rank','Access','Store type','grocery_rank')]

write.table(variable_df,'pop_type_table.txt',row.names=FALSE,sep=',')
write.table(bap_exp,'bap_grocery.txt',row.names = FALSE,sep=',')

# removing tables now that the processing is done
rm(bap_coords)
rm(bap_exp)
rm(bapcolname)
rm(grocery_shp)


# Step 4: compute percent of population w/in 1/4 mile: BufferPop17Pct

# need to use bg_grocery in denominator because gcomb only has blocks that overlap the buffer

# Dropping all blocks that are outside of the buffer to match the "All grocery" data frame with the other routes data frames.
gcomb <- gcomb[which(gcomb$AllStores == 1),]

bufferpct_function <- function(groceryCol,groceryName,groceryRank) {
  # keeping only desired grocery_type
  tdf <- gcomb[gcomb[,groceryCol]==1,]
  
  # three county region
  regionBufferPops <- colSums(tdf[,buffer_colnames])
  regionTotalPops <- colSums(bg_grocery[,bg_colnames])
  regionBufferPcts <- regionBufferPops/regionTotalPops
  
  # each county
  countyBufferPops <- as.data.table(tdf[,c('COUNTYFP10',buffer_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
  countyTotalPops <- as.data.table(bg_grocery[,c('COUNTYFP10',bg_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
  countyTotalPops <- countyTotalPops[which(countyTotalPops$COUNTYFP10 %in% unique(countyBufferPops$COUNTYFP10)),]
  countyBufferPcts <- countyBufferPops[,-1]/countyTotalPops[,-1]
  
  # three county urban/rural
  # city blocks
  tdf <- merge(tdf,city_blocks[,c('GEOID10','urban')],by='GEOID10',all.x=TRUE)
  tdf[is.na(tdf$urban),'urban'] <- 0
  
  # whole region urban/rural
  regionUrbanBufferPops <- as.data.table(tdf[,c('urban',buffer_colnames)])[,lapply(.SD,sum),by=.(urban)]
  regionUrbanTotalPops <- colSums(bg_grocery[,bg_urban_colnames])
  regionUrbanTotalPops <- rbind(regionUrbanTotalPops,regionTotalPops-regionUrbanTotalPops)
  regionUrbanBufferPcts <- regionUrbanBufferPops[,-1]/regionUrbanTotalPops
  regionUrbanTotalPops <- data.frame(regionUrbanTotalPops)
  
  # each county urban/rural
  countyUrbanBufferPops <- as.data.table(tdf[,c('COUNTYFP10','urban',buffer_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10,urban)]
  countyUrbanBufferPops <- countyUrbanBufferPops[order(countyUrbanBufferPops$COUNTYFP10,decreasing=FALSE),]
  countyUrbanBufferPops <- countyUrbanBufferPops[order(countyUrbanBufferPops$urban,decreasing=TRUE),]
  countyUrbanTotalPops <- as.data.table(bg_grocery[,c('COUNTYFP10',bg_urban_colnames)])[,lapply(.SD,sum),by=.(COUNTYFP10)]
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
  df$type <- groceryName
  df$typerank <- groceryRank
  
  df <- merge(df,variable_df,by='variable')
  
  return(df)
}

alldf <- bufferpct_function('AllStores','All stores',1)
localdf <- bufferpct_function('Groceries','Large or small groceries',2)
intercitydf <- bufferpct_function('WIC','WIC authorized stores',3)
hourlydf <- bufferpct_function('Pantries','Food pantries',4)
intermittentdf <- bufferpct_function('Convenience','Convenience stores',5)
infrequentdf <- bufferpct_function('Tobacco','Tobacco retailers',6)

df <- data.frame(rbind(alldf,localdf,intercitydf,hourlydf,intermittentdf,infrequentdf))
write.csv(df,'grocery_access.csv',row.names=FALSE)

