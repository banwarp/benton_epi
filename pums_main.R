# pums_main.R
# Peter Banwarth
# 10/3/2018

# Supervised machine learning using PUMS data to predict census-tract level tables

# load libraries
library(plyr)
library(dplyr)
library(lazyeval)
library(ggplot2)
library(tidyr)
library(reshape2)
library(data.table)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(caret)
library(pROC)
library(raster)
library(tibble)

########################
####### scripts ########
########################

setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/ACS_tables_script/processing_scripts')
source('acs_data_extraction_us.R')

# geography script
geo_function <- function(abbr,y) {
  setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists/us/geo')
  fn <- paste('g',y,'5',abbr,'.csv',sep='')
  geo_df <- read.csv(fn,header=FALSE,sep=',')
  pumas <- geo_df[!is.na(geo_df[,47]),c(5,47,50)]
  rm(geo_df)
  pumas[,3] <- as.character(pumas[,3])
  return(pumas)
}

puma_prep <- function(abbr,y,lookup_table,agg_vector,granularity,year_folder) {
  # generate list of pumas
  pumas <- geo_function(abbr,y)
  names(pumas) <- c('puma_geo','puma_num','puma_names')
  pumas[,'abbr'] <- abbr
  puma_geo <- pumas[,1]
  puma_num <- pumas[,2]
  puma_names <- pumas[,3]
  puma_abbr <- pumas[,4]
  
  # run acs data extraction
  df <- acs_function(abbr,puma_geo,puma_names,lookup_table,agg_vector,granularity,year_folder)
  names(df)[which(names(df)=='geo_names')] <- 'puma_names'
  
  
  # appending variables of interest
  df <- cbind(puma_abbr,df)
  df <- merge(pumas[,c('puma_num','puma_names')],df)
}

# acs data extraction script
acs_function <- function(abbr,geo_vector,geo_names,lookup_table,agg_vector,granularity,year_folder) {
  # extract data
  df <- acs_data_extraction_us(lookup_table, geo_vector, geo_names, agg_vector, granularity,year_folder,abbr)

  # compute percentages
  universe <- data.frame('Table_Number' = unique(df[,'Table_Number']),
                         'Universe' = c(0,rep(1,8),0,rep(2,8),3,3,3,2,3,2,3,rep(2,11)))
  pop_denom <- df[which(df[,'Table_Number']=='B01001'),c(1:(6+length(geo_vector)))]
  housing_denom <- df[which(df[,'Table_Number']=='B11001'),c(1:(6+length(geo_vector)))]
  
  meta_df <- df[,c(1:6)]
  
  df_pop <- df[which(df[,'Table_Number'] %in% c(universe[which(universe['Universe']==1,1)])),c(1:(6+length(geo_vector)))]
  df_housing <- df[which(df[,'Table_Number'] %in% c(universe[which(universe['Universe']==2,1)])),c(1:(6+length(geo_vector)))]
  df_eigen <- df[which(df[,'Table_Number'] %in% c(universe[which(universe['Universe']==3,1)])),c(1:6,(7+length(geo_vector)):(6+2*length(geo_vector)))]
  
  df_pop[,c(7:ncol(df_pop))] <- rbind.fill(apply(df_pop[,c(7:ncol(df_pop))],1,function(x) x/pop_denom[c(7:(6+length(geo_vector)))]))
  df_housing[,c(7:ncol(df_housing))] <- rbind.fill(apply(df_housing[,c(7:ncol(df_housing))],1,function(x) x/housing_denom[c(7:(6+length(geo_vector)))]))
  names(df_eigen) <- c(names(meta_df),geo_names)

  df <- rbind(pop_denom,housing_denom,df_pop,df_housing,df_eigen)
  
  # dropping totals and first category for each set, dropping Datapoint
  total_list <- which(df[,'Datapoint']=='Total:')[-1] # excepting population
  df <- df[-sort(c(total_list,total_list+1)),] # drops the +1 for collinearity
  
  df <- df[order(df[,'Table_line']),-c(1,3:6)]
  rownames(df) <- df$Table_line
  
  # transposing
  df <- as.data.frame(t(df[,-1]))
  df[,'geo_names'] <- rownames(df)
  
  return(df)
}

pop_read_function <- function(alpha,y) {
  pop_file <- paste('ss',y,'pus',alpha,'.csv',sep='')
  pop_df <- fread(pop_file,select=c('SERIALNO','POVPIP','PWGTP'))
  pop_df <- pop_df[which(pop_df[,'POVPIP'] < 185),]
  return(pop_df)
}

# reading the huge housing files
hou_read_function <- function(alpha,y) {
  hou_file <- paste('ss',y,'hus',alpha,'.csv',sep='')
  hou_df <- fread(hou_file,select=c('PUMA','ST','SERIALNO','FS'))
  return(hou_df)
}


##################### End of scripts ##########################

# Geographies list
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists/us/geo')
abbr_df <- read.csv('abbr_df.csv',sep=',',header=TRUE)

########## ACS data ############
# acs data extraction parameters
lookup_table_arg <- 'lookup_2016_pums.csv'
agg_vector_arg <- vector()
granularity_arg <- 'geographies'
year_folder_arg <- "/us2016"

# Data extraction from all geographies

### uncomment to run

acs_list <- lapply(abbr_df[,1],puma_prep,y=2016,
                lookup_table=lookup_table_arg,
                agg_vector=agg_vector_arg,
                granularity=granularity_arg,
                year_folder=year_folder_arg)

# acs data frame
acs_df <- rbind.fill(acs_list)
rm(acs_list)
acs_df <- merge(abbr_df,acs_df)


# getting land area
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/puma_shapefiles")
file_list <- list.files()
file_list <- file_list[grep(".shp$",file_list)]
shape_list <- lapply(file_list,shapefile)
puma_area_list <- lapply(shape_list,function(x) data.frame(cbind(x[,1],x[,2],x[,7])))
puma_area_df <- rbind.fill(puma_area_list)
rm(puma_area_list)
for(i in 1:3) {
 puma_area_df[,i] <- as.numeric(puma_area_df[,i])
}

names(puma_area_df) <- c('ST','puma_num','land_area')

acs_df <- merge(acs_df,puma_area_df,by=c('ST','puma_num'))

# computing population density
acs_df[,'9'] <- acs_df[,'9']/acs_df[,'land_area']
acs_df <- acs_df[,which(names(acs_df)!='land_area')]

# saving data frame for faster loading
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.table(acs_df, file="acs_pums_usa.csv", append = FALSE, sep=",", row.names=FALSE)

# # readin if table already built; comment to hide
# setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
# acs_df <- read.csv('acs_pums_usa.csv',header=T,sep=',')

########## PUMS data ###########

# set working directory
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums')

join_list <- c('a','b','c','d')

# getting unique households under 1.85*FPL
pop_list <- lapply(join_list,pop_read_function,y=16)
pop_df <- rbind.fill(pop_list)
rm(pop_list)

hou_list <- lapply(join_list,hou_read_function,y=16)
hou_df <- rbind.fill(hou_list)
rm(hou_list)

# join population and housing data
puma_df <- merge(pop_df,hou_df,by='SERIALNO')
rm(hou_df)
rm(pop_df)

puma_df[which(puma_df[,'FS']==2),'FS'] <- 0

# Generating weighted SNAP utilization variable
snap <- as.data.table(puma_df[,c('PUMA','ST','FS','PWGTP')])[,lapply(.SD,weighted.mean,w=PWGTP),by=.(PUMA,ST),.SDcols='FS']
# snap <- as.data.table(puma_df[,c('PUMA','ST','FS')])[,lapply(.SD,mean),by=.(PUMA,ST)]
names(snap) <- c('puma_num','ST','SNAP')
rm(puma_df)

# Saving snap
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.table(snap, file="snap_df.csv", append = FALSE, sep=",", row.names=FALSE)

# merging SNAP and ACS data

df <- as.data.frame(merge(snap,acs_df,by=c('puma_num','ST')))
df <- df[,c(1,3,2,4:ncol(df))]

# Saving df
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.table(df, file="snap_acs_df.csv", append = FALSE, sep=",", row.names=FALSE)
df <- df[,-c(4,5)]
# hack to fix the 'X' problem
names(df)[-c(1:3)] <- paste('X',names(df)[-c(1:3)],sep='')

# dummying state variable
df[,'ST'] <- as.factor(df[,'ST'])
STdummy <- dummyVars(~., data = df)
df <- data.frame(predict(STdummy,df))
names(df) <- gsub("X","",names(df))

# scaling and centering
# storing sd and mean for reverse-scaling predictions later
SNAP_sd <- apply(df[,-1],2,sd)
SNAP_mean <- apply(df[,-1],2,mean)
# saving sd and mean
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.csv(data.frame(SNAP_sd,SNAP_mean),file='SNAP_scaling.csv')

df[,-1] <- apply(df[,-1],2,function(x) (x-mean(x))/sd(x))


# Machine learning

library(nnet)
library(caret)
library(pROC)

set.seed(20)

train_df <- df[,-1]
inTrain <- createDataPartition(y=c(1:nrow(train_df)), p=0.75, list=FALSE)   # We wish 75% for the trainset 

train.set <- train_df[inTrain,]
test.set  <- train_df[-inTrain,]
nrow(train.set)/nrow(test.set) # should be around 3

fitControl <- trainControl(method='repeatedcv',
                           number=10,
                           repeats=5, returnResamp='none')

glmnetGrid <- expand.grid(alpha=seq(from=0,to=1,by=0.1),
                        lambda=seq(from=0.0001,to=1,length=100))

model.glmnet <- train(SNAP ~ ., data=train.set, method='glmnet', metric='RMSE',
                    trControl=fitControl,tuneGrid=glmnetGrid)

prediction.glmnet <- predict(model.glmnet, test.set[-1])
plot(prediction.glmnet,test.set$SNAP)

# saving model
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
saveRDS(model.glmnet,'model_glmnet.rds')
# model.glmnet <- readRDS('model_glmnet.rds')

# model on 100% of the data
model.glmnet.100 <- train(SNAP ~ ., data=train_df, method='glmnet', metric='RMSE',
                          trControl=fitControl,tuneGrid=glmnetGrid)
saveRDS(model.glmnet.100,'model_glmnet_100.rds')

#################################################
##### Applying model to Oregon Census Tracts ####
#################################################

# applied model
applied_model <- model.glmnet
# applied_model <- model.glmnet.100

# Reading in census tracts in Oregon
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists')
tract_vector <- scan("OR_all_census_tracts_2016.txt",sep=",")
tract_names <- scan("OR_all_census_tracts_names.txt",what="",sep=",")
agg_vector <- vector()
granularity_ct_bg <- "tracts_blocks"
year_folder <- "/year2016"
abbr <- 'or'

# building data set
or_df <- acs_function(abbr,tract_vector,tract_names,lookup_table_arg,agg_vector,granularity_ct_bg,year_folder)
or_df <- or_df[,c(ncol(or_df),1:(ncol(or_df)-1))]

# getting land area
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists/tract_shape_files")
or_land <- data.frame(shapefile('cb_2017_41_tract_500k.shp'))
or_land <- or_land[,c(5,8)]
names(or_land) <- c('geo_names','land_area')
or_df <- merge(or_df,or_land, by='geo_names')
or_df[,'land_area'] <- as.numeric(or_df[,'land_area'])
or_df[,'9'] <- or_df[,'9']/or_df[,'land_area']
or_df <- or_df[,which(names(or_df)!='land_area')]

# Dummying state variable - all Oregon
ST_df <- df[c(1:nrow(or_df)),c(which(names(df)=='ST.1'):which(names(df)=='ST.56'))]
ST_df[] <- 0
ST_df['ST.41'] <- 1

or_df <- data.frame(or_df[,1],ST_df,or_df[,-1])
names(or_df)[1] <- 'geo_names'
names(or_df) <- gsub("X","",names(or_df))

# scaling and centering
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
SNAP_scaling <- read.csv('SNAP_scaling.csv', sep=',',row.names=1)
or_df[,-1] <- data.frame(t(apply(or_df[,-1],1,function(x) (x-SNAP_scaling[-1,2])/SNAP_scaling[-1,1])))

# Drop the 9 ocean census tracts (or whichever has NA)
or_df[,'drop'] <- rowSums(or_df[,-1]) # easy way to find out the rows with NA
or_df <- or_df[-which(is.na(or_df[,'drop'])),-which(names(or_df)=='drop')]

or_prediction <- predict(applied_model,or_df[,-1])

or_prediction_rescaled <- or_prediction*SNAP_scaling[1,1]+SNAP_scaling[1,2]

# write predictions to file
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.table(data.frame(geography=as.numeric(or_df[,'geo_names']),SNAP_utilization=or_prediction_rescaled),
            'SNAP_predictions.csv',row.names=FALSE)



#################################################
##### Applying model to RHA Counties ####
#################################################

# Reading in census tracts in Oregon
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists')
county_vector <- c(14,33,34)
county_names <- c('Benton','Lincoln','Linn')
agg_vector <- vector()
granularity_ct_bg <- "geographies"
year_folder <- "/year2016"
abbr <- 'or'

# building data set
or_df <- acs_function(abbr,county_vector,county_names,lookup_table_arg,agg_vector,granularity_ct_bg,year_folder)
or_df <- or_df[,c(ncol(or_df),1:(ncol(or_df)-1))]

# getting land area
setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/geography_lists/county_shape_files")
or_land <- data.frame(shapefile('cb_2017_us_county_500k.shp'))

or_land <- or_land[which(or_land$STATEFP==41),]
or_land <- or_land[which(or_land$NAME %in% c('Benton','Lincoln','Linn')),c(6,8)]

names(or_land) <- c('geo_names','land_area')
or_df <- merge(or_df,or_land, by='geo_names')
or_df[,'land_area'] <- as.numeric(or_df[,'land_area'])
or_df[,'9'] <- or_df[,'9']/or_df[,'land_area']
or_df <- or_df[,which(names(or_df)!='land_area')]

# Dummying state variable - all Oregon
ST_df <- df[c(1:nrow(or_df)),c(which(names(df)=='ST.1'):which(names(df)=='ST.56'))]
ST_df[] <- 0
ST_df['ST.41'] <- 1

or_df <- data.frame(or_df[,1],ST_df,or_df[,-1])
names(or_df)[1] <- 'geo_names'
names(or_df) <- gsub("X","",names(or_df))

# scaling and centering
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
SNAP_scaling <- read.csv('SNAP_scaling.csv', sep=',',row.names=1)
or_df[,-1] <- data.frame(t(apply(or_df[,-1],1,function(x) (x-SNAP_scaling[-1,2])/SNAP_scaling[-1,1])))

# Drop the 9 ocean census tracts (or whichever has NA)
or_df[,'drop'] <- rowSums(or_df[,-1]) # easy way to find out the rows with NA
or_df <- or_df[-which(is.na(or_df[,'drop'])),-which(names(or_df)=='drop')]

or_prediction <- predict(model.glmnet,or_df[,-1])

or_prediction_rescaled <- or_prediction*SNAP_scaling[1,1]+SNAP_scaling[1,2]

# write predictions to file
setwd('L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/pums/output')
write.table(data.frame(geography=as.numeric(or_df[,'geo_names']),SNAP_utilization=or_prediction_rescaled),
            'SNAP_predictions.csv',row.names=FALSE)


# DONE
