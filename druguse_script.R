# druguse_script.R
# UCI repository: https://archive.ics.uci.edu/ml/datasets/Drug+consumption+%28quantified%29
# data file: https://archive.ics.uci.edu/ml/machine-learning-databases/00373/drug_consumption.data

# Peter Banwarth
# 2/27/2019

# this script uses parameter tuning in the generalized boosted model algorithm to predict drug use


##################
##### IMPORTANT: YOU DO NOT NEED TO RUN THIS SCRIPT BEFORE THE PRESENTATION. WE WILL BE RUNNING THE SCRIPT TOGETHER
##### AT THE PRESENTATION IN ORDER TO DEMONSTRATE MACHINE LEARNING IN ACTION.
##################


### Script Block 1 ######################################################################################################

# install packages
install.packages('data.table')
install.packages('caret')
install.packages'(pROC')
install.packages('plyr')
install.packages('MLmetrics')
install.packages('gbm')

######### EXTRA PACKAGES THAT MAY NEED TO BE INSTALLED ########
install.packages('curl')
install.packages('e1701')

# load packages
library(data.table)
library(caret)
library(pROC)
library(plyr)
library(MLmetrics)
library(gbm)

### Script Block 2 ######################################################################################################

# read in data file
du_website <- 'https://archive.ics.uci.edu/ml/machine-learning-databases/00373/drug_consumption.data'
du_data <- data.frame(fread(du_website,header=F,colClasses = c("numeric",rep("factor",5),rep("numeric",7),rep("character",19))))

## alternate script if you previously saved the data
# setwd('folder path')
# du_data <- read.table('filename.txt',sep=',',header=F,colClasses = c("numeric",rep("factor",5),rep("numeric",7),rep("character",19)))

### Script Block 3 ######################################################################################################

#### Data cleanup and pre-processing ####

# Add column headers
names(du_data) <- c("ID","Age","Gender","Education","Country",
                    "Ethnicity","Neuroticism","Extraversion","Openness",
                    "Agreeableness","Conscientiousness","Impulsiveness",
                    "Sensation","Alcohol","Amphetamine","Amylnitrite",
                    "Benzos","Caffeine","Cannabis","Chocolate","Cocaine",
                    "Crack","Ecstasy","Heroin","Ketamine","Legal_high"
                    ,"LSD","Methadone","Mushrooms","Nicotine","Semeron","VSA")

# relabel drug use categories. Classes:
# Never used or last used more than 1 decade ago: Nonuser
# Other response: current user
du_data[,c(14:ncol(du_data))] <- data.frame(lapply(du_data[,c(14:ncol(du_data))],function(x) {as.numeric(as.character(gsub("CL","",x)))}))
du_data[,c(14:ncol(du_data))] <- floor(du_data[,c(14:ncol(du_data))]^(1/4))

# relabel factor variables
du_data$Age <-  mapvalues(du_data$Age, from = levels(du_data$Age),
                          to = c('Age25_34','Age18_24','Age35_44','Age45_54','Age55_64','Age65+'))
du_data$Gender <-  mapvalues(du_data$Gender, from = levels(du_data$Gender),
                          to = c('Male','Female'))
du_data$Education <-  mapvalues(du_data$Education, from = levels(du_data$Education),
                             to = c('EducProfCert','EducSomeCollege','Educ18','Educ17','Educ16','Educ15','EducBach','EducMasters','EducPhD'))
du_data$Country <-  mapvalues(du_data$Country, from = levels(du_data$Country),
                             to = c('Australia','Other','NewZealand','USA','Ireland','Canada','UK'))
du_data$Ethnicity <-  mapvalues(du_data$Ethnicity, from = levels(du_data$Ethnicity),
                             to = c('EthWhiteBlack','EthWhite','EthAsian','EthBlack','EthOther','EthWhiteAsian','EthBlackAsian'))

### Script Block 4 ######################################################################################################

##### Case 1: Predicting cannabis use #####
# steps:
# partition dataset into train and test sets
# model data one time (i.e. without tuning)
# computer error in model
# model data with parameter tuning
# compute error in best model

# specify predictor variables and outcome variable: Cannabis
predictor_names <- names(du_data)[c(2:13)]
outcome_name <- 'Cannabis'

# set seed for reproducibility
set.seed(123)

### Script Block 5 ######################################################################################################

# partition dataset by randomly selecting 3/4 of the rows
partitionRows <- createDataPartition(y=c(1:nrow(du_data)), p=0.75, list=FALSE)
cannabis_ml <- du_data[partitionRows,]
cannabis_check <- du_data[-partitionRows,]

### Script Block 6 ######################################################################################################

### First model: Generalized Boosted Model without tuning
cannabis_model <- gbm(cannabis_ml[,outcome_name]~.,data=cannabis_ml[,predictor_names],distribution="bernoulli",n.trees=100,verbose=FALSE)

# inspect model
summary(cannabis_model)
print(cannabis_model)

### Script Block 7 ######################################################################################################

# predict on cannabis_check
cannabis_probs <- predict.gbm(cannabis_model,cannabis_check[,predictor_names],n.trees=100,type='response')

# visual inspection of cannabis_probs
hist(cannabis_probs)

### Script Block 8 ######################################################################################################

# recoding cannabis_predictions with most logical split
cannabis_predictions <- cannabis_probs
cannabis_predictions[cannabis_predictions<.76] <- 0
cannabis_predictions[cannabis_predictions>0] <- 1

### Script Block 9 ######################################################################################################

# accuracy of predictions
accuracy <- postResample(pred=as.factor(cannabis_predictions), obs=as.factor(cannabis_check[,outcome_name]))

### Error computation: Log loss computation evaluation
### Use absolute value of logloss
logloss <- abs(LogLoss(cannabis_probs,as.numeric(cannabis_check[,outcome_name])-1))

### Script Block 10 ######################################################################################################

### Second model: Generalized Boosted Model with tuning

# make drug use a factor for caret package
cannabis_ml[,c(14:ncol(cannabis_ml))] <- data.frame(lapply(cannabis_ml[,c(14:ncol(cannabis_ml))],function(x) {as.factor(ifelse(x==1,'yes','no'))}))
cannabis_check[,c(14:ncol(cannabis_check))] <- data.frame(lapply(cannabis_check[,c(14:ncol(cannabis_check))],function(x) {as.factor(ifelse(x==1,'yes','no'))}))

### Script Block 11 ######################################################################################################

# training settings
train_settings <- trainControl(method='cv',number=5,returnResamp = 'none',summaryFunction=twoClassSummary,classProbs=TRUE)
# tuneGrid is default: tgrid <- expand.grid(interaction.depth=c(1:3),n.trees=c(1:3)*50,shrinkage=0.1,n.minobsinnode=10)

### Script Block 12 ######################################################################################################

# run generalized boosted model
cannabis_model1 <- train(cannabis_ml[,predictor_names],cannabis_ml[,outcome_name],
                  method='gbm',
                  trControl=train_settings,
                  metric="ROC",
                  verbose=FALSE)

# inspect model
summary(cannabis_model1)

print(cannabis_model1)

### Script Block 13 ######################################################################################################

# predict on cannabis_check
cannabis_predictions1 <- predict(object=cannabis_model1,cannabis_check[,predictor_names],type='raw')

# accuracy of predictions
accuracy1 <- postResample(pred=cannabis_predictions1, obs=cannabis_check[,outcome_name])

### Error computation: Log loss computation evaluation
cannabis_probs1 <- predict(object=cannabis_model1,cannabis_check[,predictor_names],type='prob')
logloss1 <- LogLoss(cannabis_probs1$yes,as.numeric(cannabis_check[,outcome_name])-1)

### Script Block 14 ######################################################################################################

#### Compare one-time gbm with tuned gbm
print(data.frame(label=c('One-time gbm accuracy','tuned gbm accuracy','one-time gbm log loss','tuned gbm log-loss'),value=c(accuracy[1],accuracy1[1],logloss,logloss1)))

### Script Block 15 ######################################################################################################

#### Case 2: Repeat exercise for much rarer drug class: Heroin ###

table(du_data$Heroin)
outcome_name <- 'Heroin'

# partition dataset by randomly selecting 3/4 of the rows
partitionRows <- createDataPartition(y=c(1:nrow(du_data)), p=0.75, list=FALSE)
heroin_ml <- du_data[partitionRows,]
heroin_check <- du_data[-partitionRows,]

### First model: Generalized Boosted Model without tuning
heroin_model <- gbm(heroin_ml[,outcome_name]~.,data=heroin_ml[,predictor_names],distribution="bernoulli",n.trees=100,verbose=FALSE)

### Script Block 16 ######################################################################################################

# inspect model
# summary(heroin_model)
# print(heroin_model)

# predict on heroin_check
heroin_probs <- predict.gbm(heroin_model,heroin_check[,predictor_names],n.trees=100,type='response')

# visual inspection of heroin_probs
hist(heroin_probs)

# recoding heroin_predictions from visual inspection
heroin_predictions <- heroin_probs
heroin_predictions[heroin_predictions<.15] <- 0
heroin_predictions[heroin_predictions>0] <- 1

# accuracy of predictions
accuracy <- postResample(pred=as.factor(heroin_predictions), obs=as.factor(heroin_check[,outcome_name]))

# Log loss computation evaluation
# Use absolute value
logloss <- abs(LogLoss(heroin_probs,as.numeric(heroin_check[,outcome_name])-1))

### Script Block 17 ######################################################################################################

### Second model Generalized Boosted Model with tuning

# make drug use a factor for caret package
heroin_ml[,c(14:ncol(heroin_ml))] <- data.frame(lapply(heroin_ml[,c(14:ncol(heroin_ml))],function(x) {as.factor(ifelse(x==1,'yes','no'))}))
heroin_check[,c(14:ncol(heroin_check))] <- data.frame(lapply(heroin_check[,c(14:ncol(heroin_check))],function(x) {as.factor(ifelse(x==1,'yes','no'))}))

# training settings
train_settings <- trainControl(method='cv',number=5,returnResamp = 'none',summaryFunction=twoClassSummary,classProbs=TRUE)
# tuneGrid is default: tgrid <- expand.grid(interaction.depth=c(1:3),n.trees=c(1:3)*50,shrinkage=0.1,n.minobsinnode=10)

# run generalized boosted model
heroin_model1 <- train(heroin_ml[,predictor_names],heroin_ml[,outcome_name],
                   method='gbm',
                   trControl=train_settings,
                   metric="ROC",
                   verbose=FALSE)

# inspect model
# summary(heroin_model1)
# print(heroin_model1)

### Script Block 18 ######################################################################################################

# predict on heroin_check
heroin_predictions1 <- predict(object=heroin_model1,heroin_check[,predictor_names],type='raw')

# accuracy of predictions
accuracy1 <- postResample(pred=heroin_predictions1, obs=heroin_check[,outcome_name])

# Log loss computation evaluation
heroin_probs1 <- predict(object=heroin_model1,heroin_check[,predictor_names],type='prob')
logloss1 <- LogLoss(heroin_probs1$yes,as.numeric(heroin_check[,outcome_name])-1)

### Script Block 19 ######################################################################################################

#### Compare one-time gbm with tuned gbm
print(data.frame(label=c('One-time gbm accuracy','tuned gbm accuracy','one-time gbm log loss','tuned gbm log-loss'),value=c(accuracy[1],accuracy1[1],logloss,logloss1)))
