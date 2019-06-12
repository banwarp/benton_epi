# acs_data_extraction
# Peter Banwarth
# 12/13/2016

# This function reads in a slightly modified sequence_table_number_lookup_5yr file from the ACS and uses it
# to read and aggregate specific data lines from the ACS ftp sequence data tables

# rm(list = ls())

# steps:
# read in lookup_table
# identify relevant data lines, including universe totals
# initalize empty data frame for storing results
# populate data frame with data for discrete geographies
# aggregate non-flagged data lines to regions
# compute percentages for non-flagged data lines
# compute margins of error for aggregated and percentage data
# compute relative standard error for reliability
# set flagged aggregates and percentages to -99
# return data frame

# arguments:
# lookup_table_name.  A file name in quotes with .csv of the table used to determine which rows are processed
# geo_vector.  A vector of geography codes used to determine which geographies are included and agreggated
# geo_names.  A vector of geography names used for column headings, including the aggregated regions
# agg_vector.  A vector identifying which of the geographies to aggregate by listing the start and end
#              positions for each set of geographies
# granularity.  geographies or tracts_blocks
# year_folder.  Which year to get the data from

acs_data_extraction <- function(lookup_table_name, geo_vector, geo_names, agg_vector,granularity,year_folder) {

  # Subroutines used in this function
  # Euclidean norm for vectors
  euclidnorm <- function(x){
    return(sqrt(sum(x^2)))
  }
  
  # Percent margin of error function
  pct_moe_function <- function(input_vector){
    a <- input_vector[1]
    b <- input_vector[2]
    sa <- input_vector[3]
    sb <- input_vector[4]
    if(is.na(max(input_vector))) {
      return(max(input_vector))
    } else if(b==0) {
      return(NA)
    } else if(a/b == 1) {
      return(sa/b)
    } else if(sa^2-(((a/b)^2))*sb^2 < 0) {
      return((1/b)*sqrt(sa^2+((a/b)^2*sb^2)))
    } else {
      return((1/b)*sqrt(sa^2-((a/b)^2*sb^2)))
    }
  }
  
  # Reliability function
  reliability_function <- function(moe_vector, est_vector,flag_true) {
    if(flag_true==TRUE) {
      level_2_matrix <- matrix(data=2,nrow=nrow(moe_vector), ncol=length(geo_vector))
      level_1_matrix <- matrix(data=1,nrow=nrow(moe_vector), ncol=length(geo_vector))
    } else {
    level_2_matrix <- matrix(data=2,nrow=nrow(moe_vector), ncol=2*(agg_geo_end-dis_geo_start+1))
    level_1_matrix <- matrix(data=1,nrow=nrow(moe_vector), ncol=2*(agg_geo_end-dis_geo_start+1))
    }
    
    reliability <- moe_vector/(z_score*est_vector+.0000001)
    reliability <- as.matrix(reliability)*is.finite(as.matrix(reliability))
    reliability_2 <- pmin(2*floor(2*reliability),level_2_matrix,na.rm=TRUE)
    reliability_1 <- pmin(floor((10/3)*reliability),level_1_matrix, na.rm=TRUE)
    
    return(pmax(reliability_2,reliability_1))
  }
  ### End subroutines section
  
  setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/lookup_tables")
  lookup_table <- read.csv(lookup_table_name, sep=",", header=TRUE)
  
  # creating vector with row numbers of relevant data lines
  rel_data <- which(grepl(1,lookup_table$Relevant) & lookup_table$Line_number > 0)
  
  # User-error check:  Append the indices of all "universe totals" if they were not marked as relevant
  total_data <- rel_data
  total_data[]<- 0
  
  for(i in rel_data) {
    td_set <- which(rel_data == i)
    total_table <- lookup_table[which(lookup_table$Sequence.Number == lookup_table$Sequence.Number[[i]]),]
    total_table <- total_table$Table_line[which(total_table$Start.Position == lookup_table$Start.Position[[i]])]
    total_index <- as.numeric(total_table[3])
    if(lookup_table$Line_number[total_index]==0) {
      total_data[td_set] <- total_index +1
    } else total_data[td_set] <- total_index
  }
  
  total_data <- unique(total_data)
  rel_data <- sort(unique(c(rel_data,total_data)))
  ### End user-error check
  
  # Vector of percent names
  geo_percent <- paste(geo_names, c(rep("_pct", length(geo_names))), sep="")
  # Vector of margin of error names
  moe_names <- paste(geo_names, c(rep("_moe", length(geo_names))), sep="")
  # Vector of percent margin of error names
  pct_moe_names <- paste(geo_names, c(rep("_pct_moe", length(geo_names))), sep="")
  # Vector of reliability names
  est_rel_names <- paste(geo_names, c(rep("_rel", length(geo_names))), sep="")
  # Vector of percent reliability names
  pct_rel_names <- paste(geo_names, c(rep("_pct_rel", length(geo_names))), sep="")
  
  # setting start and end indices
  meta_start <- 1
  meta_end <- 6
  dis_geo_start <- 7
  dis_geo_end <- meta_end+length(geo_vector)
  agg_geo_start <- meta_end+length(geo_vector)+1
  agg_geo_end <- meta_end+length(geo_names)
  geo_pct_start <- meta_end+length(geo_names)+1
  geo_pct_end <- meta_end+2*length(geo_names)
  est_moe_start <- meta_end+2*length(geo_names)+1
  est_moe_end <- meta_end+2*length(geo_names)+length(geo_vector) # only includes discrete geographies
  pct_moe_start <- meta_end+3*length(geo_names)+1
  pct_moe_end <- meta_end+4*length(geo_names)
  est_rel_start <- meta_end+4*length(geo_names)+1
  est_rel_end <- meta_end+5*length(geo_names)
  pct_rel_start <- meta_end+5*length(geo_names)+1
  pct_rel_end <- meta_end+6*length(geo_names)
  
  z_score <- 1.645
  
  #initializing an empty data frame to store results
  results_df <- as.data.frame(matrix(nrow=length(rel_data),ncol=pct_rel_end))
  names(results_df) <- c("Table_Number", "Table_line", "Table_Title", "Universe","Datapoint","U_total_line",
                         geo_names, geo_percent, moe_names, pct_moe_names, est_rel_names, pct_rel_names)
  
  # reading in relevant tables and aggregating data
  current_table_name <- ""
  
  # setting folder for geographies or tracts and blocks
  if(granularity =="geographies") {
    setwd(paste("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/acs_ftp/geographies",year_folder,sep=""))  
  } else if(granularity == "tracts_blocks") {
    setwd(paste("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/acs_ftp/tracts_blocks",year_folder,sep=""))
  }

  # Begin loop
  for(i in rel_data) {
    rowset <- which(rel_data==i) #determines which row of the results_df to use
    
    data_table_name <- as.character(lookup_table$table_name[[i]]) # extract file name
    margin_table_name <- as.character(lookup_table$margin_name[[i]]) # extract moe file name
    
    # test if file name has changed.  If so, read in the new table; otherwise skip
    if(data_table_name != current_table_name) {
      data_table <- read.csv(data_table_name, sep=",", header=FALSE)
      indices <- seq.int(nrow(data_table))
      na_list <- which(is.na(match(data_table[,6],geo_vector)))
      data_lines <- indices[-na_list]
      data_lines <- data_lines[order(geo_vector)]
      data_table <- data_table[data_lines,]
      
      margin_table <- read.csv(margin_table_name, sep=",", header=FALSE)
      margin_table <- margin_table[data_lines,]
      
      current_table_name <- data_table_name
    }
    
    #identify universe
    uni_table <- lookup_table[which(lookup_table$Sequence.Number == lookup_table$Sequence.Number[[i]]),]
    uni_table <- uni_table$Table.Title[which(uni_table$Start.Position == lookup_table$Start.Position[[i]])]
    
    # record table title, universe, and data point
    results_df[rowset,1] <- as.character(lookup_table$Table.ID)[[i]]
    results_df[rowset,2] <- lookup_table$Table_line[[i]]
    results_df[rowset,3] <- as.character(uni_table[[1]])
    results_df[rowset,4] <- as.character(uni_table[[2]])
    results_df[rowset,5] <- as.character(lookup_table$Table.Title)[[i]]

    # transfer data from data_table to results data frame
    results_df[rowset,c(dis_geo_start:dis_geo_end)] <- as.numeric(as.character(data_table[,lookup_table$Line_number[i]]))
    results_df[rowset,c(est_moe_start:est_moe_end)] <- as.numeric(as.character(margin_table[,lookup_table$Line_number[[i]]]))
    
    # marking row with universe total for percentages
    total_table <- lookup_table[which(lookup_table$Sequence.Number == lookup_table$Sequence.Number[[i]]),]
    total_table <- total_table$Table_line[which(total_table$Start.Position == lookup_table$Start.Position[[i]])]
    total_index <- as.numeric(total_table[3])
    if(lookup_table$Line_number[total_index]==0) {
      results_df[rowset,6] <- which(results_df$Table_line == total_table[3]+1)
    } else results_df[rowset,6] <- which(results_df$Table_line == total_table[3])
    
  }
  
  ### End of loop
  
  # flagging non-aggregatable data
  flag_vector <- which(grepl("MEDIAN",results_df[,3]) | grepl("QUINTILE", results_df[,3]) | grepl("PER CAPITA", results_df[,3]))
  
  # Aggregating geographies to regions when there ARE flagged variables
  if(length(flag_vector) > 0) {
    if(length(agg_vector) > 0){
      agg_length <- length(agg_vector)/2
      for(i in 1:agg_length) {
        k <- 2*i - 1
        results_df[-flag_vector,agg_geo_start+k-1] <- rowSums(results_df[-flag_vector,c((meta_end+agg_vector[k]):(meta_end+agg_vector[k+1]))])
        results_df[-flag_vector,est_moe_start+length(geo_vector)+k-1] <- apply(results_df[-flag_vector,c((est_moe_start+agg_vector[k]-1):(est_moe_start+agg_vector[k+1]-1))],1,euclidnorm)
      }
    }
      
    # Computing percentages
    results_df[-flag_vector,c(geo_pct_start:geo_pct_end)] <- results_df[-flag_vector,c(dis_geo_start:agg_geo_end)]/results_df[results_df[-flag_vector,6],c(dis_geo_start:agg_geo_end)]
    
    # Computing percent margins of error
    for (i in 1:length(geo_names)) {
      j <- meta_end+i
      comp_matrix <- cbind(results_df[-flag_vector,j],results_df[results_df[-flag_vector,6],j],
                       results_df[-flag_vector,est_moe_start+i-1],
                       results_df[results_df[-flag_vector,6],est_moe_start+i-1])
      results_df[-flag_vector,pct_moe_start+i-1] <- apply(comp_matrix,1,pct_moe_function)
    }
    
    # Computing reliability scores for non-flag vector variables
    results_df[-flag_vector,c(est_rel_start:pct_rel_end)] <- reliability_function(results_df[-flag_vector,c(est_moe_start:pct_moe_end)],
                                                                                  results_df[-flag_vector,c(dis_geo_start:geo_pct_end)],FALSE)
    
    # Computing reliability scores for flag vector variables
    results_df[flag_vector,c(est_rel_start:(est_rel_end-(agg_geo_end-dis_geo_end)))] <- reliability_function(results_df[flag_vector,c(est_moe_start:est_moe_end)],
                                                                                                               results_df[flag_vector,c(dis_geo_start:dis_geo_end)],TRUE)
    # -99 for medians, quintiles, per capita
    results_df[flag_vector,c(agg_geo_start:geo_pct_end,(pct_moe_start-length(agg_vector)/2):pct_moe_end,(pct_rel_start-length(agg_vector)/2):pct_rel_end)] <- -99
    
  } else { 
    
    ### Aggregating geographies to regions when there are NO flagged variables
    if(length(agg_vector) > 0) {
      agg_length <- length(agg_vector)/2
      for(i in 1:agg_length) {
        k <- 2*i - 1
        results_df[,agg_geo_start+i-1] <- rowSums(results_df[,c((meta_end+agg_vector[k]):(meta_end+agg_vector[k+1]))])
        results_df[,geo_pct_end+length(geo_vector)+i] <- apply(results_df[,c((geo_pct_end+agg_vector[k]):(geo_pct_end+agg_vector[k+1]))],1,euclidnorm)
      }
    }
      
    # Computing percentages
    results_df[,c(geo_pct_start:geo_pct_end)] <- results_df[,c(dis_geo_start:agg_geo_end)]/results_df[results_df[,6],c(dis_geo_start:agg_geo_end)]
    
    # Computing percent margin of error
    for (i in 1:length(geo_names)) {
      j <- meta_end+i
      comp_matrix <- cbind(results_df[,j],results_df[results_df[,6],j],
                           results_df[,est_moe_start+i-1],
                           results_df[results_df[,6],est_moe_start+i-1])
      results_df[,pct_moe_start+i-1] <- apply(comp_matrix,1,pct_moe_function)
    }
    
    # Computing reliability scores
    results_df[,c(est_rel_start:pct_rel_end)] <- reliability_function(results_df[,c(est_moe_start:pct_moe_end)], results_df[,c(dis_geo_start:geo_pct_end)],FALSE)
  }
  
  return(results_df)
}
