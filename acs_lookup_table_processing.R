# lookup_table_processing
# Peter Banwarth
# 12/13/2016

# arguments
# estub:  The string that identifies year, 1-year or 5-year estimates, and state for estimates
# mstub:  The string that identifies year, 1-year or 5-year estimates, and state margins of error
# lk_input: The raw data table file name (string ending in .csv):  "rawdatatable.csv"
# lk_output: The processed data table file name (string ending in .csv):  "lk_tble_template_year.csv"

acs_lookup_table_processing <- function(estub,mstub,lk_input,lk_output) {

  library(data.table)
  
  setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/raw_lookup_tables")
  
  # read in lookup table
  lookup_table <- read.csv(lk_input, sep=",", header=TRUE, colClasses = c("Sequence.Number"="character"))
  ### IMPORTANT!  CONFIRM THAT NUMBER OF ROWS IN READ TABLE IS CORRECT ###
  
  # keep only desired columns
  lookup_table <- lookup_table[,c(2:5,8)]
  
  # create table name column
  lookup_table[,6] <- paste(c(rep(estub, length(lookup_table[,2]))), lookup_table[,2], c(rep("000.txt", length(lookup_table[,2]))), sep="")
  lookup_table[,7] <- paste(c(rep(mstub, length(lookup_table[,2]))), lookup_table[,2], c(rep("000.txt", length(lookup_table[,2]))), sep="")
 
  # create table line column and empty columns for indexing
  lookup_table <- cbind(lookup_table,seq.int(nrow(lookup_table)),c(rep(0,nrow(lookup_table))),c(rep(0,nrow(lookup_table))))
  
  # convert sequence number to numeric
  lookup_table[,2] <- as.numeric(lookup_table[,2])
  
  # identify which start positions are not blank
  na_list <- which(is.na(lookup_table[,4]))
  indices <- seq.int(nrow(lookup_table))
  num_list <- indices[-na_list]

  # create indexing columns - start position and line number
  lookup_table[num_list,9]<-lookup_table[num_list,4]
  # lookup_table[num_list_no_med_type+2,10]<-lookup_table[num_list_no_med_type,4]
  # lookup_table[num_list_med_type+3,10] <- lookup_table[num_list_med_type,4]

  # Sequencing line numbers
  for(i in num_list[-length(num_list)]) {
    j <- num_list[which(num_list == i)+1]
    lookup_table[c((i+1):(j-1)),9] <- lookup_table[i,4]
    if(j-i>2) {
      lookup_table[c((i+2):(j-1)),10] <- lookup_table[i,9]+seq.int(j-i-2)-1
    }
  }
  
  # Finishing last set since it doesn't work within confines of for loop
  lookup_table[num_list[length(num_list)]:nrow(lookup_table),9] <- lookup_table[num_list[length(num_list)],4]
  lookup_table[(num_list[length(num_list)]+2):nrow(lookup_table),10] <- lookup_table[num_list[length(num_list)],4]+seq.int(nrow(lookup_table)-num_list[length(num_list)]-1)-1
  
  # Fixing "--" indices
  med_type_list <- which(grepl("*--$",lookup_table$Table.Title))
  
  for(i in med_type_list){
    table_id <- lookup_table$Table.ID[i]
    start_pos <- lookup_table[i,9]
    alt_indices <- lookup_table[which(lookup_table$Table.ID==table_id & lookup_table[,9]==start_pos),8]
    end_pos <- alt_indices[length(alt_indices)]
    lookup_table[i,10] <- 0
    lookup_table[c((i+1):end_pos),10] <- lookup_table[c((i+1):end_pos),10]-1
  }
  
  # Fixing Quintile indices - This needs to be checked if a new Lookup Table is downloaded
  if(length(which(grepl("Quintile",lookup_table$Table.Title))) > 0) {
    quintile_indices <- which(grepl("Quintile",lookup_table$Table.Title))
    lookup_table[quintile_indices[c(1,6,12)],10] <- 0
    lookup_table[c(quintile_indices[-c(1,6,12)],quintile_indices[c(5,11,17)]+1),10] <- lookup_table[c(quintile_indices[-c(1,6,12)],quintile_indices[c(5,11,17)]+1),10]-1
  }
  
  # reorder columns and assign names
  lookup_table <- lookup_table[c(8,1,2,6,7,10,9,5)]
  names(lookup_table) <- c("Table_line","Table.ID","Sequence.Number","table_name","margin_name","Line_number","Start.Position","Table.Title")
  
  # save processed table
  setwd("L:/Health/Epidemiology/Banwarth_Epi/ACS_tables/lookup_tables")
  write.table(lookup_table, file=lk_output, sep=",",row.names=FALSE)
}
