#this script identifies the species of interest in our dataset that have sufficient replication
#this is identical to the 03B_identify_species_of_interest.R script, though runs it on the species
#that haven't already gone through the invasive species filter. it outputs a .csv species_interest.csv
#representing the species we may want to run occupancy models on

#load libraries
library(data.table)

#now lets check to see which ones no longer pass our threshold 
#files output from the coordinatecleaner script (prior to invasive occurence filters)
setwd('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data')
all_files <- list.files(pattern = "\\.csv$")
total_files <- length(all_files)
total_files
yr_breaks <- c(1900, 1965, 1980, 1995, 2005, 2015, 2025)
yr_labels <- c("pre_1965", "1966_1980", "1981_1995", "1996_2005", "2006_2015", "2016_2025")

process_file <- function(i) {
  
  file_path <- all_files[i]
  
  pct_done <- (i / total_files) * 100
  cat(sprintf("\rProcessing file %d of %d (%.1f%% done)", i, total_files, pct_done))
  
  dt <- fread(file_path, select = c("year", "decimalLongitude"))
  
  dt[, year := suppressWarnings(as.numeric(year))]
  dt[, decimalLongitude := suppressWarnings(as.numeric(decimalLongitude))]
  
  dt <- dt[!is.na(year) & !is.na(decimalLongitude)]
  
  binned_years <- cut(dt$year, breaks = yr_breaks, labels = yr_labels, right = TRUE)
  counts <- as.list(table(binned_years))
  
  counts$filename <- file_path 
  
  return(counts)
}

results_list <- lapply(seq_along(all_files), process_file)

sum.df_alt <- rbindlist(results_list, fill = TRUE)
print(sum.df_alt)

sum.df_alt2 <- data.frame(sum.df_alt[,-7])



sum(rowSums(sum.df_alt2>1)>=6)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=5)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=4)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=3)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=2)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=1)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=0)/nrow(sum.df_alt2)

sum(rowSums(sum.df_alt2)>50)


sum(rowSums(sum.df_alt2>1)>=6 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=5 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=4 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)#this is the sweet spot-67.3%
sum(rowSums(sum.df_alt2>1)>=3 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=2 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=1 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)
sum(rowSums(sum.df_alt2>1)>=0 & rowSums(sum.df_alt2)>50)/nrow(sum.df_alt2)

#this is the list of species of interest for the next steps of the project
#this includes speces that have at least 2 records in at least 4 of the 6 time periods, and have more than 50 records total
species_interest <- all_files[rowSums(sum.df_alt2>1)>=4 & rowSums(sum.df_alt2)>50]
write.csv(data.frame(Species_file_name=species_interest), '/scratch/kricks/GBIF_store/processed_Rhizobia/species_interest.csv', row.names=F)
