#this script iteratively runs coordinatecleaner on all of our data
#this is a standard step in dealing with gbif data, removing out likely errors
#note, it as been set up to run as a batched array, and is deployed by the following bash script: 01A_coordinate_cleaner_array_bash.sh
#we run 300 species at a time, in chunks across multiple arrays run simultaneously. This significantly speeds up the process

#######################################################################################
#setup array logic & chunking
#######################################################################################
#read the task id from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

#define chunk parameters i.e, how many species are run in each array
#here we're running 300 species in a chunk
chunk_size <- 300
#list all the files that have been downloaded from gbif
#path to output from the 02_data_split.sh script in Step1_data_download section
all.sp <- list.files("/scratch/kricks/GBIF_store/processed_Rhizobia/subset_data")
total_species <- length(all.sp) 

#determine start and end indices for this specific chunk
start_index <- ((task_id - 1) * chunk_size) + 1
end_index <- min(task_id * chunk_size, total_species)

#if trying to run specific set of species, change this code
#start_index <- 7770 
#end_index <- 7800

###############################################################################################
###############################################################################################
###############################################################################################
#load in libraries
library(dplyr)
library(rgbif)
library(sf)
library(data.table)
library(CoordinateCleaner)
library(rnaturalearth)

#path to output from the 02_data_split.sh script in 01_data_download section
setwd("/scratch/kricks/GBIF_store/processed_Rhizobia/subset_data")
all.sp <- list.files()
for(i in start_index:end_index){
  temp.d <- read.csv(all.sp[i]) #read in data and massage to a usable format
  temp.d$decimalLatitude <- as.numeric(temp.d$decimalLatitude)
  temp.d$decimalLongitude <- as.numeric(temp.d$decimalLongitude)
  temp.d <- temp.d[!is.na(temp.d$decimalLatitude) & !is.na(temp.d$decimalLongitude), ]
  
  #if csv file is empty, skip and move on
  if(nrow(temp.d) == 0) {
    print(paste("Skipping empty file:", all.sp[i]))
    next
  }
  
  #coordinate cleaner command
  #determines which coordinates should be flagged, based on classic problems
  #can edit to change flags
  #the following are fairly standard 
  #removes occurences at the exact centroid of a country, that have equal lat/long, that are located the gbif or instituations
  #are located in oceans or lakes, are at zero lat/long, or are outliers (based on being over 3000 km from the closest sample)
  flags <- clean_coordinates(x = temp.d, 
                             lon = "decimalLongitude", 
                             lat = "decimalLatitude",
                             tests = c("centroids", "equal", "gbif", "institutions", "seas", "zeros", "outliers"),
                             outliers_method = 'distance',
                             outliers_td = 3000,
                             verbose=F)
  
  #remove those flagged by coordinate cleaner
  #if the entire species has been flagged, skip to the next line
  dat_fl <- temp.d[flags$.summary,]
  if(nrow(dat_fl) == 0) {
    print(paste("All points flagged, skipping:", all.sp[i]))
    next
  }
  
  #additionally remove species that map back to antartica
  #load world map
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  points_sf <- st_as_sf(dat_fl, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  
  #assign each point to a continent and remove antartic samples
  cont.list <- st_drop_geometry(st_join(points_sf, world_map["continent"], join = st_nearest_feature))
  dat_fl <- dat_fl[is.na(cont.list$continent) | cont.list$continent != "Antarctica", ]
  
  #export data
  write.csv(dat_fl, paste("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/", all.sp[i], 
                          sep=""),row.names=FALSE)
  
  print(i)
}


