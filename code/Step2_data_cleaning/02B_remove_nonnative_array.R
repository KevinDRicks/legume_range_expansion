#this script iteratively runs to prune out invasive occurences from the datasets
#note, it as been set up to run as a batched array, and is deployed by the following bash script: 02A_remobe_nonnative_bash.sh
#the basic logic is we compare each occurence against what is noted at invasive or native from the WCVP dataset
#however, this can be a bit fuzzy and can be based on different regional descriptions of invasive vs native
#for example, in going through the data, one species was considered native in british columbia, but invasive in washington state
#we've developed a custom script to address this, considering areas adjacent to native regions as native as well

#######################################################################################
#setup array logic & chunking
#######################################################################################
#read the task id from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

#define chunk parameters i.e, how many species are run in each array
#here we're running 10 species in a chunk
chunk_size <- 10
#read in a list of species files to look at
#this is a path to the output from the 01B_coordinate_cleaner_array_format.R script from the Step2_data_cleaning section
my_species <- list.files('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/')
length(my_species)
total_species <- length(my_species) # 8721 separate species

#determine start and end indices for this specific chunk
start_index <- ((task_id - 1) * chunk_size) + 1
end_index <- min(task_id * chunk_size, total_species)


#if trying to run specific set of species, change this code
#start_index <- 7770 
#end_index <- 7800


#load in libraries
library('rgbif')
library('readr')
library('dplyr')
library('rWCVP')
library('sf')
library('rnaturalearth')
library('units')
library('terra')

#extract out species names based on gbif titles for the small subset that will be run in this part of the array
#these will be used in refrence to the WCVP database
my_keys <- gsub(".csv","",my_species[c(start_index:end_index)])
results_list <- lapply(my_keys, function(k) {
  name_usage(key = k)$data
})
extracted_names <- lapply(results_list, function(df) {
  df[, c("key", "canonicalName")]
})
final_df <- data.frame(do.call(rbind, extracted_names))
colnames(final_df) <- c('gbif_key','gbif_name')
final_df$key <- my_keys
###############################################################################################
###############################################################################################
###############################################################################################
#read in the file determing the adjacent regions, and conduct some minor data massaging
#the region_adj.csv is a custom file we've built. It describes regions that are adjacent to each other
#to determine where native and invasive may breakdown
mapper.region <- read.csv('/home/kricks/codebases/range_expansion/Nitrogen_fixation/region_adj.csv',fileEncoding = "latin1")
rownames(mapper.region) <- mapper.region[,1]
mapper.region <- mapper.region[,-1]

#mirror this on the diagonal
mat <- as.matrix(mapper.region)
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
mapper.region_full <- as.data.frame(mat)
rownames(mapper.region_full) <- gsub("\u00A0", " ", rownames(mapper.region_full))

#download and read the official TDWG Level 2 Regions directly from GitHub
#will be directly useful for determining a species native range 
tdwg_regions <- st_read(
  "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level2.geojson",
  quiet = TRUE
)

#read in kew garden datasets
local_wcvp_names <- read_delim("/scratch/kricks/GBIF_store/wcvp_raw/wcvp_names.csv", delim = "|", guess_max = 100000)
local_wcvp_dist <- read_delim("/scratch/kricks/GBIF_store/wcvp_raw/wcvp_distribution.csv", delim = "|", guess_max = 100000) %>%
  mutate(plant_name_id = as.character(plant_name_id))


#match species names against the database 
matched_names <- wcvp_match_names(names_df =final_df,
                                  name_col = "gbif_name",
                                  wcvp_names = local_wcvp_names)
#table(table(matched_names$gbif_key))#some species might have more than one match in here
#table(matched_names$gbif_key)[which(table(matched_names$gbif_key)>1)]

#extract list of native areas for each species
#run in a loop for each of the 10 species in this array
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
skipped_for_manuel_edit <- NULL
for(i in 1:nrow(my_species_range_list)){
  #subset to species of interest
  temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key[i],]
  
  #extract distributions
  species_distributions <- temp.sp %>%
    #Convert the wcvp_accepted_id to character so it matches the distribution data
    mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
    inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
    select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
  
  #specifically extract native regions
  regions <- unique(species_distributions[species_distributions$introduced==0,]$area)  
  if(length(regions)>0){
    
    #now extract all the regions adjacent to these native areas that could 
    #realistically be within dispersal range while not being invasive
    region_adj <- regions
    for(j in 1:length(regions)){
      adj.region <- rownames(mapper.region_full)[mapper.region_full[rownames(mapper.region_full) == regions[j], ]=='Yes']
      region_adj <- c(region_adj,adj.region)
    }
    region_adj <- unique(region_adj)
    
    #need to rename central Asia as our database is named middle asia
    region_adj[region_adj=='Central Asia'] <- 'Middle Asia'
    
    #extracts polygons of these regions
    species_polygons <- tdwg_regions %>%
      filter(LEVEL2_NAM %in% region_adj[])
    
    
    polys_vect <- vect(species_polygons)
    polys_buffered <- buffer(polys_vect, width = 100000)
    
    polys_merged <- aggregate(polys_buffered)
    
    #read in actual data of interest
    latlong.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",my_species_range_list$gbif_key[i],".csv"))
    
    points_vect <- vect(latlong.d, 
                        geom = c("decimalLongitude", "decimalLatitude"), 
                        crs = "EPSG:4326")
    
    is_inside <- rowSums(relate(points_vect, polys_merged, "intersects")) > 0
    
    subset.d <- latlong.d[is_inside, ]  
    #plot(polys_merged)
    #plot(latlong.d$decimalLatitude~latlong.d$decimalLongitude)
    
    if(nrow(subset.d) == 0) {
      print(paste("All points removed, skipping:", all.sp[i]))
      next
    }
    
    #output data
    write.csv(subset.d, paste("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/", 
                              my_species_range_list$gbif_key[i],".csv",sep=""), row.names=FALSE)
    
  }
  if(length(regions)==0){
    skipped_for_manuel_edit <- c(skipped_for_manuel_edit,
                                 my_species_range_list$gbif_key[i])
  }
  print(i)
  
}

#some of these fail. in the next script, 02C, we manually inspect this small list of failures
if(length(skipped_for_manuel_edit)>0){
  write.csv(skipped_for_manuel_edit, 
            paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/skipped_files/", task_id,".csv"),
            row.names=FALSE)  
  
}




