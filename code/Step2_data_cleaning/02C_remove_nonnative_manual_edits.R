#a small number of species fail in the 02B_remove_nonnative_array.R have failed. We do a manual inspection on these species

#load libraries
library('readr')
library('dplyr')
library('rWCVP')
library('sf')
library('rnaturalearth')
library('units')
library('terra')
library('rgbif')


#read in the file determing the adjacent regions, and conduct some minor
#data massaging
mapper.region <- read.csv('/home/kricks/codebases/range_expansion/Nitrogen_fixation/region_adj.csv',fileEncoding = "latin1")
rownames(mapper.region) <- mapper.region[,1]
mapper.region <- mapper.region[,-1]

#mirror this on the diaganol
mat <- as.matrix(mapper.region)
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
mapper.region_full <- as.data.frame(mat)
rownames(mapper.region_full) <- gsub("\u00A0", " ", rownames(mapper.region_full))


#Download and read the official TDWG Level 2 Regions directly from GitHub
#will be directly useful for determining a species native range 
tdwg_regions <- st_read(
  "https://raw.githubusercontent.com/tdwg/wgsrpd/master/geojson/level2.geojson",
  quiet = TRUE
)

#read in key garden datasets
local_wcvp_names <- read_delim("/scratch/kricks/GBIF_store/wcvp_raw/wcvp_names.csv", delim = "|", guess_max = 100000)
local_wcvp_dist <- read_delim("/scratch/kricks/GBIF_store/wcvp_raw/wcvp_distribution.csv", delim = "|", guess_max = 100000) %>%
  mutate(plant_name_id = as.character(plant_name_id))











my_species1 <- list.files('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/')
my_species2 <- list.files('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/')


length(my_species1)
length(my_species2)


fail.files <- list.files('/scratch/kricks/GBIF_store/processed_Rhizobia/skipped_files/')

fail.vec <- NULL
for(i in 1:length(fail.files)){
  fail.vec <- c(fail.vec,unlist(read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/skipped_files/',fail.files[i]))))
}




#failed species
species.int <- list.files("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data")
cleaned.files <- list.files("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2")
missing.d <- setdiff(species.int,cleaned.files)
missing.d


################################################################################################################
#sp1
################################################################################################################
i<-2943160
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data


extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#all invasive-remove
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/Rhizobia/cleaned_data_step2/",missing.d[i]),row.names=F)

################################################################################################################
#sp2
################################################################################################################
i<-2943167
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data
results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp3
################################################################################################################
i<-2943171
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#europe is invasive-manual edit
temp.d <- temp.d[temp.d$decimalLongitude > 50,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp4
################################################################################################################
i<-2943179
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp5
################################################################################################################
i<-2943184
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#invasive
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp6
################################################################################################################
i<-2943186
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#manual edit
temp.d <- temp.d[temp.d$decimalLongitude>50,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp7
################################################################################################################
i<-2955682
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#failed to match wcvp
#online search says native to asia
#manual edit
temp.d <- temp.d[temp.d$decimalLongitude > 50,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp8
################################################################################################################
i<-2956441
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#not a real species, garden hybrid, remove
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp9
################################################################################################################
#i<-2972031
#maps::map("world");axis(1);axis(2)
#temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
#points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')
#
#results_list <- name_usage(key = i)$data
#
#extracted_names <- cbind(results_list$key,results_list$canonicalName)
#
#final_df <- data.frame(extracted_names)
#colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i
#
#matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
#matched_names
#my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
#temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
##extract distributions
#species_distributions <- temp.sp %>%
#  #Convert the wcvp_accepted_id to character so it matches the distribution data
#  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
#  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
#  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
#species_distributions
#fuzzy mutch with only 1 individual, and failures
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp10
################################################################################################################
i<-2972034
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp11
################################################################################################################
i<-2972043
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#somewhat a mess, human use, remove
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp12
################################################################################################################
i<-2972060
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp13
################################################################################################################
i<-2972061
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failure
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp14
################################################################################################################
i<-2972062
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp15
################################################################################################################
i<-2972067
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failed
#remove out of africa samples
temp.d <- temp.d[temp.d$decimalLongitude < 60,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp16
################################################################################################################
i<-2972073
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failure
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp17
################################################################################################################
i<-2972077
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failure
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp18
################################################################################################################
i<-2978045
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual wcvp edit
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp19
################################################################################################################
i<-3908215
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual wcvp edit 
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp20
################################################################################################################
i<-3908221
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual wcvp edit
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp21
################################################################################################################
i<-3908225
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual edit
temp.d <- temp.d[temp.d$decimalLongitude < -20,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)

################################################################################################################
#sp22
################################################################################################################
i<-3913178
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual edit
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp23
################################################################################################################
i<-3913180
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual wcvp edit
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp24
################################################################################################################
i<-3913185
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#manual edit to remove invasive, missing from wcvp
temp.d <- temp.d[temp.d$decimalLongitude < -10,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp25
################################################################################################################
i<-3913190
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#missing from wcvp
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp26
################################################################################################################
i<-3968870
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#removed, few point
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp27
################################################################################################################
i<-5343326
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#only invasive
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp28
################################################################################################################
i<-5343327
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp29
################################################################################################################
i<-5343332
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
temp.d <- temp.d[temp.d$decimalLongitude > 25,]
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp30
################################################################################################################
i<-5343342
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp31
################################################################################################################
i<-5343357
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp32
################################################################################################################
i<-5343362
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp33
################################################################################################################
i<-5343364
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp34
################################################################################################################
i<-5343368
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp35
################################################################################################################
i<-5343371
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp36
################################################################################################################
i<-5611691
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp37
################################################################################################################
i<-5612275
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#remove
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp38
################################################################################################################
i<-5628840
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#missing wcvp
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp39
################################################################################################################
i<-5631191
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp40
################################################################################################################
i<-5631918
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#unclear why this failed
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp41
################################################################################################################
i<-8213866
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data

extracted_names <- cbind(results_list$key,results_list$canonicalName)

final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
matched_names
my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
species_distributions
#remove artificial hybrid
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)


################################################################################################################
#sp42
################################################################################################################
i<-8347080
maps::map("world");axis(1);axis(2)
temp.d <- read.csv(paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data/",i,".csv"))
points(decimalLatitude~decimalLongitude,temp.d,pch=16,col='firebrick')

results_list <- name_usage(key = i)$data


extracted_names <- cbind(results_list$key,results_list$canonicalName)
final_df <- data.frame(extracted_names)
colnames(final_df) <- c('gbif_key','gbif_name');final_df$gbif_key <- i

#matched_names <- wcvp_match_names(names_df =final_df,name_col = "gbif_name",wcvp_names = local_wcvp_names)
#matched_names
#my_species_range_list <- data.frame(final_df,Home_range=NA, Invasive_range=NA)
#temp.sp <- matched_names[matched_names$gbif_key==my_species_range_list$gbif_key,]
#extract distributions
#species_distributions <- temp.sp %>%
  #Convert the wcvp_accepted_id to character so it matches the distribution data
#  mutate(wcvp_accepted_id = as.character(wcvp_accepted_id)) %>%
#  inner_join(local_wcvp_dist, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
#  select(gbif_key, wcvp_name, area_code_l3, area = region, introduced, extinct)
#species_distributions
#read.csv(paste0('/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/',missing.d[i]),row.names=F)
#fuzzy match causes failure, but its in the right location
write.csv(temp.d,paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/",i,".csv"),row.names=F)

