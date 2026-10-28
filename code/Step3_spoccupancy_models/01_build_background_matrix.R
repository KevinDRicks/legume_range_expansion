#this script creates a grid across the global that measures sampling intensity over time 
#of legumes in the gbif dataset. This can used as a measure of "effort" put into sampling 
#at any one location over time, which is an important component in the detection sub-model 
#in the spoccupancy framework.

#this script can be slow, so it can either be submitted as a  slurm script
#for ease, we simply used an interactive node on the computing cluster (see below for set up)
#salloc --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=20G --time=04:00:00 
#module load r
#R

.libPaths("/scratch/kricks/Rlibs/Juptyer_R/")
#load in libraries
library(data.table)
library(terra)
library(dplyr)
library(sf)

#path to the output from the 01B_coordinate_cleaner_array_format.R script from Step2_data_cleaning
setwd("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data")

#define spatial grid using sf
#grid is encompasses whole globe, with cells sized at ~100 km wide
globe_bbox <- st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90), crs = 4326)
global_grid_sf <- st_make_grid(globe_bbox, cellsize = c(0.9, 0.9)) %>%
  st_sf(geometry = .) %>%
  mutate(cell_id = row_number())

#global grid cell id
global_grid_sf$cell_id <- paste("cell",1:nrow(global_grid_sf),sep="_")

#extract centroids of grids
grid_centroids <- st_centroid(global_grid_sf) %>% st_coordinates() %>% as.data.frame()
colnames(grid_centroids) <- c("center_lon", "center_lat")
row.names(grid_centroids) <- paste("cell",1:nrow(grid_centroids),sep="_")


#create master dataframe to store aggregate sampling effort
master_effort <- data.frame(cell_id = integer(), Bin = numeric(), SubBin = numeric(), N = integer())

#list all species CSVs
csv_files <- list.files()
#build an empty list to fill with data from each species
effort_list <- vector("list", length(csv_files))
#loop through all species to sum the total number of occurences in each cell and time bin
for (f in 1:length(csv_files)) {
  
  #read in species i
  sp_data <- read.csv(csv_files[f]) %>%
    select(any_of(c("decimalLongitude", "decimalLatitude", "year"))) %>%
    filter(!is.na(year) & year >= 1900)
  
  #put occurences into 6 temporal bins, and then 5 year subbins within each of those
  #these represent temporal bins in which we characterize occupancy
  #remove occurences that have outlier coordinates and have no time stamp
  sp_data <- sp_data %>%
    mutate(
      
      SubBin = floor((year - 1) / 5) * 5 + 1,
      Bin = case_when(
        year < 1966 ~ 1,
        year >= 1966 & year < 1981 ~ 2,
        year >= 1981 & year < 1996 ~ 3,
        year >= 1996 & year < 2006 ~ 4,
        year >= 2006 & year < 2016 ~ 5,
        year >= 2016 ~ 6,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(Bin)) %>%
    filter(decimalLongitude >= -180 & decimalLongitude <= 180 &
             decimalLatitude >= -90 & decimalLatitude <= 90)
  sp_data <- data.frame(sp_data)
  
  #sp_data[sp_data$Bin==6,]
  #table(sp_data[sp_data$year >= 1900 & sp_data$year < 1966,]$SubBin)
  
  #for each occurrence, map back to our grid
  if(nrow(sp_data)>0){
    
    #convert points to sf object
    points_sf <- st_as_sf(sp_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    #use st_intersects to get a list of indices
    intersections <- st_intersects(points_sf, global_grid_sf)
    
    #extract only the first match for every point
    sp_data$cell_id <- global_grid_sf$cell_id[sapply(intersections, function(x) if (length(x) > 0) x[1] else NA)]
    
    file_effort <- sp_data %>%
      filter(!is.na(cell_id)) %>%
      group_by(cell_id, Bin, SubBin) %>%
      summarise(N = n(), .groups = "drop")
    
    #pace the resulting dataframe into the list
    effort_list[[f]] <- file_effort
  }
  
  print(f)
}

#join all species together
master_effort <- bind_rows(effort_list) %>%
  group_by(cell_id, Bin, SubBin) %>%
  summarise(N = sum(N), .groups = "drop")
master_effort


#convert the row names of your centroids dataframe into a proper column
grid_centroids$cell_id <- rownames(grid_centroids)

#join the coordinates onto effort matrix
master_effort_coords <- master_effort %>%
  left_join(grid_centroids, by = "cell_id")

head(master_effort_coords)

#output data
write.csv(data.frame(master_effort_coords),"/scratch/kricks/GBIF_store/processed_Rhizobia/global_effort.csv",row.names=F)



#we can visualize sampling effort here
master_effort <- read.csv("/scratch/kricks/processed_Rhizobia/Rhizobia/global_effort.csv")
library(dplyr)
library(ggplot2)
library(maps)

#aggregate effort data
global_aggregated_effort <- master_effort %>%
  mutate(cell_num = as.numeric(gsub("cell_", "", cell_id))) %>%
  group_by(cell_num) %>%
  summarise(Total_Sampling_Intensity = sum(N), .groups = "drop")

#create global grid
res <- 0.9
ncol_global <- 360 / res
nrow_global <- 180 / res

#generates index grid
full_grid <- data.frame(
  cell_num = 1:(ncol_global * nrow_global),
  col_idx = rep(1:ncol_global, times = nrow_global),
  row_idx = rep(1:nrow_global, each = ncol_global)
) %>%
  mutate(
    #convert to center coordinates
    center_lon = (col_idx - 1) * res + (res / 2) - 180,
    center_lat = (row_idx - 1) * res + (res / 2) - 90
  )

#merge effort matrix onto the full grid and fill NA with 0
map_data <- full_grid %>%
  left_join(global_aggregated_effort, by = "cell_num") %>%
  mutate(Total_Sampling_Intensity = ifelse(is.na(Total_Sampling_Intensity), 0, Total_Sampling_Intensity))


#generate map
#set log scale as data is highly skewed.
library(ggplot2)
world_map <- borders("world", colour = "gray30", size = 0.2)

sampling_map <- ggplot() +
  geom_raster(data = map_data, 
              aes(x = center_lon, y = center_lat, fill = Total_Sampling_Intensity)) +
  scale_fill_viridis_c(option = "plasma", 
                       trans = "log10",
                       na.value = "black", 
                       name = "Observations",
                       labels = scales::label_comma()) +
  world_map +
  coord_sf(expand = FALSE) + 
  labs(x = "Longitude", y = "Latitude")
sampling_map





