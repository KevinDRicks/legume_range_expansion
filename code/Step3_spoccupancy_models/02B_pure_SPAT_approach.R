#this script iteratively builds occupancy models for each species. note, it as been set up
#to run as a batched array, and is deployed by the following bash script: 02A_pure_SPAT_array.sh


#######################################################################################
#load in necessary packages
#######################################################################################
library(spOccupancy)
library(coda)
library(ggplot2)
library(lubridate)
library(tidyr)
library(data.table)
library(dplyr)
library(sf)
library(CoordinateCleaner)
library(terra)
library(splines)


#######################################################################################
#load background effort matrix and list of species
#######################################################################################
#read in the background effort matrix-path to the output from the 01_build_background_matrix.R script
bg_effort <- data.frame(fread("/scratch/kricks/GBIF_store/processed_Rhizobia/global_effort.csv")) 

#bg_effort matrix should have unique cell identifier 
valid_cells <- unique(bg_effort$cell_id)
#read in species of interest that we want to process-output from the the 03B_identify_species_of_interest_step2.R script
species.int <- read.csv("/scratch/kricks/GBIF_store/processed_Rhizobia/species_interest_step2.csv")

#######################################################################################
#setup array logic
#######################################################################################
#read in the task ID from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

#######################################################################################
#each script will only do a single species
#######################################################################################
#this is the path to species data, generated in the 02B_remove_nonnative_array.R script
setwd("/scratch/kricks/GBIF_store/processed_Rhizobia/cleaned_data_step2/")
species <- task_id
species.id <- gsub(".csv","",species.int[species,])

#read in species data
#place occurences into same temporal bins used in building background effort database
#we only want occurences with known lat/long coordinates and data post 1900
working.d <- read.csv(species.int[species,]) %>%
  select(any_of(c("decimalLongitude", "decimalLatitude", "year"))) %>%
  filter(!is.na(year) & year >= 1900)

#reshape data, placing occurrences in temporal bins as well as subbins, based on 5 year time spans
working.d <- working.d %>%
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
  filter(decimalLongitude >= -180 & decimalLongitude <= 180 & #remove out occurences not mapping to normal coordinates
           decimalLatitude >= -90 & decimalLatitude <= 90)
working.d <- data.frame(working.d)
dim(working.d)



###############################################################################
###############################################################################
#recreate the fixed global grid that was used to make the background sampling matrix 
#we're going to assign all of our points to this coordinate grid
globe_bbox <- st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90), crs = 4326)
global_grid_sf <- st_make_grid(globe_bbox, cellsize = c(0.9, 0.9)) %>%
  st_sf(geometry = .) %>%
  mutate(cell_id = row_number())

global_grid_sf$cell_id <- paste("cell",1:nrow(global_grid_sf),sep="_")

#extract centroids of grids
grid_centroids <- st_centroid(global_grid_sf) %>% st_coordinates() %>% as.data.frame()
colnames(grid_centroids) <- c("center_lon", "center_lat")
row.names(grid_centroids) <- paste("cell",1:nrow(grid_centroids),sep="_")

#cells in grid_centroids should match the cells in your background effort matrix
#look at small subset
grid_centroids[row.names(grid_centroids) %in% bg_effort$cell_id[c(1:10)],]
bg_effort[c(1:10),]

#assign cell ids to our observations
#first, convert our observed data into an sf object
points_sf <- st_as_sf(working.d, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
#use st_intersects to get a list of which cells on the global grid that these below to
intersections <- st_intersects(points_sf, global_grid_sf)
#a small subset of samples can match to multiple grid cells (falling right on the edge)
#we arbitrarily assign these samples instead to the first matched grid cell
#importantly, this is how samples were processed in building the effort matrrix
working.d$cell_id <- global_grid_sf$cell_id[sapply(intersections, function(x) if (length(x) > 0) x[1] else NA)]



###############################################################################
#crop our background effort matrix to a buffered region surrounding our points
###############################################################################
#convert observations to spatial points and project to metric (EPSG:8857)
points_sf <- st_as_sf(working.d, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
points_proj <- st_transform(points_sf, crs = 8857)

#create buffer area around all points (1000 km)
#because we are in EPSG:8857, dist is in meters
buffer_m <- 1000000 
point_buffers_proj <- st_buffer(points_proj, dist = buffer_m)

#union all those overlapping buffered regions into one single polygon
final_bg_poly_proj <- st_union(point_buffers_proj)


#convert centroids to sf, and project them to EPSG:8857 to match the polygon
grid_centroids_sf <- st_as_sf(grid_centroids, coords = c("center_lon", "center_lat"), crs = 4326)
grid_centroids_proj <- st_transform(grid_centroids_sf, crs = 8857)
grid_centroids_proj$cell_id <- rownames(grid_centroids)

#find which centroids fall inside the unioned buffers
intersections <- st_intersects(grid_centroids_proj, final_bg_poly_proj)
inside_buffer_indices <- which(lengths(intersections) > 0)

#extract the list of valid background cells for this specific species
valid_species_cells <- grid_centroids_proj$cell_id[inside_buffer_indices]

#subset the global background effort matrix and centroids to include only those in this map
working.bg_effort <- bg_effort[bg_effort$cell_id %in% valid_species_cells, ]
working.grid_centroids <- grid_centroids[rownames(grid_centroids) %in% valid_species_cells, ]

#there should be no new cellids in this species as opposed to the working background effort matrix
#as it should just be a susbet of the background matrix
#we see here that everything in the working species is in the background effort matrix
setdiff(unique(working.d$cell_id),unique(working.bg_effort$cell_id)) 



###############################################################################
###############################################################################
#to make things simple, we're going to remove all cells from the working.grid_centroids that are in the water, 
#using coordinate cleaner, unless they also have an observation from the background matrix in them, 
#in which case we'll keep them and just have to deal with them
#these are likely just coastal sites, whose centroid mean fell out of the expected 
flags <- clean_coordinates(x = data.frame(working.grid_centroids,species='A'), 
                           lon = "center_lon",lat = "center_lat",tests = c("seas"))
head(flags)

#these are the cells with observations in them
cells_w_obs <- row.names(working.grid_centroids[which(row.names(working.grid_centroids) %in% unique(working.bg_effort$cell_id)),])

#these are the cells that are coded as being on land
cells_on_land <- row.names(working.grid_centroids[which(flags$.sea == TRUE),])

#as occasionally, some of our observations will be coastal, the center of this grid may put them in the water
#we're still going to want to keep those cells
#so we combine these two lists of cell, and subset our grid centroids to include only these cells
valids_cells <- unique(c(cells_w_obs,cells_on_land))
working.grid_centroids <- working.grid_centroids[which(row.names(working.grid_centroids) %in% valids_cells),]

#we still want to initially focus on those that have observations
working.grid_centroids_obs <- working.grid_centroids[which(row.names(working.grid_centroids) %in% cells_w_obs),]



#######################################################################################
#download and process environmental data 
#######################################################################################
#set working directory to the environmental data that's been download, which we can read in
#we're extracting out all the environmental variables for each of cell and time period
#note however, that in this iteraction of the pipeline we actually don't include these 
#data in the models. We include them here though so they can be easily piped into others
#workflows
setwd("/scratch/kricks/GBIF_store/env_data")
altitude <- rast("global_elevation_0.9deg.tif") #load in altitude data
site_locs <- working.grid_centroids

#convert to terra vector
site_vect <- vect(site_locs, geom = c("center_lon", "center_lat"), crs = "EPSG:4326")


#writing a function that extracts values from a raster for a set of points, 
#and if there are NAs, it buffers the points and extracts the mean value from the buffer
#this is because there are points in our dataset that are coded as water, likely because 
#they're coastal and the grain of our analyses
env.raster <- altitude
buffer.width <- 100000 #can adjust
site.locations <- site_locs

initial.precip.1 <- terra::extract(altitude, site_vect)
site_covariates_all <- cbind(site_locs, initial.precip.1)
bad_sites <- site_covariates_all[is.na(site_covariates_all$precip_bin1_map), ]

VarID <- "altitude"

env_extract <- function(env.raster, site.locations, buffer.width, VarID) {
  #convert points to terra vector
  site_vect <- vect(site.locations, geom = c("center_lon", "center_lat"), crs = "EPSG:4326")
  
  #extract raster values
  env.values <- terra::extract(env.raster, site_vect)
  
  site_covariates_all <- data.frame(site.locations, env.values[,2])
  colnames(site_covariates_all)[3] <- VarID
  #identify points with NA values
  bad.points <- site_covariates_all[is.na(site_covariates_all[, 3]), ]
  
  if(nrow(bad.points) > 0) {
    #conver NA points to vector
    bad_vect <- vect(bad.points, geom = c("center_lon", "center_lat"), crs = "EPSG:4326")
    points_vect_buffered <- terra::buffer(bad_vect, width = buffer.width)
    
    #extract out using buffered polygons
    env_values_buffered <- terra::extract(
      env.raster, 
      points_vect_buffered, 
      fun = mean, 
      na.rm = TRUE
    )
    
    #update the main dataframe
    site_covariates_all[is.na(site_covariates_all [, 3]), 3] <- env_values_buffered[, 2]
  }
  return(site_covariates_all)
}


all.env <- data.frame(env_extract(altitude, site_locs, buffer.width, "Altitude"))

#sites that don't have altitude estimates i just set to 0, as these are likely coastal
#and at sea level
if(sum(is.na(all.env$Altitude))>0){
  all.env$Altitude[is.na(all.env$Altitude)] <- 0
}

all.env$AbsLat <- abs(all.env$center_lat)

#extract environmental covariate just for location with observations
all.env_observed <- all.env[which(row.names(all.env) %in% unique(rownames(working.grid_centroids_obs))),]
dim(all.env_observed)

######################################################################################
#build detection/nondetection matrices for each time period
#this is for both our species of interest, as well as overall effort
#######################################################################################
#standardize cell_id factors across the board based on the working grid
target_cells <- row.names(working.grid_centroids_obs)
working.bg_effort$cell_id <- factor(working.bg_effort$cell_id, levels = target_cells)
working.d$cell_id <- factor(working.d$cell_id, levels = target_cells)

#number of time bins
n_bins <- 6

#build aligned matrices for all bins
matrices_list <- lapply(1:n_bins, function(b) {
  #subset data for the current time Bin
  bg_sub <- working.bg_effort[working.bg_effort$Bin == b, ]
  d_sub  <- working.d[working.d$Bin == b, ]
  #define subbin for this temporal bin
  universal_subbins <- sort(unique(c(bg_sub$SubBin, d_sub$SubBin)))
  
  #apply these  levels to both subsets as factors
  bg_sub$SubBin <- factor(bg_sub$SubBin, levels = universal_subbins)
  d_sub$SubBin  <- factor(d_sub$SubBin, levels = universal_subbins)
  
  #create the background effort matrix
  effort_mat <- xtabs(N ~ cell_id + SubBin, data = bg_sub)
  
  #create species detection matrix
  species_mat <- table(d_sub$cell_id, d_sub$SubBin)
  species_mat[species_mat > 0] <- 1 #convert to presence/absence
  
  #return matrices
  return(list(effort = effort_mat, species = species_mat))
})

#transform list into your global environment variable names
for (b in 1:n_bins) {
  assign(paste0("counts_matrix_", b), matrices_list[[b]]$effort)
  assign(paste0("species_counts_matrix_", b), matrices_list[[b]]$species)
}

#sanity check :)
matrices_list[[1]]$effort[matrices_list[[1]]$species>0]
sum(matrices_list[[1]]$species[matrices_list[[1]]$effort==0])


######################################################################################
#spOccupancy requires arrays-we transform our data into that fomr here
#######################################################################################
#extract the effort and species matrices from the list
effort_list <- lapply(matrices_list, function(x) x$effort)
species_list <- lapply(matrices_list, function(x) x$species)

#convert 0 effort to NA in the effort and species matricies  
n_bins <- length(matrices_list)
for (i in 1:n_bins) {
  #identify where effort is 0
  zero_effort <- effort_list[[i]] == 0
  
  #replaces these indices with NA in both species and effort matrices
  species_list[[i]][zero_effort] <- NA
  effort_list[[i]][zero_effort] <- NA
}

#find the maximum number of columns (subbins) across all time bins
max_cols <- max(sapply(species_list, ncol))

#creating a helper function to pad a matrix with NAs up to max_cols
#this may feel silly, but we need the same level of subbins for the occupancy models
#so for example here, the earliest time bin (1900-1965) has more subbins then the other 
#time bins (all subbins are 5 years). but when input into spoccupancy there needs to be 
#equal number of subbins so we just add a bunch of NAs
pad_with_na <- function(mat, max_c) {
  current_c <- ncol(mat)
  if (current_c < max_c) {
    #create a matrix of NAs for the missing columns
    na_padding <- matrix(NA, nrow = nrow(mat), ncol = (max_c - current_c))
    mat <- cbind(mat, na_padding)
  }
  return(mat)
}

#apply the padding function to all matrices in the lists
c(dim(species_list[[1]])[2],dim(species_list[[2]])[2],dim(species_list[[3]])[2],
  dim(species_list[[4]])[2],dim(species_list[[5]])[2],dim(species_list[[6]])[2])
species_list_padded <- lapply(species_list, pad_with_na, max_c = max_cols)
effort_list_padded  <- lapply(effort_list, pad_with_na, max_c = max_cols)
c(dim(species_list_padded[[1]])[2],dim(species_list_padded[[2]])[2],dim(species_list_padded[[3]])[2],
  dim(species_list_padded[[4]])[2],dim(species_list_padded[[5]])[2],dim(species_list_padded[[6]])[2])

#convert these arrays and reshape dimensions to be in the right order
species_array <- simplify2array(species_list_padded)
species_array <- aperm(species_array, c(1, 3, 2))
species_array[,2,]

count_array <- simplify2array(effort_list_padded)
count_array <- aperm(count_array, c(1, 3, 2))
count_array[,2,]

#check to make sure it looks right :)
dim(species_array)
dim(count_array)

dim(count_array[,2,])

#######################################################################################
#merge data into one big structure
#######################################################################################
#creating a matrix of observed occupancy over time

#here are matrices for altitude and absolute latitude. These represent environmental correlates
#that we could include in the models. We ultimatley didn't include these here, but we include 
#code here for others if need 
alt_matrix <- all.env_observed[,rep("Altitude",6)]
lat_matrix <- all.env_observed[,rep("AbsLat",6)]


#the dimensions should be: sites x time bins
#years needs to be a matrix too 
#create matrix from scatch with same dimensions as alt and lat
site.length <- nrow(alt_matrix)
years_matrix <- matrix(rep(1:6, each = site.length), nrow = site.length, ncol = 6)

occ.list.dynamic <- list(alt = alt_matrix,     
                         lat = lat_matrix, 
                         years = years_matrix,   
                         ID = rownames(all.env_observed))
dim(occ.list.dynamic$temp)
dim(occ.list.dynamic$precip)
dim(occ.list.dynamic$years)
length(occ.list.dynamic$ID)


#######################################################################################
#create custom projection system, that scale to species of interest
#everything is centered on the centroid of this species, with distance in km
#relative to this centroid
#######################################################################################
coords_sf <- st_as_sf(data.frame(all.env_observed), coords = c("center_lon", "center_lat"), crs = 4326)

#finde the geographic center of this specific subset's bounding box
#(st_union with s2 geometry enabled safely handles dateline wrapping)
centroid <- st_centroid(st_union(coords_sf))
center_coords <- st_coordinates(centroid)

#build a custom LAEA projection centered on this species' data
custom_crs <- paste0("+proj=laea +lat_0=", center_coords[,"Y"], 
                     " +lon_0=", center_coords[,"X"], 
                     " +datum=WGS84 +units=m +no_defs")

#transform coordinates 
coords_proj <- st_transform(coords_sf, crs = custom_crs)

#convert to km for spOccupancy scaling
metric_coords <- st_coordinates(coords_proj) / 1000 


#################################################################
#scale all variables
#################################################################
#scaling
scale_params <- data.frame(Variable = c("alt", "lat","years"),
                           Mean = c(mean(as.matrix(occ.list.dynamic$alt), na.rm = TRUE), 
                                    mean(as.matrix(occ.list.dynamic$lat), na.rm = TRUE),
                                    mean(as.matrix(occ.list.dynamic$years), na.rm = TRUE)),
                           
                           SD = c(sd(as.matrix(occ.list.dynamic$alt), na.rm = TRUE), 
                                  sd(as.matrix(occ.list.dynamic$lat), na.rm = TRUE),
                                  sd(as.matrix(occ.list.dynamic$years), na.rm = TRUE)))

occ.list.dynamic$alt <- (occ.list.dynamic$alt - scale_params$Mean[1]) / scale_params$SD[1]
occ.list.dynamic$lat <- (occ.list.dynamic$lat - scale_params$Mean[2]) / scale_params$SD[2]
occ.list.dynamic$years <- (occ.list.dynamic$years - scale_params$Mean[3]) / scale_params$SD[3]  


log_transform_count_array <- log(count_array)
mean_effort <- mean(log_transform_count_array, na.rm = TRUE)
sd_effort <- sd(log_transform_count_array, na.rm = TRUE)
scaled_count_array <- (log_transform_count_array - mean_effort) / sd_effort


merged.data.dynamic <- list(occ.covs = occ.list.dynamic,
                            y = species_array,
                            det.covs = list(effort = scaled_count_array),
                            coords = as.matrix(all.env_observed[, c(1, 2)]))

merged.data.dynamic$coords <- metric_coords



#######################################################################################
#model building
#######################################################################################
library(splines)

#this is the occupancy formula
my.occ.formula <- ~ years  
#this is the detection formula. we add on an exponential here to allow for saturation of GBIF effort
#in predicting detection 
my.det.formula <- ~ effort + I(effort^2)

#calculate min and max distances to set realistic priors for spatial decay (phi)
#this is standard approach in the spoccupancy pipeline
dist.mat <- dist(merged.data.dynamic$coords)
min.dist <- min(dist.mat[dist.mat > 0]) 
max.dist <- max(dist.mat)               

z.inits <- apply(merged.data.dynamic$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))

#spatial parameters inits
my.inits <- list(beta = 0, 
                 alpha = 0, 
                 z = z.inits, 
                 sigma.sq = 1, #spatial variance
                 phi = 3 / mean(dist.mat), #spatial decay
                 sigma.sq.t = 1.5, 
                 rho = 0.2) 

#spatial parameters priors
my.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                  alpha.normal = list(mean = 0, var = 2.72), 
                  sigma.sq.ig = c(2, 1),#spatial variance prior
                  phi.unif = c(3 / max.dist, 3 / min.dist),#spatial decay prior
                  sigma.sq.t.ig = c(2, 0.5), 
                  rho.unif = c(-1, 1))




#exponential model with 10 nearest neighbors and ar1 = T
cov.model <- 'exponential'
n.neighbors <- 10
ar1 <- TRUE

#how many iterations to run through 
n.batch <- 3720
batch.length <- 25
n.burn <- 30000 
n.thin <- 15

#run model
out.sp <- stPGOcc(occ.formula = my.occ.formula, 
                  det.formula = my.det.formula, 
                  data = merged.data.dynamic, 
                  inits = my.inits, 
                  priors = my.priors, 
                  cov.model = cov.model,#spatial
                  n.neighbors = n.neighbors,#spatial (NNGP)
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  verbose = TRUE, 
                  ar1 = ar1,
                  n.report = 200,
                  n.burn = n.burn, 
                  n.thin = n.thin,
                  n.omp.threads = 6,
                  n.chains = 3) 


#there's a lot of data in these model output that. we don't need it all, and with ~6000 species we won't have the space to save all of these model
#we clear these details out the of the model object before saving the model
out.sp$z.samples    <- NULL
out.sp$psi.samples  <- NULL
out.sp$like.samples <- NULL

out.sp$p.samples     <- NULL
out.sp$y.rep.samples <- NULL

#save model output and a summary of the output
save_path_rds <- paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/Output/Pure_SPAT_rds/model_", species.id, ".rds")
saveRDS(out.sp, file = save_path_rds)

save_path_summary <- paste0("/scratch/kricks/GBIF_store/processed_Rhizobia/Output/Pure_SPAT_model/model_", species.id, ".txt")
capture.output(summary(out.sp), file = save_path_summary)

