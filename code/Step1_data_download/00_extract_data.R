#this script inputs data from Werner et al (2014) and Kates et al (2024), and collects 
#all legumes with a known nitrogen fixation status. It then extracts a GBIF ID key and 
#submits it for download download. 

#we have to do this in three parts: 1) loading in werner data and aligning to GBIF IDs,
#2) loading in Kates data and aligning to GBIF IDs, and 3) reconciling differences between
#the two datasets


#load in packages
library(rgbif)
library(data.table)

######################################################################################################
#load in and align species taxonomy from Werner et al (2014)
######################################################################################################
setwd('/home/kricks/codebases/range_expansion/publication_scripts/data_download')
data.nf <- read.csv('updated_werner.csv')
data.nf <- data.nf[data.nf$Legume=='Yes',]
table(data.nf$NFC)
table(data.nf$Data_fixing)

#extract the GBIF backbone
matches_df <- name_backbone_checklist(name = data.nf$Species,kingdom='Plantae')
#matches_df <- name_backbone_checklist(name = data.nf$Species,kingdom='Plantae', bucket_size = 100,sleep = 3)#sometimes the above command fails and requires the alternative slower approach

#not all of these match back to a species
table(matches_df$rank)

#extract only those that can be called to species
subset.nf <- data.nf[matches_df$rank=='SPECIES',]
subset.match <- matches_df[matches_df$rank=='SPECIES',]

sum(matches_df$rank=='SPECIES')
table(matches_df$rank)

#a couple may be misspelled but likely okay
table(subset.match$matchType)

#now add on the key for those species
subset.nf$Key <- subset.match$speciesKey

#this is how many unique species there actually are
length(sort(table(subset.match$speciesKey),decreasing = T))

duplicate.taxa <- sort(table(subset.match$speciesKey),decreasing = T)
duplicate.taxa <- duplicate.taxa[duplicate.taxa>1]#659 with duplicate names 

#example of some of the duplicates
subset.nf[subset.nf$Key=='2951499',]
subset.nf[subset.nf$Key=='2951499',]$Species

#do these duplicates have the same assigned fixation status?
duplicate.status <- NULL
for(i in 1:length(duplicate.taxa)){
  st.vec <- table(subset.nf[subset.nf$Key==names(duplicate.taxa)[i],]$Data_fixing)
  yes.st <- st.vec[names(st.vec)=='Yes']
  no.st <- st.vec[names(st.vec)=='No']
  if(length(yes.st)==0){
    yes.st <- 0
  }
  if(length(no.st)==0){
    no.st <- 0
  }
  duplicate.status <- rbind(duplicate.status,
                            c(names(duplicate.taxa)[i],yes.st,no.st))
}
duplicate.status <- data.frame(duplicate.status)
colnames(duplicate.status) <- c('Key','Yes','No')
duplicate.status[,2] <- as.numeric(duplicate.status[,2])
duplicate.status[,3] <- as.numeric(duplicate.status[,3])
dup.remove <- duplicate.status[duplicate.status$Yes>0 & duplicate.status$No>0,]$Key
#this is a small number that should be removed as its unclear their status

dim(subset.nf)
subset.nf <- subset.nf[!(subset.nf$Key %in% dup.remove),]
dim(subset.nf)

unq.sp <- unique(subset.nf$Key)
subset.nf_new <- NULL
for(i in 1:length(unq.sp)){
  subset.nf_new <- rbind(subset.nf_new,
                         subset.nf[subset.nf$Key==unq.sp[i],][1,])#extract the first row. Other rows may be need
}
dim(subset.nf_new)


######################################################################################################
#load in and align species taxonomy from Kates et al (2024)
######################################################################################################
data.nf <- read.csv('kates_dat.csv')
data.nf <- data.nf[data.nf$Family=="Leguminosae",]

#two different sorts of names
sum(data.nf$scientificName_Voucher!=data.nf$Taxon.in.Tree)

#Between the voucher name and the tree name, there are multiple mismatches in names (1353)
#I split them into a mismatch and match sets, and will work with them separately, 
#until i can rejoin them
mismatch.df <- data.nf[data.nf$scientificName_Voucher!=data.nf$Taxon.in.Tree,]
matching.df <- data.nf[data.nf$scientificName_Voucher==data.nf$Taxon.in.Tree,]

#for those that match, it doesn't matter which label i use
matching.key <- name_backbone_checklist(name = matching.df$Taxon.in.Tree,kingdom='Plantae')
matching.key.df <- data.frame(matching.key)

#turn this into a form that can be use down the line
matching.key.df <- data.frame(Key=matching.key.df$speciesKey,
                              Species=matching.df$Taxon.in.Tree,
                              State=matching.df$State)
#remove NA calls from keys, of which there is a small number (15. )
matching.key.df <- matching.key.df[!is.na(matching.key.df$Key),]

#Now we interrogate the mismatches 
#I first evaluate if these mismatched names produce different GBIF calls
#extract species keys from mismatching names
mismatch.tt<- name_backbone_checklist(name = mismatch.df$Taxon.in.Tree,kingdom='Plantae')
mismatch.vch <- name_backbone_checklist(name = mismatch.df$scientificName_Voucher,kingdom='Plantae')
mismatch.tt.df <- data.frame(mismatch.tt)
mismatch.df.vch <- data.frame(mismatch.vch)

#create workable df for analysis
combined.species.key <- data.frame(tt_key=mismatch.tt.df$speciesKey,
                                   vch_key=mismatch.df.vch$speciesKey,
                                   tt_label=mismatch.df$Taxon.in.Tree,
                                   vch_label=mismatch.df$scientificName_Voucher,
                                   State=mismatch.df$State)

#Between all of these, every species received at least one speciesKey between the two potential names
combined.species.key <- combined.species.key[!(is.na(combined.species.key$tt_key) & is.na(combined.species.key$vch_key)),]
sum(is.na(combined.species.key$tt_key))#the tree had 6 missing values
sum(is.na(combined.species.key$vch_key))#the vouchers had 91 missing values

#I'll now disregard samples that only have 1 successful assignment
#those with only 1 assignment can be joined with the rest of the analyses
combined.species.key.one_success <- combined.species.key[is.na(combined.species.key$tt_key) | 
                                                           is.na(combined.species.key$vch_key),]
combined.species.key.two_success <- combined.species.key[!(is.na(combined.species.key$tt_key) | 
                                                             is.na(combined.species.key$vch_key)),]

#of these, 925 result in the same call, and 331 have different calls
sum(combined.species.key.two_success$tt_key==combined.species.key.two_success$vch_key)
sum(combined.species.key.two_success$tt_key!=combined.species.key.two_success$vch_key)

#successful matches, it doesn't matter what label we use, we'll just include it
matching_key <- combined.species.key.two_success[combined.species.key.two_success$tt_key==combined.species.key.two_success$vch_key,]
mismatching_key <- combined.species.key.two_success[combined.species.key.two_success$tt_key!=combined.species.key.two_success$vch_key,]


#I'll dig a bit into those with different calls, but maybe not too much
#i've go through these by hand, and about half of them are synonyms. This can be seen in their
#gbif ranges, which are completley overlapping. This just means data could be combined together
#however, some are clearly not synonyms, with completely separate ranges. Unsure what to do with those 
#for now, will not use any of these


#clean up different subsets so they can be pooled.
#renamed species with only one successful label from the mismatches
d1 <- combined.species.key.one_success[!is.na(combined.species.key.one_success$tt_key),c(1,3,5)]
d2 <- combined.species.key.one_success[!is.na(combined.species.key.one_success$vch_key),c(2,5,5)]
colnames(d1) <- c('Key','Species','State');colnames(d2) <- c('Key','Species','State')
matching_key <- matching_key[,c(1,3,5)];colnames(matching_key) <- c('Key','Species','State')
#mismatched_new <- mismatching_key[,c("tt_key","tt_label","State")];colnames(mismatched_new) <- c('Key','Species','State')
  
#merge ids  
cleaned_names <- rbind(matching.key.df,
                       d1,d2,matching_key)

#remove unknown state
cleaned_names <- cleaned_names[cleaned_names$State!="?",]


#~7000 unique names
length(unique(cleaned_names$Key))
table((cleaned_names$State))

#a small number are repeats
#I'm going to check that those repeats have the same estimated state

duplicate.taxa <- sort(table(cleaned_names$Key),decreasing = T)
duplicate.taxa <- duplicate.taxa[duplicate.taxa>1]#343 with duplicate names 

#do these duplicates have the same assigned fixation status?
duplicate.status <- NULL
for(i in 1:length(duplicate.taxa)){
  st.vec <- table(cleaned_names[cleaned_names$Key==names(duplicate.taxa)[i],]$State)
  yes.st <- st.vec[names(st.vec)=='Present']
  no.st <- st.vec[names(st.vec)=='Absent']
  if(length(yes.st)==0){
    yes.st <- 0
  }
  if(length(no.st)==0){
    no.st <- 0
  }
  duplicate.status <- rbind(duplicate.status,
                            c(names(duplicate.taxa)[i],yes.st,no.st))
}
duplicate.status <- data.frame(duplicate.status)
colnames(duplicate.status) <- c('Key','Present','Absent')
duplicate.status[,2] <- as.numeric(duplicate.status[,2]);duplicate.status[,3] <- as.numeric(duplicate.status[,3])
duplicate.status[duplicate.status$Present>0 & duplicate.status$Absent>0,]
#All of them had the exact same esimated state


#i'm going to condense all the species into a single table
unq.sp <- unique(cleaned_names$Key)
new_clean <- NULL
for(i in 1:length(unq.sp)){
  new_clean <- rbind(new_clean,
                     cleaned_names[cleaned_names$Key==unq.sp[i],][1,])#extract the first row. Other rows may be need
}

table((cleaned_names$State))
table((new_clean$State))



######################################################################################################
#reconcile data between the werner and kates datasets 
######################################################################################################
#only 50% overalp with the werner dataset,
int.sp <- (intersect(subset.nf_new$Key,new_clean$Key))
length(int.sp)

overlap.df <- NULL
for(i in 1:length(int.sp)){
  overlap.df <- rbind(overlap.df,c(subset.nf_new[subset.nf_new$Key==int.sp[i],]$Data_fixing,
                                   new_clean[new_clean$Key==int.sp[i],]$State))
}


table(overlap.df[overlap.df[,1]=='No',2])
table(overlap.df[overlap.df[,1]=='Yes',2])
#small number of differences between datasets. These can be flagged for later 

colnames(new_clean) <- c('Key','Species','NF_Status_Kates')
colnames(subset.nf_new)[c(1,6,12)] <- c('Species','NF_Status_Werner','Key')
combine_datasets <- rbind(data.frame(new_clean,NF_Status_Werner=NA),
                          data.frame(subset.nf_new[,c(12,1)],NF_Status_Kates=NA,NF_Status_Werner=subset.nf_new[,c(6)]))

combine_datasets[combine_datasets$NF_Status_Kates=='Absent' & 
                   !is.na(combine_datasets$NF_Status_Kates=='Absent'),]$NF_Status_Kates <- 'Nonfixer'
combine_datasets[combine_datasets$NF_Status_Kates=='Present' & 
                   !is.na(combine_datasets$NF_Status_Kates=='Present'),]$NF_Status_Kates <- 'Fixer'

combine_datasets[combine_datasets$NF_Status_Werner=='No' & 
                   !is.na(combine_datasets$NF_Status_Werner=='No'),]$NF_Status_Werner <- 'Nonfixer'
combine_datasets[combine_datasets$NF_Status_Werner=='Yes' & 
                   !is.na(combine_datasets$NF_Status_Werner=='Yes'),]$NF_Status_Werner <- 'Fixer'

combine_datasets$NF_Status_Estimate <- NA

unique_key <- unique(combine_datasets$Key)

combine_thin <- NULL
for(i in 1:length(unique_key)){
  temp.d <- combine_datasets[combine_datasets$Key==unique_key[i],]
  if(nrow(temp.d)==1){
    est.vec <- c(temp.d$NF_Status_Kates,temp.d$NF_Status_Werner)
    temp.d$NF_Status_Estimate <- est.vec[!is.na(est.vec)]
    combine_thin <- rbind(combine_thin,temp.d)
  }
  if(nrow(temp.d)==2){
    est.vec <- c(temp.d$NF_Status_Kates,temp.d$NF_Status_Werner)
    est.vec <- est.vec[!is.na(est.vec)]
    if(est.vec[1]==est.vec[2]){
      temp.d$NF_Status_Kates <- est.vec
      temp.d$NF_Status_Werner <- est.vec
      temp.d$NF_Status_Estimate <- est.vec
      combine_thin <- rbind(combine_thin,temp.d[1,])
    }
    
    if(est.vec[1]!=est.vec[2]){
      temp.d[1,c(3,4)] <- est.vec
      temp.d$NF_Status_Estimate <- 'Mixed'
      combine_thin <- rbind(combine_thin,temp.d[1,])
    }
  }
  print(i)
}

table(combine_thin$NF_Status_Estimate)


max(table(combine_thin$Key))
head(combine_thin)


#extract all of the gbif names as well
keys <- combine_thin$Key
names_list <- lapply(keys, function(k) {
  name_usage(key = k)$data$canonicalName
})

# Convert to a simple character vector
combine_thin$gbif_name <- unlist(names_list)

#output a description of all the species identified
write.csv(combine_thin,"gbif_key_merge.csv",row.names = F)
combine_thin <- read.csv("gbif_key_merge.csv")
nrow(combine_thin)



######################################################################################################
#inputting gbif download
######################################################################################################

#we want the species key, because this is the currently accepted name, and can pull in other records
#that are synonyms that have been combined. Useage key might be just for a single variety, or subspecies
#or previously accepted synonym
#occ_download(pred_in("SPECIES_KEY",as.list(combine_thin$Key)),format = "SIMPLE_CSV",pred("hasCoordinate", TRUE),
#             user='kricks',pwd='#htv_i4ZebuFqv2',email='kevin.d.ricks@gmail.com')
occ_download_wait('0011505-260806074905277')
#
#occ_download_import(as.download("0008695-260226173443078.zip"))





