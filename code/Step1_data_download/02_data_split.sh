#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-freder19
#SBATCH --job-name=my-job
#SBATCH --output=/scratch/kricks/GBIF_store/rhizobia_subset.out
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=750G

#the entirity of the gbif dataset was download as a single csv file
#this is entirely too large to deal with at any point
#this script runs a fairly simple r script that efficiently opens 
#this file and then splits it by the separate gbif species keys
#note: this is about the computational limit for our cluster
#if expanding to a larger dataset, a more computaitonally efficient way 
#of splitting these data is necessary



cd /home/kricks/scratch/GBIF_store/processed_Rhizobia/
rm -rf subset_data
mkdir subset_data

module load r


Rscript - <<EOF
library(data.table)
library(bit64)
temp.data <- fread("gbif_extracted/0011505-260806074905277.csv", quote = "")
invisible(temp.data[, fwrite(data.table(speciesKey = .BY$speciesKey, .SD),
                              paste0("subset_data/", .BY$speciesKey, ".csv")),
                    by = speciesKey])


EOF