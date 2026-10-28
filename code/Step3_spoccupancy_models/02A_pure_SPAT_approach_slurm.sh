#!/bin/bash
#SBATCH --job-name=Pure_SPAT_occ_models
#SBATCH --output=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/Pure_SPAT_logs/occ_model_%A_%a.out
#SBATCH --error=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/Pure_SPAT_logs/occ_model_%A_%a.err
#SBATCH --time=16:00:00                 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6               
#SBATCH --mem=30G                       
#SBATCH --array=1-6000%300          #each species is a single task, 300 going at once

#this script runs an R script, 02B_pure_SPAT_approach.R, to iteratively build occupancy models
#on each species
#note, this is setup to run as a batched array, to run in a timely manner
#note, while 16 hours and 30 G works for the vast majority of species, for a small fraction 
#of the species, this will not be enough memory or node time
#these will have to go back and be rerun with more memory or time

#load r module
module load r

#passes the Slurm array task ID directly into the R script
Rscript /home/kricks/codebases/range_expansion/publication_scripts/spoccupancy_models/02B_pure_SPAT_approach.R $SLURM_ARRAY_TASK_ID
