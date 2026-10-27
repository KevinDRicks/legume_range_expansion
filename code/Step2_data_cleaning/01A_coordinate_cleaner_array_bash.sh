#!/bin/bash
#SBATCH --job-name=occ_models
#SBATCH --output=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/clean_%A_%a.out
#SBATCH --error=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/clean_%A_%a.err
#SBATCH --time=02:00:00                 # 2 hours should safely cover a set of 300 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4               
#SBATCH --mem=20G                       # Adjust this based on your peak memory usage
#SBATCH --array=1-30%30          # 30 total tasks, run a maximum of 30 simultaneously

#this script runs an R script, 01B_coordinate_cleaner_array_format.R, to iteratively run coordinatecleaner
#on all of our data
#note, this is setup to run as a batched array, to run in a timely manner



#load r module
module load r

# Pass the Slurm array task ID directly into the R script
Rscript /home/kricks/codebases/range_expansion/publication_scripts/data_cleaning/01B_coordinate_cleaner_array_format.R $SLURM_ARRAY_TASK_ID