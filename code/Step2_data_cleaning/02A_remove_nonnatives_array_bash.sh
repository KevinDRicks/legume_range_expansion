#!/bin/bash
#SBATCH --job-name=nonnative_clean
#SBATCH --output=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/nonnative_clean_%A_%a.out
#SBATCH --error=/home/kricks/codebases/range_expansion/Nitrogen_fixation/slurm_submission_models/nonnative_clean_%A_%a.err
#SBATCH --time=01:00:00                 # 1 hour should safely cover a set of 10 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G                       # Adjust this based on your peak memory usage
#SBATCH --array=1-900%300          # 900 total tasks, run a maximum of 10- simultaneously
#SBATCH --cpus-per-task=1         # R is single-threaded


#this script runs an R script, 02B_remove_nonnative_array.R, to prune out invasive occurences from the datasets
#note, this is setup to run as a batched array, to run in a timely manner


#load r module
module load r

# Pass the Slurm array task ID directly into the R script
Rscript /home/kricks/codebases/range_expansion/publication_scripts/data_cleaning/02B_remove_nonnative_array.R $SLURM_ARRAY_TASK_ID