#!/bin/bash                       # Use the bash shell interpreter
#SBATCH -J run_singler            # Give the job a name
#SBATCH --time=2:00:00            # Request 2 hours
#SBATCH -n 1                      # Request 1 core
#SBATCH -N 1                      # Request 1 node
#SBATCH --mem=10Gb                # Request 10 Gb
#SBATCH --output=%j.out           # Write the job output to a file prefixed by the job number
#SBATCH --error=%j.err            # Write the job error to a file prefixed by the job number
 
module purge                      # Remove loaded modules
module load R/4.0.0               # Load R module

Rscript --no-save 06_singler_cell.R  # Run the script 