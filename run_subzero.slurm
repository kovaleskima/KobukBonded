#!/bin/bash

#SBATCH --job-name=Subzero_Kobuk_Bonds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=makov@uw.edu

#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=50G
#SBATCH --gpus=0
#SBATCH --time=07-00:00:00 

#SBATCH --chdir=/mmfs1/home/makov/KobukBonded
#SBATCH --export=all
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.out

# Modules to use (optional).
module load matlab

# Your programs to run.
matlab -nodisplay -r "SubzeroClean;quit;"
