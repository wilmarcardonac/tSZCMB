#!/bin/sh
#SBATCH --cpus-per-task=12
#SBATCH --job-name=tSZ-CMB
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out

srun ./tSZ