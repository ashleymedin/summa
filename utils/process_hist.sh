#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=0-01:00
#SBATCH --job-name=HIST
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  geo-stack/2023a

python hist_per_GRU.py rmse
python hist_per_GRU.py maxe
python hist_per_GRU.py kgem
