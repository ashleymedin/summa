#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=1-00:00
#SBATCH --job-name=CONCAT
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  geo-stack/2023a

python concat_groups_split_summa.py sundials_1en8
