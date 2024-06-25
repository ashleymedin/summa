#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00
#SBATCH --job-name=DIFF
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out

module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  geo-stack/2023a

python check_bit_4_bit_withTol.py /home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-sundials_be_noJac/ /home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/summa-be/ 0.1 > out.txt
