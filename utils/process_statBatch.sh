#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=0-05:00
#SBATCH --job-name=STATB
#SBATCH --mail-user=gwu479@usask.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rpp-kshook
#SBATCH --output=/home/avanb/TestScripts/output/slurm-%A_%a.out
echo "$SLURM_ARRAY_TASK_ID"

# ----------------------------------------------------------------------------------------------
# RUN WITH:
# sbatch --array1-[number of jobs] [script name]
# sbatch --array=1-200 process_statBatch.sh
# sbatch --array=201-201 process_statBatch.sh
# ----------------------------------------------------------------------------------------------

module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  geo-stack/2023a

python timeseries_to_statistics.py sundials_1en6 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be1 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be32 $SLURM_ARRAY_TASK_ID 200
python timeseries_to_statistics.py be16 $SLURM_ARRAY_TASK_ID 200
