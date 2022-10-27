#!/bin/sh
#SBATCH --partition=AllNodes
#SBATCH --job-name=sko11
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --output=output-npart_ko1

##ls /local/storage
srun singularity exec  /home/dongwj/root.sif ./command_1.sh $SLURM_ARRAY_TASK_ID
