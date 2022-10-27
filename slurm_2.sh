#!/bin/sh
#SBATCH --partition=AllNodes
#SBATCH --job-name=sko12
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=output-npart_ko2

##ls /local/storage
srun singularity exec  /home/dongwj/root.sif ./command_2.sh $SLURM_ARRAY_TASK_ID
