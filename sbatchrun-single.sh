#!/bin/bash

# Job name
#SBATCH --job-name=hello-slurm
# Number of Nodes
#SBATCH --nodes=1
# Set partition
#SBATCH --partition=single
# Number of CPU-cores per task
#SBATCH --cpus-per-task=1

echo "single tree"
echo $1
/work/kulagins/scheduler-a2/$1


