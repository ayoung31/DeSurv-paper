#!/bin/bash
#SBATCH --job-name=desurv_pipeline
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

# Wait until system load drops (non-Slurm R processes finish)
echo "Waiting for system load to drop below 4..."
while [ $(awk '{print int($1)}' /proc/loadavg) -ge 4 ]; do
    echo "$(date): Load is $(cat /proc/loadavg | cut -d' ' -f1), waiting..."
    sleep 60
done

echo "$(date): Load dropped, starting pipeline..."
cd /home/naimrashid/Downloads/DeSurv-paper
Rscript -e 'targets::tar_make()'
