#!/bin/bash

#SBATCH --mem=8g
#SBATCH --cpus-per-task 11
#SBATCH -t 02:00:00

module load r/4.4.0

Rscript --max-connections=512 temp_analyses/run_wtx.R "$1" "$2" "$3" "$4"