#!/bin/bash --login
#SBATCH --job-name=TEST1
#SBATCH --output=/net/work/wimmers/operations_oppenheimer/logs/log.txt
#SBATCH --time=150:00:00
#SBATCH --mem=40g
#SBATCH --partition=smp
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT,END
#SBATCH --mail-user=zoerner@tu-berlin.de

export JULIA_DEPOT_PATH="/net/work/wimmers/juliapkg"

module add julia/1.7.3

julia src/general_v2.jl



