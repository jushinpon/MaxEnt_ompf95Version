#!/bin/sh
#SBATCH --output=maxent.out
#SBATCH --job-name=maxent
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --ntasks-per-node=12
rm -rf ./output
mkdir output
threads=`lscpu|grep "^CPU(s):" | sed 's/^CPU(s): *//g'`
export OMP_NUM_THREADS=$threads
/opt/maxent/maxent.x
