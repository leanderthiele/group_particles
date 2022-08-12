#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=5G
#SBATCH -t 02:00:00

date

module load anaconda3/2021.11
module load hdf5/gcc/1.10.6

python -u CAMELS_meanT.py SIMBA LH 028 1

date
