#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --partition=debug
#SBATCH --exclusive
#SBATCH --output=outputfile.txt
#SBATCH --error=outputfile.txt

module load matlab/2019b
matlab -r SprayDet_Matlab
