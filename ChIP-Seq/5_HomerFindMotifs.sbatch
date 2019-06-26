#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=4:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod

module load homer

#input is a peaks file
findMotifsGenome.pl $1 mm10 ${1}motif 
