#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=2:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short

module load homer

s=$(basename $1)

echo $1

findPeaks ${1}tagdir/ -i 180609_ChIPAc/LK-input_S16.samtagdir -style super -o ${1}tagdir/${s}Super.txt -superSlope -1000 -L 2 -typical ${1}normalpeaks.txt




