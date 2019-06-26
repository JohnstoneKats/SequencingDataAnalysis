#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=2:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
# SBATCH --partition=prod_short


module load bismark

echo ${1}

#input file is the output from bismark alignment

bismark_methylation_extractor --bedgraph ${1}

