#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=2:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod_short

 
module load fastqc

s=$(basename $1)
d=$(dirname $1)

echo ${1}

### can input as a list of files

fastqc -o ${d} -f fastq --noextract -t 8 ${1}

