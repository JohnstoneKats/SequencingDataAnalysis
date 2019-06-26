#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=6:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod

module load fastqc

echo ${1}

### can input as a list of files

fastqc -o ${d} -f fastq --noextract -t 8 ${1}

