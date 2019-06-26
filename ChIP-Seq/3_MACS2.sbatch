#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short



module load macs

s=$(basename $1)

echo $s
echo $1 


macs2 callpeak -g mm -f BAM -t $1 -c Sample.sorted.rmdup.bam --cutoff-analysis --outdir ${1}_macs -n ${s}peaks

# or for broad peaks use
# macs2 callpeak -g mm -f BAM -t $1 -c Sample.sorted.rmdup.bam --broad -q 0.01 --cutoff-analysis --outdir ${1}BROADq01_macs2 -n ${s}BROADq01_peaks

