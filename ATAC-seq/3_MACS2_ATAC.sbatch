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
x=`echo $s | cut -d "." -f 1`

#echo $s
echo $x
echo $1


macs2 callpeak -g mm -f BAM  --call-summits --nomodel -t ${1} --outdir ${x}_macs2 -n ${x}_peaks
