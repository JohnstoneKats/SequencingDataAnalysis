#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=2:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod_short


module load macs


s=$(basename $1)

echo $s
echo $1


macs2 callpeak -g mm -f BAMPE -t $1 -c $2 --nomodel --gsize 2.82e8 --extsize 200 --outdir ${s}_macs2 -n ${s}_normpeaks
