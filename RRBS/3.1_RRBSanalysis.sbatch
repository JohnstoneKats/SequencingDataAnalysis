#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod

module load bowtie2
module load bismark

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`

echo $s
echo $x

bismark --unmapped --multicore 4 /data/reference/indexes/mouse/mm10/bowtie2/ -1 ${t}/${x}_R1_001_val_1.fq.gz -2 ${t}/${x}_R2_001_val_2.fq.gz 

