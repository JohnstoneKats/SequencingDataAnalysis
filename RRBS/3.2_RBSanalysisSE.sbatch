#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=6:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

module load bowtie2

module load bismark

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`

echo $s
echo $x

bismark --multicore 4  /data/reference/indexes/mouse/mm10/bowtie2/ $1

