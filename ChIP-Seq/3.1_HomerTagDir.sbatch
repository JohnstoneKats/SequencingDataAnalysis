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
module load samtools

s=$(basename $1)

echo $1


samtools view -h $1 > ${1}.sam
makeTagDirectory ${s}tagdir/ ${1}.sam -format sam
findPeaks ${1}tagdir -style factor -o auto -i inputtagdir/
pos2bed.pl -o ${1}tagdir/peaks.bed ${1}tagdir/peaks.txt
awk '{print "chr"$0}' ${1}tagdir/peaks.bed > ${1}tagdir/peaks.bed.chr.bed
annotatePeaks.pl ${1}tagdir/peaks.bed.chr.bed mm10 > ${1}tagdir/peaksannotated.txt


