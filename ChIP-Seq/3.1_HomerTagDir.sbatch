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

#input as sam file, so first convert processed bam to sam
samtools view -h $1 > ${1}.sam
makeTagDirectory ${s}tagdir/ ${1}.sam -format sam

#need to specify tag directory of input 
findPeaks ${1}tagdir -style factor -o auto -i inputtagdir/
pos2bed.pl -o ${1}tagdir/peaks.bed ${1}tagdir/peaks.txt

#add "chr" if using reference files without chr
awk '{print "chr"$0}' ${1}tagdir/peaks.bed > ${1}tagdir/peaks.chr.bed

#annotate peaks with closest gene
annotatePeaks.pl ${1}tagdir/peaks.chr.bed mm10 -annStats ${1}tagdir/annstats.txt  > ${1}tagdir/peaksannotated.txt


