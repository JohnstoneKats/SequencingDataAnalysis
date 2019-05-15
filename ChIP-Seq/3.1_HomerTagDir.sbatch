#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### If experimenting: use srun --pty /bin/bash

 
module load homer
module load samtools

s=$(basename $1)
#x=`echo $s | cut -d "_" -f 1-2`

#echo $s
#echo $x
echo $1

### uses homer to make tagdir, call peaks, then annotate peaks ### (also adds chr to bed file if using cluster reference) ### 
samtools view -h $1 > ${1}.sam
makeTagDirectory ${s}tagdir/ ${1}.sam -format sam
#findPeaks ${1}tagdir -style factor -o auto -i 181009_zhengChIP/1213K-Pooled-Input_S20_R1_001.sam 
#pos2bed.pl -o ${1}tagdir/peaks.bed ${1}tagdir/peaks.txt
#awk '{print "chr"$0}' ${1}tagdir/peaks.bed > ${1}tagdir/peaks.bed.chr.bed
#annotatePeaks.pl ${1}tagdir/peaks.bed.chr.bed mm10 > ${1}tagdir/peaksannotated.txt


