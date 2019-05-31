#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=20:00:00   
# SBATCH --output=%j.stdout
# SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### If experimenting: use srun --pty /bin/bash

 
module load cutadapt

s=$(basename $1)
x=`echo $s | cut -d "_" -f 1-2`
t=$(dirname $1)

#echo $s
#echo $x
echo ${1}

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 30,30 -m 30 -o ${t}/${x}_R1.trimmed.fastq.gz -p ${t}/${x}_R2.trimmed.fastq.gz  ${t}/${x}_R1.fastq.gz ${t}/${x}_R2.fastq.gz
 


#makeTagDirectory ${1}tagdir/ $1 
#findPeaks ${1}tagdir/ -style super -o auto -superSlope -1000 -L 2 -typical ${1}normalpeaks.txt



#findPeaks ${1}tagdir -style factor -o auto -i 171127ChIPAc/Input_S6.samtagdir 
#pos2bed.pl -o ${1}.bed ${1}tagdir/peaks.txt
#awk '{print "chr"$0}' ${1}tagdir/peaks.bed > ${1}tagdir/peaks.bed.chr.bed
#annotatePeaks.pl ${1}tagdir/peaks.bed.chr.bed mm10 > ${1}tagdir/peaksannotated.txt


