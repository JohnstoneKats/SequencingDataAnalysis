#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### If experimenting: use srun --pty /bin/bash

#s=$(basename $1) 
#x=`echo $s | cut -d "_" -f 1-2`
echo ${1}
module load deeptools

computeMatrix scale-regions -S 190430_Pol2/*.bw -R ref/hg19HighTHP1exp.bed -out 190430_Pol2/highexpgenes.gz -b 1000 -a 1000

#bamCoverage -p 8 -e 225 --normalizeUsingRPKM -b ${1} -o ${1}.bw

#computeMatrix reference-point -S EGR2/SRA/SRR4345670.sam.sorted.bam.rmdup.bamatac.bw EGR2/SRA/SRR4345670.sam.sorted.bam.rmdup.bamatac.bw EGR2/egr2.bw -R SRR4345670_macs2/SRR4345670_peaks_summits.bed -out EGR2/atacpeakcentred.gz -b 1000 -a 1000

#computeMatrix reference-point -S EGR2/egr2.bw  -R mm10TSSmodifiedmerged.bed -out 190430_Pol2/allgenes.gz -b 1000 -a 1000

plotHeatmap -m 190430_pol2/highexpgenes.gz -out 190430_Pol2/highexpheatmaptest.pdf   

#s=$(basename $1)
#x=`echo $s | cut -d "_" -f 1-2`

#echo $s
#echo $x
#echo $1




