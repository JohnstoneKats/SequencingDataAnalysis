#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short


echo ${1}
module load deeptools

bamCoverage -p 8 -e 225 --normalizeUsingRPKM -b ${1} -o ${1}.bw

computeMatrix scale-regions -S 190430_Pol2/*.bw -R ref/hg19HighTHP1exp.bed -out 190430_Pol2/highexpgenes.gz -b 1000 -a 1000

#computeMatrix reference-point -S EGR2/egr2.bw  -R mm10TSSmodifiedmerged.bed -out 190430_Pol2/allgenes.gz -b 1000 -a 1000

plotHeatmap -m 190430_pol2/highexpgenes.gz -out 190430_Pol2/highexpheatmaptest.pdf   



