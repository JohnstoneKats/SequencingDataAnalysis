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

computeMatrix reference-point  -S ${1}.bw -R peaks.bed -out ${1}peakcentered.gz -b 1000 -a 1000

#or scale to mm10 reference genes
#computeMatrix scale-regions -S ${1}.bw  -R mm10reference.bed -out ${1}scaletogenes.gz -b 1000 -a 1000

plotHeatmap -m ${1}peakcentered.gz -out ${1}peakcenteredHeatmap.pdf  

#change heatmap colours + other options
#plotHeatmap -m ${1}peakcentered.gz -out ${1}peakcenteredHeatmap.pdf --colorList 'white, purple' --whatToShow 'heatmap and colorbar' --kmeans 5 --zMin 0 --zMax 100 --samplesLabel SAMPLE1 --outFileSortedRegions  SampleKmeans5.txt




