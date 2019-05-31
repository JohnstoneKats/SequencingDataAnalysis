#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod


module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a /./data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o 180514_stefanRNA/180514countsUNSTRANDED.txt /./researchers/stefan.bjelosevic/RNA-QuantSeq/180514_NB501056_0130_AH73NNBGX7/*/Bam/*.sorted.bam
