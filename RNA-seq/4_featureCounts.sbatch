#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
#SBATCH --partition=prod_short


module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o $(1)counts.txt $(1)

