#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.name@petermac.org
#SBATCH --partition=prod

### input is list of .bam filenames generated automatically by seqliner or manually by hisat
### ensure correct reference genome is used
### change your.name (above) to you own email to receive alerts

module load subread

featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a /./data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o counts.txt ${1}
