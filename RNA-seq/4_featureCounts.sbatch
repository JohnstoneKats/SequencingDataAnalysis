#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=12:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash


module load subread
#featureCounts -t exon -g gene_id -F 'GTF/GFF' -p -O -M -B –C -s 0 -a /./researchers/stephin.vervoort/HG19_ERCC_hisat2/Homo_sapiens.GRCh37.73_ERCC92.gtf -o stephinINTS6RNA/Stephincounts.txt stephinINTS6RNA/*.sorted.bam
featureCounts -t exon -g gene_id -F 'GTF/GFF' -O -M -s 0 -a /./data/reference/gtf/Mus_musculus.GRCm38.90.gtf  -o 180514_stefanRNA/180514countsUNSTRANDED.txt /./researchers/stefan.bjelosevic/RNA-QuantSeq/180514_NB501056_0130_AH73NNBGX7/*/Bam/*.sorted.bam

