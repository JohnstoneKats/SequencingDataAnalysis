#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=4:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
#SBATCH --partition=prod


 
module load rna-seqc
module load java

s=$(basename $1)
d=$(dirname $1)

echo ${1}
echo ${d}


RNA-SeQC -s ${2}  -t /data/reference/indexes/mouse/ensembl_GRCm38.73/bowtie2/Tophat2_transcriptome_data/Mus_musculus.GRCm38.73.gtf  -o ${d} -r /data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -rRNA /data/reference/indexes/mouse/ensembl_GRCm38.73/bowtie2/Tophat2_transcriptome_data/RNASeq-QC/annotation/mm10_rRNA.list
