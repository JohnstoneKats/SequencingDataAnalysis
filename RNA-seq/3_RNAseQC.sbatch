#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=15:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org

#SBATCH --partition=prod

module load rna-seqc
module load java

RNA-SeQC -s "Test|m6a/Fastq/190222-X1A_L007/KSL-1_combined.rmdup.bam|test"  -t ref/mm10genes.gtf -o m6a/Fastq/RNAseqc -r ./../../../data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -BWArRNA ref/mm10rRNA.fa ${1}




