#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00   
#SBATCH --output=%j.stdout
#SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod_short

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

##module load hisat2/2.0.4 
module load macs
### hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

s=$(basename $1)
###x=`echo $s | cut -d "_" -f 1-2`

echo $s
#echo $x
echo $1 


macs2 callpeak -g mm -f BAM -t $1 -c 181104_MK26RNA/SRA/SRR2132490.fastq.gz.sam.sorted.bam.rmdup.bam --extsize 225 --broad -q 0.01 --broad-cutoff 0.1 --cutoff-analysis --outdir ${1}BROADq01q1macs2 -n ${s}BROADq01q1peaks

#macs2 callpeak -g mm -f BAM -t $1 -c 181104_MK26RNA/SRA/SRR2132490.fastq.gz.sam.sorted.bam.rmdup.bam --extsize 225 --nomodel --outdir ${1}NOMODmacs2 -n ${s}NOMODpeaks

