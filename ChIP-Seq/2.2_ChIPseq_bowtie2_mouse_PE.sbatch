#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=128G
#SBATCH --time=8:00:00   
#SBATCH --output=%j.stdout 
#SBATCH --error=%j.stderr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

module load bowtie2
module load samtools
module load igvtools

###hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

#file=$(ls | grep fastq.gz$)

s=$(basename $1)
t=$(dirname $1)
x=`echo $s | cut -d "_" -f 1-2`
echo $s
echo $x
###echo $file

bowtie2 -p 32 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -1 ${t}/${x}_R1_001.fastq.gz -2 ${t}/${x}_R2_001.fastq.gz -S ${t}/${x}.sam


#s=$(basename $1)
#x=`echo $s | cut -d "_" -f 1-2`

#echo $s
#echo $x
###echo $file
###echo "hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2"
#hisat2 -p 8 -x research/HG19_ERCC_hisat2/Homo_sapiens_GRCh37.73.dna.toplevel_ERCC92 -U $1 -S ${1}_hisat2.sam
#bowtie2 -p 8 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -U $1 -S ${1}.sam
##module load hisat2/2.0.4 
### hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

###s=$(basename $1)
###x=`echo $s | cut -d "_"`

#echo $s
#echo $x
#echo "hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2"
###hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2
###hisat2 -p 8 -x /data/reference/indexes/human/hg19/hisat2/grch37/Homo_sapiens.GRCh37.73.dna.toplevel -U $file -S ${file}_hisat2
samtools view -@ 8 -Sbo ${t}/${x}.sam.bam  ${t}/${x}.sam
samtools sort -@ 8 -o .${t}/${x}sam.sorted.bam ${t}/${x}.sam.bam
samtools rmdup -s ${t}/${x}.sam.sorted.bam ${t}/${x}.sam.sorted.bam.rmdup.bam
samtools index ${t}/${x}.sam.sorted.bam.rmdup.bam ${t}/${x}.sam.sorted.bam.rmdup.bam.bai
igvtools count -z 5 -w 25 -e 225 ${t}/${s}.sam.sorted.bam.rmdup.bam ${t}/${s}.sam.sorted.bam.rmdup.bam.tdf mm10


#bowtie2 -p 8 -x /data/reference/indexes/mouse/mm10/bowtie2/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}.sam
