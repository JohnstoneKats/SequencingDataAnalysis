#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=12:00:00   
# SBATCH --output=%j.stdout
# SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### Run this command from the folder with the sbatch files
### for dir in `ls`; do sbatch trimmed_anbalysis.sbatch $dir; done
### indices files: /data/reference/indexes/mouse/mm10/hisat2

### If experimenting: use srun --pty /bin/bash

 
module load homer


s=$(basename $1)
x=`echo $s | cut -d "_" -f 1`

#echo $s
#echo $x
echo $1

findMotifsGenome.pl $1 mm10 171127ATAC/$(x)motif -preparsedDir 171127ATAC/

#macs2 callpeak -g mm -f BAM -t $1 -c 171012_ChIP-seq/InputMK_S25_R1_001.fastq.gz_bowtie2.sam.sorted.bam.rmdup.bam --nomodel --extsize 225 --broad-cutoff 1e-6 --broad --outdir ${1}_macs2 -n ${1}_broadpeaks
#echo "hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2"
###hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2
###hisat2 -p 8 -x /data/reference/indexes/human/hg19/hisat2/grch37/Homo_sapiens.GRCh37.73.dna.toplevel -U $file -S ${file}_hisat2
#samtools view -@ 8 -Sbo ${1}.bam $1
#samtools sort -@ 8 -o ${1}.sorted.bam ${1}.bam
#samtools rmdup -s $file ${file}.rmdup.bam
#samtools index ${1}.sorted.bam ${1}.sorted.bam.bai
#igvtools count --includeDuplicates --strands read -z 5 -w 25 ${1}.sorted.bam ${1}sorted.bam.tdf hg19
