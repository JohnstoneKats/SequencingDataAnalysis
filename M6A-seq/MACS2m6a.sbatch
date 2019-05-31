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

##module load hisat2/2.0.4 
module load macs
### hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

s=$(basename $1)
###x=`echo $s | cut -d "_" -f 1-2`

echo $s
#echo $x
echo $1


macs2 callpeak -g mm -f BAMPE -t $1 -c $2 --nomodel --gsize 2.82e8 --extsize 200 --outdir ${s}_macs2 -n ${s}_normpeaks
#echo "hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2"
###hisat2 -p 8 -x /data/reference/indexes/mouse/mm10/hisat2/grcm38/Mus_musculus.GRCm38.dna.toplevel -U $file -S ${file}_hisat2
###hisat2 -p 8 -x /data/reference/indexes/human/hg19/hisat2/grch37/Homo_sapiens.GRCh37.73.dna.toplevel -U $file -S ${file}_hisat2
#samtools view -@ 8 -Sbo ${1}.bam $1
#samtools sort -@ 8 -o ${1}.sorted.bam ${1}.bam
#samtools rmdup -s $file ${file}.rmdup.bam
#samtools index ${1}.sorted.bam ${1}.sorted.bam.bai
#igvtools count --includeDuplicates --strands read -z 5 -w 25 ${1}.sorted.bam ${1}sorted.bam.tdf hg19
