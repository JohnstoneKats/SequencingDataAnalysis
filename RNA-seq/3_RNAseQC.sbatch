#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=20:00:00   
# SBATCH --output=%j.stdout
# SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@petermac.org
# SBATCH --workdir=$1
#SBATCH --partition=prod

### If experimenting: use srun --pty /bin/bash

 
module load rna-seqc
module load java

s=$(basename $1)
#x=`echo $s | cut -d "_" -f 1-2`
d=$(dirname $1)
#echo $s
#echo $x
echo ${1}
echo ${d}

#o# uses homer to make tagdir (from sam) , call peaks, then annotate peaks ### (also adds chr to bed file if using cluster reference) ### 

RNA-SeQC -s "Test|m6a/Fastq/190222-X1A_L007/KSL-1_combined.rmdup.bam|test"  -t ref/mm10genes.gtf -o m6a/Fastq/RNAseqc -r ./../../../data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -BWArRNA ref/mm10rRNA.fa ${1}




