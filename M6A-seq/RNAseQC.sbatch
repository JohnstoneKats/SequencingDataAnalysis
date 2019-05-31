#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=2:00:00   
# SBATCH --output=%j.stdout
# SBATCH --error=%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madison.kelly@petermac.org
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

RNA-SeQC -s ${2}  -t /data/reference/indexes/mouse/ensembl_GRCm38.73/bowtie2/Tophat2_transcriptome_data/Mus_musculus.GRCm38.73.gtf  -o ${d} -singleEnd -r /data/reference/indexes/mouse/mm10/fasta/Mus_musculus.GRCm38.dna.toplevel.fa -rRNA /data/reference/indexes/mouse/ensembl_GRCm38.73/bowtie2/Tophat2_transcriptome_data/RNASeq-QC/annotation/mm10_rRNA.list




#makeTagDirectory ${1}tagdir/ $1 
#findPeaks ${1}tagdir/ -style super -o auto -superSlope -1000 -L 2 -typical ${1}normalpeaks.txt



#findPeaks ${1}tagdir -style factor -o auto -i 171127ChIPAc/Input_S6.samtagdir 
#pos2bed.pl -o ${1}.bed ${1}tagdir/peaks.txt
#awk '{print "chr"$0}' ${1}tagdir/peaks.bed > ${1}tagdir/peaks.bed.chr.bed
#annotatePeaks.pl ${1}tagdir/peaks.bed.chr.bed mm10 > ${1}tagdir/peaksannotated.txt


