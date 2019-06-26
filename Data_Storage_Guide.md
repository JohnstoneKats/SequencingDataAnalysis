
# Organising and storing your next-gen seq data

*What data do you need to keep and where?*

During data processing a number of intermediate files are generated including: trimmed fastq, unsorted sam, bam with dulicates etc. It is important to remove unnessecary intermediates to reduce your digital footprint. Many of the scripts in this github incude the removal of unnessecary files, however, ensure you also remove any additional files you may create. Most files should be stored in your researchers folder(researchers/user.name), not your home directory (homer/uname).

For each sample, it is a good idea to keep the **raw data file** (i.e. fastq file), at least one **processed data file** (.bam or .bed) and the **code** used to generate the figures so they can be reproduced from these files. 

As a general guide, keep a record of where each of the below files are located for each experiment (including runID):

### ChIP-Seq
1) Fastq files (These will be archived in bioinf/archive so no need to store them in your home/researchers folder) 
2) sorted bam with duplicates removed (.sorted.rmdup.bam) 
3) Index for above bam (.sorted.rmdup.bam.bai) 
4) Peak call files (.bed)
5) IGV viewable count files (.tdf)
6) sbatch scripts used to generate above files and/or figures

Additional analysis files can include:
- tag directories (from Homer)
- motif results files (from Homer)
- bigwigs, and matrices (from deeptools)
- any plots generated directly (.pdf, .png etc.)

Files which can be removed: 
- .sam files including sorted sam generated during alignment

### RNA-Seq 
1) Original Fastq files (These will be archived in bioinf/archive so no need to store them in your home/researchers folder) 
2) sorted bam (.bam) 
3) Index for above bam (.bam.bai) 
4) counts table (.txt)
5) sbatch script used to generate counts file
6) R scripts used to analyse counts and generate figures

Additional analysis files can include:
- Voom normalised counts table
- Results table 

Files which can be removed: 
- .sam files including sorted sam

## GEO dataset Submission
ALL data from a project should be uploaded to GEO upon acceptance to a journal, or even before. It may be beneficial to fill in the metadata file for each sequencing run as the analysis is performed. 
More comprehensive instructions can be found [here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html) 

In short, for GEO submission you require: 

**1) Raw fastq files**
- Use the original fastq in the run archives

**2) Processed data files**
- can include : .tdf, peaks.bed, .bigWig, counts.txt etc... 

**3) Metadata file**
- found [here](https://www.ncbi.nlm.nih.gov/geo/info/examples/seq_template_v2.1.xls) 


*What if you have multiple sequencing "types" for a single project?* 
Fill in a metadata file for each individual -seq, e.g. 1 metadata sheet for ChIP-seq and 1 for ATAC-seq etc. 
and upload all together in the same folder on the GEO server. Email the [GEO staff](geo@ncbi.nlm.nih.gov) to let them know you would like all of these samples linked in the same project (or "superseries"). The GEO staff are generally very helpful. 
