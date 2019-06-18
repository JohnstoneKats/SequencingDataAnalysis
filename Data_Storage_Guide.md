
# Organising and storing your next-gen seq data

*What data do you need to keep and where?*

During data processing you will end up with a number of intemediate files including trimmed fastq, unsorted sam, bam with dulicates etc. 
It is important to remove unnessecary intermediates to reduce your digital footprint. 

It is a good idea to keep the most important files for your analysis, and a record of how each figure can be reproduced from these files. 
Ensure you keep a record of where each of the below files are located for each experiment (including runID)

### ChIP-Seq
1) Original Fastq files (These will be archived in bioinf/archive so no need to store them in your home/researchers folder) 
2) sorted bam with duplicates removed (.sorted.rmdup.bam) 
3) Index for above bam (.sorted.rmdup.bam.bai) 
4) Peak call files (.bed)
5) IGV viewable count files (.tdf)

Additional analysis files can include:
- tagdirectories (from Homer)
- bigwigs, and matrices (from deeptools)
- any plots generated directly (.pdf, . png etc.)


### RNA-Seq 
1) Original Fastq files (These will be archived in bioinf/archive so no need to store them in your home/researchers folder) 
2) sorted bam with duplicates removed (.sorted.rmdup.bam) 
3) Index for above bam (.sorted.rmdup.bam.bai) 
4) Peak call files (.bed)
5) IGV viewable count files (.tdf)

## GEO dataset SUbmission
ALL data from a project should be uploaded to GEO upon acceptance to a journal, or even before. It may be beneficial to fill in the metadata file for each sequencing run. 
More comprehensive instructions can be found [here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html) 

In short, for GEO submission you require: 

1) Raw fastq file
- Use the original fastq in the run archives

2) processed data file
- can include : .tdf, peaks.bed, .bigWig etc... 

3) metadata file
- found [here](https://www.ncbi.nlm.nih.gov/geo/info/examples/seq_template_v2.1.xls) 

* What if you have multiple sequencing "types" for a single project? * 
Fill in a metadata file for each individual -seq, e.g. 1 for ChIP-seq 1 for ATAC-seq etc. 
and upload all together in the same folder on the GEO server. Email the GEO staff to let them know you would like all of these samples linked in the same project (or "superseries")
