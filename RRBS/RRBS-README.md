Reduced representation bisulfite sequencing analysis - currently only for mouse (mm10)
1. fastQC 

2. Trimming
(Currently PE only) 

3. Alignment

	3.1 Paired-end alignment - default bismark settings

	3.2 Single-end alignment - default bismark settings

	3.3 Paired-end alignment, saves unaligned/discordant pairs, realigns leftover reads as SE reads - default bismark settings

4. Make bedgraph for IGV visualisation - using Bismarks methylation extractor

5. Differential methylation analysis in R with Methylkit, generates all results files, objects and basic QC plots 

6. Figure generation, requires objects generated in (5.) outputs predominantly .pdf files
