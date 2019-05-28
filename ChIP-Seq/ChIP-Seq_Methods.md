# Methods
General methods for publicaiton purposes

## Low input ChIP


KL and KSL cells were sorted from BcorΔE9-10 (n=3), BcorWT (n=2), BcorΔE9-10KrasG12D (n=3) and KrasG12D (n=3) mice. 1-5 x 105 cells were crosslinked with 1% formaldehyde and 5% FCS in PBS at room temperature for 10 minutes, then the reactions were quenched with 1.25M glycine for 5 minutes. Crosslinked cells were washed with ice cold dilution buffer (20mM Tris pH8, 150mM NaCl, 2mM EDTA, 1%Triton-X) then resuspended in nuclear extraction buffer (50mM TRIS pH7.5, 150mM NaCl, 5mM EDTA, 0.5%NP-40, 1% TritonX-100). Cells were washed in dilution buffer again before being lysed in ChIP lysis buffer (20mM TRIS pH7.5, 150mM NaCl, 0.5M EDTA, 1%NP-40, 0.3% SDS). Chromatin was fragmented to ~150-300bp by sonicating lysates in a Covaris Ultrasonicator using 100l microtubules (Covaris). Immunoprecipitation was performed overnight at 4oC in 50% dilution and 50% lysis buffer with magnetic protein A and G Dynabeads (Life technologies) pre-coupled to 1g antibody per sample (antibodies listed in Supplementary Table 2). Beads were then washed with wash buffer 1 (20mM TRIS pH8, 500mM NaCl, 2mM EDTA, 1% Triton-X, 0.1% SDS) followed by wash buffer 2 (220mM TRIS pH 8, 250mM LiCl, 2mM EDTA, 0.5% NP-40, 0.5% deoxycholate) and then TE buffer (10mM TRIS pH7.5, 1mM EDTA) on ice. Beads were then resuspended in reverse crosslinking buffer (1% SDS, 100mM NaHCO3, 200mM NaCL, 300ug/mL proteinase K) and incubated at 55°C for 30 minutes followed by 65°C for 30 minutes. DNA was then purified using Zymo ChIP clean and concentrate kit according to manufacturer’s instructions (Zymo Research) and quantified using Qubit dsDNA HS Assay Kit (Thermo Fisher Scientific). Indexed libraries were prepared using KAPA Hyper Prep Kit for Illumina platforms (Kapa Biosystems) and the SeqCap Adapter Kit (Roche) following vendor’s instructions. Library QC and quantification was performed using D1000 high sensitivity screen tape with 4200 TapeStation Instrument (Agilent Technologies) and size selected for between 200bp and 500bp using a Pippin Prep system (Sage Science). Libraries were pooled and sequenced with 75bp single end sequencing to a depth of 15-20 x 106 reads per sample on a NextSeq (Illumina). bcl2fastq (v2.17.1.14) was used for de-multiplexing. The Fastq files generated by sequencing were aligned to the mouse reference genome (GRCm38/mm10) using bowtie2 (v2.3.3). Samtools (v1.4.1) was used for manipulation of SAM and BAM files. MACS (v2.1.1) was used for traditional peak-calling. R package Csaw (v1.10.0) was used to quantify regions enriched for histone marks, which were then associated with the closest TSS using Bedtools (v2.26). Superenhnacer analysis was performed using Homer (v4.8). Browser viewable TDF files were generated using IGVTools (v2.3.95) and ChIP-Seq tracks were visualized using IGV (v2.3.55). Graphics were generated using deeptools(v2.5.3) and R packages pheatmap(v1.0.8), plotly(v4.7.1) and ggplot2(v2.2.1).

From [Bcor loss perturbs myeloid differentiation and promotes leukaemogenesis](https://www.nature.com/articles/s41467-019-09250-6) 