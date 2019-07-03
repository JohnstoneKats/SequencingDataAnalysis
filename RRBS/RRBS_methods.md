## RRBS Methods
*altered From Maddy's Thesis*

Genomic DNA was extracted from FACS sorted tumour cells using Gentra Puregene Blood kit
(Qiagen, cat: 158467) as per manufacturer’s instructions. Resulting DNA was quantified using
a Qubit 4 fluorometer (Thermo Fisher) with dsDNA high sensitivity kit (Thermo Fisher, cat:
Q32854) according to manufacturer’s protocol.
2.4.8.2 RRBS MSP1 Digestion
A mastermix of MSP1 digest components was prepared according to the number of samples as
follows:

CutSmart Buffer 3ul , NEB, cat: B7204
Msp I (20,000 U/ml) 1ul,  NEB, cat: R0106S
H2O 3ul

7 μl of mastermix was added to 100ng of gDNA diluted in 23 μl H2O in a PCR strip tube.
DNA was digested in a thermocycler at 37°C overnight with a heated lid set to 47°C.


For end repair and poly-a tailing, a master mix was prepared

Klenow (3’à5’ exo-) 1ul NEB cat: m0212
dNTPs (10mM dATP, 1mM dCTP, 1 mM dGTP) 1ul  Promega, cat: U1330

2 μl of above mastermix was added to each MSP1 digested sample and the mix was incubated
in a thermocycler for 20minutes at 30°C for end repair, then immediately followed by 20
minutes at 37°C for poly a-tailing.

The reaction was purified using a 2x volume of AMPure XP beads (Beckman Coulter, cat:
A63881) according to manufacturer’s instructions. In short, beads were mixed with the above
reaction in a1.5ml Eppendorf tube incubated at room temperature for 15-30 minutes. A
DynaMag-2 Magnet (ThermoFisher, cat: 12321D) was used to separate the magnetic beads
from the supernatant and the beads were washed with 150 μl 80%(v/v) ethanol twice. After
completely removing ethanol the beads were dried for up to 10 minutes. DNA was eluted from
the beads with 20 μl of elution buffer (EB, Qiagen, [10 mM Tris-Cl, pH 8.5]).


To allow for multiplexing of samples, methylated True-seq nano LT adapters (Illumina, cat:
20015964) were ligated to each sample DNA (see manufacturers website for adapter
sequences). 2μl of 1/20 diluted adapter was added to each sample in a PCR strip tube then 8μl
of ligation mastermix (below) was added and samples were mixed.

T4 DNA Ligase Buffer 3ul NEB cat: N0202, 
T4 DNA Ligase 1ul  NEB cat: N0202, 
H2O 4

The samples were then incubated at 16°C overnight in a thermocycler with heated lid.
Following incubation, samples were purified with 2x volume of AMPure XP beads as above
and eluted in 11μl of EB.
To quantify adapter ligation efficiency of each sample for appropriate pooling, qPCR was
performed with 1/4000 dilution of eluate using 1.1 and 2.1 Illumina library quantification primers (below) and a SensiFAST SYBR High-ROX kit as described in section 2.4.4 above.
Fold change of each sample relative to the sample with the lowest amount of adapter ligated
DNA (highest Ct value) was calculated. 10μl of the sample with the lowest concentration and
a relatively equal amount of DNA from each sample (as calculated above) was pooled together
in a single 1.5ml Eppendorf tube. AMPure XP beads were then used as described above to
reduce the total volume to 20μl.


qPCR primer 1.1 AATGATACGGCGACCACCGAGAT, 
qPCR primer 2.1 CAAGCAGAAGACGGCATACGA

Bisulfite conversion was performed with EZ DNA methylation-gold Kit (Zymo Research, cat:
D5005) according to manufacturer’s protocol. In short, samples were bisulfite converted over
3 hours with the provided conversion reagent, then bound to a provided column via
centrifugation. DNA bound to the column was then desulfonated and washed then eluted with
100μl EB pre-warmed to 55°C.


The resulting bisulfite converted samples were PCR amplified by adding 100μl of the
following 2x mastermix to 100μl of sample and incubating in a thermocycler for the below
conditions.

10x Pfu Turbo Cx Buffer 20ul,  Agilent Technologies, cat: 600410, 
dNTP (10mM each) 5ul  Promega, cat: U1511, 
F+R TruSeq primer mix (2.5uM) 20ul KAPA, cat: KK2623, 
Pfu Turbo Cx DNA Polymerase 4ul  Agilent Technologies, cat: 600410, 
H2O 51ul, 

Initial Denaturation 95°C 2 min,
95°C 30s,
15 cycles 65°C 30s,
72°C 45s,
Final extension 72°C 7 min,
Hold 4°C

PCR products were then purified with two rounds of AMPure XP beads, as described above.
120 μl of AMPure XP beads (1.2x) was added to 200μl of PCR products, washed twice and
eluted in 40μl. Then, 60ul AMPure XP beads (1.5x) were added and washed and the final
library was eluted in 25μl EB.
QC and library quantification was performed by diluting 1μl of the final library 1/10 and using
a D1000 high sensitivity screen tape with 2200 TapeStation Instrument (Agilent Technologies).
The average size of resulting RRBS libraries was 350bp, with low levels of primer dimers
(~150bp). Size selection for fragments between 200bp and 700bp was performed using a Pippin
Prep (Sage Science).


Libraries were sequenced 75bp PE on the Nextseq500 to a depth of 10-20 x106 reads per
sample. Sequencing reads were de-multiplexed using bcl2fastq (v2.17.1.14) and low-quality
reads Q<30 were removed. Fastqc (v 0.11.5) was used for quality control.


Fastq files were trimmed using trimgalore (v0.4.4) with the “--rrbs” option. Bismark (v0.18.1,
Babraham bioinformatics(Krueger & Andrews, 2011)) and bowtie2 (v 2.3.3) was used to align
the Fastq files to a bisulfite converted mm10 reference genome and assess conversion rate. The
resulting BAM files with methycall strings were then converted to MethylKit text files using
methylkit (v1.4.1) which were then annotated using bedtools (v2.26) or R package
GenomicFeatures(v1.28.5). Further analysis was achieved using R packages methylkit (v1.4.1)
genomation (v1.11.3) and GenomicRanges (v1.30.1) and visualised using ggplot2 (v2.2.1) and
illustrator.
