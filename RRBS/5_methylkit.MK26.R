sample.id<-list("A1","A2","A3","A4","B7","B8","B9","B10")
library(methylKit)
library(genomation)
library(genomation)

##load in bam files and write to methylkit txt files for easy access later#####

sample.list<-list.files(path = "RRBS/", pattern = "*.sort.bam$", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

PROCESS_BAM_FILES = TRUE


if( PROCESS_BAM_FILES ){
  for(sample_name in sample.list){ print(sample_name);
  my.methRaw = processBismarkAln( location = paste("RRBS/",sample_name,sep = ""),
                             sample.id=paste(sample_name), assembly="mm10", 
                            read.context="CpG", save.folder="/Users/kellymadison/Documents/R/MK26RNAseq/RRBS" )}}

##read back in methylkit files 
methylkit.txt.list<-list.files( path = "RRBS/",pattern = "_CpG.txt$", all.files = FALSE,
           full.names = F, recursive = F,
           ignore.case = FALSE, include.dirs = F, no.. = T)

files = as.list( paste("RRBS/",methylkit.txt.list,sep = ""))

myobjall = methRead(files, sample.id=list("DOX1","DOX2","DOX3","DOX4","UT4","UT1","UT2","UT3"),
                    assembly="mm10", treatment=c(1,1,1,1,0,0,0,0)
)
#add chr for downstream analysis 
for (i in 1:8) {myobjall[[i]]$chr<-paste("chr",myobjall[[i]]$chr,sep = "") }

getMethylationStats(myobjall[[8]],plot=T,both.strands=FALSE)
getCoverageStats(myobjall[[7]],plot=TRUE,both.strands=FALSE)

filtered.myobj=filterByCoverage(myobjall,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
filtered.myobj
methall = unite(filtered.myobj, destrand=TRUE, min.per.group=2L)

getCorrelation(na.omit(methall),method = "pearson",plot = F)
corr<-getCorrelation(na.omit(methall),method = "pearson",plot = F)
corr<-read.table("MK26RRBS corr.txt",header = T, row.names = T)
pheatmap(corr)

clusterSamples(methall, dist="correlation", method="ward", plot=TRUE)
PCASamples(methall, screeplot=TRUE)
PCASamples(methall)
PCASamples(na.omit(methall))

library(rgl)
#####all cs annotation #####
gene.obj=readTranscriptFeatures("RRBS/mm10genes.bed")
allCann=annotateWithGeneParts(as(methall,"GRanges"),gene.obj)
allCTSS<-getAssociationWithTSS(allCann)
genomation::getTargetAnnotationStats(allCann,percentage=TRUE,precedence=TRUE)

genomation::plotTargetAnnotation(allCann,precedence=TRUE,
                     main="all c annotation")

genomation::getFeatsWithTargetsStats(allCann,percentage=TRUE)


cpg.obj=readFeatureFlank("RRBS/CpG_islands_mm10.bed", feature.flank.name=c("CpGi","shores")  )

AllCcpg=annotateWithFeatureFlank(as(methall,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
genomation::getTargetAnnotationStats(AllCcpg,percentage=TRUE,precedence=TRUE)
methall
genomation::plotTargetAnnotation(allCann,
                     main="")

genomation::getFeatsWithTargetsStats(AllCcpg,percentage=TRUE)





####calculate Diff C ######
#myDiff=calculateDiffMeth(methall,overdispersion="MN",test="Chisq",mc.cores=4)
myDiffuncorrected=calculateDiffMeth(methall,mc.cores=4)


#myDiff20=getMethylDiff(myDiff,difference=25,qvalue=0.05)

myDiff20unc=getMethylDiff(myDiffuncorrected,difference=25,qvalue=0.05)


diffMethPerChr(myDiffuncorrected,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=25)

bedgraph(myDiff20unc, file.name = "MK26dmc.bed", col.name = "meth.diff", unmeth = FALSE,
         log.transform = FALSE, negative = FALSE, add.on = "")

#####annotate diff Cs ####
annotateWithGeneParts(as(myDiff20unc,"GRanges"),gene.obj)

diffAnn=annotateWithGeneParts(as(myDiff20unc,"GRanges"),gene.obj)
TSSass<-getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

genomation::plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")

genomation::getFeatsWithTargetsStats(diffAnn,percentage=TRUE)


hypo = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                     type = "hypo",mc.cores = 4)
hyper = getMethylDiff(myDiffuncorrected, difference = 25, qvalue = 0.05, 
                      type = "hyper",mc.cores = 4)

hypoAnn<-annotateWithGeneParts(as(hypo,"GRanges"),gene.obj)
getTargetAnnotationStats(hypoAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(hypoAnn,precedence=TRUE,
                     main="hypo c annotation")

hyperAnn<-annotateWithGeneParts(as(hyper,"GRanges"),gene.obj)
getTargetAnnotationStats(hyperAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(hyperAnn,precedence=TRUE,
                     main="hyper c annotation")

######cpgi######


diffCpGann=annotateWithFeatureFlank(as(myDiff20unc,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

genomation::plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="")

getFeatsWithTargetsStats(diffCpGann,percentage=TRUE)

######-------####
CpGicounts=regionCounts(filtered.myobj,cpg.obj$CpGi)

getCoverageStats(CpGicounts[[2]],plot=TRUE,both.strands=FALSE)

methCpGi = unite(CpGicounts, destrand=TRUE, min.per.group=2L)

clusterSamples(methCpGi, dist="correlation", method="ward", plot=TRUE)
PCASamples(methCpGi)
mperc.meth.CpGi=percMethylation(methCpGi)
boxplot(mperc.meth.CpGi)
myDiffuncorrectedCpGi=calculateDiffMeth(methCpGi,mc.cores=4)
diffMethPerChr(myDiffuncorrectedCpGi,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

cpgsig<- as.data.frame(mperc.meth.CpGi[myDiffuncorrectedCpGi$qvalue<0.05,])
boxplot(cpgsig)
#pheatmap(promsig,border_color = NA,show_rownames = F,scale = "row")

myDiff25CpG=getMethylDiff(myDiffuncorrectedCpGi,difference=25,qvalue=0.05)

diffCpGAnn=annotateWithGeneParts(as(myDiff25CpG,"GRanges"),gene.obj)
TSSassdiffCpG<-getAssociationWithTSS(diffCpGAnn)
a
Diff25CpGAnn<-cbind(myDiff25CpG,TSSassdiffCpG)

library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                    ensemblRedirect = FALSE))
genes <- Diff25CpGAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25CpG<-merge(Diff25CpGAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
MK26RNA<-sigsorted
diffCpGgenes<-merge (GenenamesDiff25CpG, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
diffCpGgenesshort<-diffCpGgenes[diffCpGgenes$P.Value<0.05,]
plot(diffCpGgenesshort$meth.diff,diffCpGgenesshort$logFC)
abline(lm(diffCpGgenes$logFC~diffCpGgenes$meth.diff), col="blue") # regression line (y~x)
cor.test(diffCpGgenes$meth.diff,diffCpGgenes$logFC)



########summarise counts over promoters######
#for (i in 1:8) {myobjall[[i]]$chr<-paste("chr",myobjall[[i]]$chr,sep = "") }
promoters=regionCounts(filtered.myobj,gene.obj$promoters)

getCoverageStats(promoters[[2]],plot=TRUE,both.strands=FALSE)

methproms = unite(promoters, destrand=TRUE, min.per.group=2L)

getCorrelation(methproms,plot = F)
clusterSamples(methproms, dist="correlation", method="ward", plot=TRUE)
PCASamples(methproms, screeplot=TRUE)
PCASamples(methproms)
mperc.meth.proms=percMethylation(methproms)
boxplot(mperc.meth.proms)

myDiffuncorrectedproms=calculateDiffMeth(methproms,mc.cores=4)

promsig<- as.data.frame(mperc.meth.proms[myDiffuncorrectedproms$qvalue<0.05,])
boxplot(promsig)

boxplotshort<-as.data.frame(t(promsig))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)), aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
compare_means(value~group, boxplotmelt)

means <- aggregate(value~group, boxplotmelt, mean) 

pheatmap(promsig,border_color = NA,show_rownames = F,scale = "row")

myDiff25proms=getMethylDiff(myDiffuncorrectedproms,difference=25,qvalue=0.01)
diffMethPerChr(myDiffuncorrectedproms,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

#diffAnn=annotateWithGeneParts(as(myDiff25proms,"GRanges"),gene.obj)
#TSSass<-as.data.frame(getAssociationWithTSS(diffAnn))

#Diff25promAnn<-cbind(as(diffAnn,"GRanges"),TSSass)

diffpromAnn=annotateWithGeneParts(as(myDiff25proms,"GRanges"),gene.obj)
TSSassdiffprom<-getAssociationWithTSS(diffpromAnn)

Diff25promAnn<-cbind(myDiff25proms,TSSassdiffprom)


mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25promAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25prom<-merge(Diff25promAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
diffpromgenes<-merge (GenenamesDiff25prom, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
diffpromgenesshort<-diffpromgenes[diffpromgenes$P.Value<0.05,]
plot(diffpromgenes$meth.diff,diffpromgenes$logFC,pch= 20)
abline(lm(diffpromgenes$logFC~diffpromgenes$meth.diff), col="red") # regression line (y~x)
cor.test(diffpromgenesshort$meth.diff,diffpromgenesshort$logFC)

bedgraph(myDiff25proms, file.name = "MK26DMprom.bed", col.name = "meth.diff", unmeth = FALSE,
         log.transform = FALSE, negative = FALSE, add.on = "")

####summarise over exons#####


unique.exon<-unique(gene.obj$exons)

exons=regionCounts(filtered.myobj,unique.exon)
getCoverageStats(exons[[7]],plot=TRUE,both.strands=FALSE)

methexons = unite(exons, destrand=TRUE, min.per.group=2L)

clusterSamples(methexons, dist="correlation", method="ward", plot=F)
PCASamples(methexons)

mperc.meth.exons=percMethylation(methexons)
boxplot(mperc.meth.exons)
exsig<- as.data.frame(mperc.meth.exons[myDiffuncorrectedexons$qvalue<0.05,])
boxplot(exsig)

boxplotshort<-as.data.frame(t(exsig))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)), aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated exons")
compare_means(value~group, boxplotmelt)

pheatmap(na.omit(exsig),border_color = NA,show_rownames = F)

myDiffuncorrectedexons=calculateDiffMeth(methexons,mc.cores=4)
myDiff25exons=getMethylDiff(myDiffuncorrectedexons,difference=25,qvalue=0.01)
diffMethPerChr(myDiff25exons,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

diffexonAnn=annotateWithGeneParts(as(myDiff25exons,"GRanges"),gene.obj)
TSSassdiffexon<-getAssociationWithTSS(diffexonAnn)

Diff25exonsAnn<-cbind(myDiff25exons,TSSassdiffexon)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25exonsAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25exons<-merge(Diff25exonsAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
diffexongenes<-merge (GenenamesDiff25exons, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
diffexongenesshort<-diffexongenes[diffexongenes$P.Value<0.05,]
plot(diffexongenesshort$meth.diff,diffexongenesshort$logFC,pch= 20)
abline(lm(diffexongenesshort$logFC~diffexongenesshort$meth.diff), col="red") # regression line (y~x) 

cor.test(diffexongenesshort$meth.diff,diffexongenesshort$logFC)

bedgraph(myDiff25exons, file.name = "MK26DMexon.bed", col.name = "meth.diff", unmeth = FALSE,
         log.transform = FALSE, negative = FALSE, add.on = "")

#####summarise accross introns ######

unique.intron<-unique(gene.obj$introns)

introns=regionCounts(filtered.myobj,unique.intron)
getCoverageStats(introns[[7]],plot=TRUE,both.strands=FALSE)

methintrons = unite(introns, destrand=TRUE, min.per.group=2L)

clusterSamples(methintrons, dist="correlation", method="ward", plot=F)
PCASamples(methintrons)

mperc.meth.introns=percMethylation(methintrons)
boxplot(mperc.meth.introns)
intsig<- as.data.frame(mperc.meth.introns[myDiffuncorrectedintons$qvalue<0.05,])
boxplot(intsig)

boxplotshort<-as.data.frame(t(intsig))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)), aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated introns")
compare_means(value~group, boxplotmelt)
means <- aggregate(value~group, boxplotmelt, mean) 


#pheatmap(na.omit(exsig),border_color = NA,show_rownames = F)

myDiffuncorrectedintons=calculateDiffMeth(methintrons,mc.cores=4)
myDiff25introns=getMethylDiff(myDiffuncorrectedintons,difference=25,qvalue=0.01)
diffMethPerChr(myDiff25introns,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)

diffintronAnn=annotateWithGeneParts(as(myDiff25introns,"GRanges"),gene.obj)
TSSassdiffintron<-getAssociationWithTSS(diffintronAnn)

Diff25intronsAnn<-cbind(myDiff25introns,TSSassdiffintron)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25intronsAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25introns<-merge(Diff25intronsAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
diffintrongenes<-merge (GenenamesDiff25introns, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
diffintrongenesshort<-diffintrongenes[diffintrongenes$P.Value<0.05,]
plot(diffintrongenesshort$meth.diff,diffintrongenesshort$logFC,pch= 20)
abline(lm(diffintrongenesshort$logFC~diffintrongenesshort$meth.diff), col="red") # regression line (y~x) 

cor.test(diffexongenesshort$meth.diff,diffexongenesshort$logFC)


#summarise accross LT-HSC associated enhancers (from lara-asatio)#
LTenhancers<-read.table(file="LT-HSCenhancers.txt")
colnames(LTenhancers)<-c("chr","start","end")
LTenhancers<-as(LTenhancers,"GRanges")
LTenhancersMeth<-regionCounts(filtered.myobj,LTenhancers)

LTenhancersMethall<-unite(LTenhancersMeth, destrand=TRUE, min.per.group=2L)

clusterSamples(LTenhancersMethall, dist="correlation", method="ward", plot=T)
PCASamples(LTenhancersMethall)

mperc.meth.LTen=percMethylation(LTenhancersMethall)
boxplot(mperc.meth.LTen)

boxplotshort<-as.data.frame(t(mperc.meth.LTen))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)) , aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated LTHSC enhancers")
compare_means(value~group, data=subset(boxplotmelt, !is.na(value)))


pheatmap(na.omit(mperc.meth.LTen),scale = "row")

DiffLTen=calculateDiffMeth(LTenhancersMethall,mc.cores=4)
diffMethPerChr(DiffLTen,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)
LTsig<- as.data.frame(mperc.meth.LTen[DiffLTen$pvalue<0.05&DiffLTen$meth.diff>10,])
pheatmap(na.omit(LTsig),border_color = NA,show_rownames = F, scale = "row")
boxplot(LTsig)

myDiff25LT=getMethylDiff(DiffLTen,difference=25,qvalue=0.01)

##annotate

difLTAnn=annotateWithGeneParts(as(myDiff25LT,"GRanges"),gene.obj)

TSSassdiffLT<-getAssociationWithTSS(difLTAnn)

Diff25LTAnn<-cbind(myDiff25LT,TSSassdiffLT)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25LTAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25lt<-merge(Diff25LTAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
diffltgenes<-merge (GenenamesDiff25lt, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
diffLTgenesshort<-diffltgenes[diffltgenes$P.Value<0.05,]
plot(diffLTgenesshort$meth.diff,diffLTgenesshort$logFC,pch= 20)
abline(lm(diffLTgenesshort$logFC~diffLTgenesshort$meth.diff), col="red") # regression line (y~x) 

cor.test(diffLTgenesshort$meth.diff,diffLTgenesshort$logFC)


#####accross MP enhancers #######
CMP.SE<-read.table(file="CMP-SE.txt",header =T)
colnames(CMP.SE)<-c("chr","start","end")
CMPse<-as(CMP.SE,"GRanges")
CMPeMeth<-regionCounts(filtered.myobj,CMPse)

CMPseMethall<-unite(CMPeMeth, destrand=TRUE, min.per.group=2L)

clusterSamples(CMPseMethall, dist="correlation", method="ward", plot=T)
PCASamples(CMPseMethall)

mperc.meth.CMPs=percMethylation(CMPseMethall)
boxplot(mperc.meth.CMPs)
#pheatmap(na.omit(mperc.meth.CMPs),scale = "row")

DiffCMP=calculateDiffMeth(CMPseMethall,mc.cores=4)
diffMethPerChr(DiffCMP,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=25)
CMPse<- as.data.frame(mperc.meth.LTen[DiffLTen$pvalue<0.05&DiffLTen$meth.diff>10,])
pheatmap(na.omit(LTsig),border_color = NA,show_rownames = F, scale = "row")
boxplot(LTsig)
myDiff25CMP=getMethylDiff(DiffCMP,difference=25,qvalue=0.01)

difCMPAnn=annotateWithGeneParts(as(myDiff25CMP,"GRanges"),gene.obj)

TSSassdifcmp<-getAssociationWithTSS(difCMPAnn)

Diff25CMPAnn<-cbind(myDiff25CMP,TSSassdifcmp)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl",host = "www.ensembl.org",
                                                     ensemblRedirect = FALSE))
genes <- Diff25CMPAnn$feature.name
G_list <- getBM(filters= "refseq_mrna", attributes= c("refseq_mrna","external_gene_name"),values=genes,mart= mart)

GenenamesDiff25cmp<-merge(Diff25CMPAnn,G_list,by.x="feature.name",by.y="refseq_mrna")

#integrate with RNA? 
diffcmpgenes<-merge (GenenamesDiff25cmp, MK26RNA, by.x= "external_gene_name", by.y ="row.names")
#diffexongenesshort<-diffexongenes[diffexongenes$P.Value<0.05,]
plot(diffcmpgenes$meth.diff,diffcmpgenes$logFC,pch= 20)
abline(lm(diffcmpgenes$logFC~diffcmpgenes$meth.diff), col="red") # regression line (y~x) 

cor.test(diffcmpgenes$meth.diff,diffcmpgenes$logFC)


#####average meth percentage accross diff gene feaatures #### 
install.packages("rowr")
library(rowr)
mperc.meth.av<-cbind.fill(rowMeans(mperc.meth.proms[,5:8]),rowMeans(mperc.meth.exons[,5:8]),
           rowMeans(mperc.meth.CpGi[,5:8]),rowMeans(mperc.meth.introns[,5:8]),
           rowMeans(mperc.meth.LTen[,5:8])
)
colnames(mperc.meth.av)<-c("proms","exons","CpGi","introns","HSCenhancers")
boxplot(log(mperc.meth.av),pch = 20,cex = 0.4)

######perc meth#######
mperc.meth=percMethylation(methall)
allmeth<-getData(methall)
locperc.meth<-cbind(allmeth[,1:3],mperc.meth)


diffAnn=annotateWithGeneParts(as(locperc.meth,"GRanges"),gene.obj)
TSSass<-getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

anlocperc.meth<-merge(TSSass,locperc.meth,by.x = "target.row",by.y = "row.names")

ordered<-as.data.frame(mperc.meth[order(myDiffuncorrected$qvalue),])
allsig<- as.data.frame(mperc.meth[myDiffuncorrected$qvalue<0.001&abs(myDiffuncorrected$meth.diff)>25,])

boxplot(allsig[,c(5:8,1:4)])

kmeans(allsig,4)
pheatmap(ordered[1:100,],border_color = NA,show_rownames = F)
pheatmap(allsig, border_color = NA,show_rownames = F)

install.packages("dendextend")
library(dendextend)

my_hclust_gene <- hclust(dist(allsig), method = "complete")
list <- as.data.frame(cutree(tree = as.dendrogram(my_hclust_gene), k = 6))

my_gene_col <- data.frame(cluster = ifelse(test = list == 1, yes = "cluster 1", no = "cluster 2"))
my_gene_col$cluster<-paste("cluster",list$`cutree(tree = as.dendrogram(my_hclust_gene), k = 6)`,sep = " ")
my_gene_col$cluster<-factor(my_gene_col$cluster)
my_gene_col<-as.data.frame(my_gene_col[,-1])

pheatmap(allsig, border_color = NA,show_rownames = F,annotation_row = my_gene_col)
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)


###ggplot boxplot###
boxplotshort<-as.data.frame(t(allsig))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                                 levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data = boxplotmelt, aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
compare_means(value~group, data=subset(boxplotmelt, !is.na(value)))
t.test 
#######boxplot all methylat##
library(ggpubr)
library(reshape2)
library(ggplot2)
library(ggpubr)
boxplotshort<-as.data.frame(t(as.data.frame(mperc.meth)))
boxplotshort$group<-c("DOX","DOX","DOX","DOX","UT","UT","UT","UT")
boxplotshort$group <- factor(boxplotshort$group,
                             levels = c("UT","DOX"),ordered = TRUE)
boxplotmelt<- melt(boxplotshort, id.var = "group")

ggplot(data=subset(boxplotmelt, !is.na(value)) , aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) +
  theme_light() + scale_fill_brewer(palette="Set1")+ xlab("") +
  ylab("percentage methylation")+ggtitle("Differentially methylated cytosines")
boxplot(mperc.meth)
compare_means(value~group, boxplotmelt,method = "wilcox.test")
  
######circos ######

library(BSgenome)
library("BSgenome.Mmusculus.UCSC.mm10")

######altered ###### works faster once myDIff is already made#####
chrom.length = seqlengths(Mmusculus)  
chr.len = chrom.length[grep("_|chrM|chrX|chrY", names(chrom.length), invert = T)] 

ideoDMC(myDiffuncorrected, chrom.length = chr.len, difference = 20, qvalue = 0.05,circos = T, title = "test", hyper.col = "red", hypo.col = "blue")


ideoDMC <- function(myDiff, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  hypo = getMethylDiff(myDiff, difference = difference, qvalue = qvalue, 
                           type = "hypo",mc.cores = 4)
  hyper = getMethylDiff(myDiff, difference = difference, qvalue = qvalue, 
                            type = "hyper",mc.cores = 4)
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per,pruning.mode = "coarse") = seqlevels(myIdeo)
  seqlengths(g.per)=(chr.len)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po,pruning.mode = "coarse") = seqlevels(myIdeo)
  seqlengths(g.po)=(chr.len)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", size = 0.1, aes(x = midpoint, y = meth.diff, color = id), radius = 25, trackWidth = 30) + scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames),vjust = 0, radius = 55, trackWidth = 7) + labs(title = "Dox vs UT")
    
  } else {
    
    p <- ggplot()+  layout_karyogram(myIdeo, cytoband = FALSE)
    p  + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                          aes(x = midpoint, 
                              y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                           hypo.col)) + labs(title = title)
    # new alternative commented out
    #autoplot(c(g.po, g.per), layout = "karyogram", geom = "point", size = 0.65, 
    #aes(x = midpoint,y = meth.diff, color = id))   scale_colour_manual(values = c(hyper.col, 
    #                                                                                        hypo.col))  labs(title = title)
    
  }
}


#####perc meth plot##### 

