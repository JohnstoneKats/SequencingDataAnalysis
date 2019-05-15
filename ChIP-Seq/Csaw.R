library(csaw)
library(pheatmap)
library(ggpubr)
library(reshape2)
library(GenomicRanges)
library(genomation)

bamfiles <- list.files( pattern = ".bam$")
data_all <- windowCounts(bamfiles[c(1:5,21)], ext=200, width=1000,bin=T)

#filter out low logcpm values
data_keep<- aveLogCPM(asDGEList(data_all)) >= 0
data_all <- data_all[data_keep,]

#filter out anything which doesn't surpass input (with FC of >3)
input<-data_all [,6]
IP<-data_all [,c(1:5)]
filter.stat<-filterWindows(IP, input, type="control",prior.count =5,norm.fac=list(IP, input))
keep <- filter.stat$filter > log2(2)
data<-IP[keep,]


#totals and design
totals <- data@colData$totals

#extract coordinates from filtered data
x<-data@rowRanges

#make sample info 
SI<-matrix(data = NA,nrow=5,ncol=2)
SI[,1]<-c(1:5)
SI[,2]<-c("B","B","B","A","A")
colnames(SI)<-c("Sample","Treatment")
SI<-as.data.frame(SI)


#create concatination of group for design matrix
group <- paste(SI$Bcor,sep=".")
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
cont.matrix <- makeContrasts(AvsB =B-A,
                             levels=design)



norm <- normOffsets(data, lib.sizes=totals)
y <- asDGEList(data, norm.factors=norm)

y <- estimateDisp(y, design,trend.method = "loess")

fit <- glmQLFit(y, design, robust=TRUE)

results <- glmQLFTest(fit,contrast = cont.matrix)
Results<-results$table

adj.counts<- cpm(y, log=TRUE)

colnames(adj.counts)<-c("B1","B2","B4","A1","A2")

plotMDS(adj.counts, top=100000,col = c("red","red","red","blue","blue"))

adj.counts_sortedbypvalue <- adj.counts[order(Results$PValue),]

merged<-merge(results,adj.counts,by= "row.name")
rownames(merged) <- merged$Row.names
merged$Row.names <- NULL

Ranges<-data.frame(data@rowRanges)
mergedlocation <- merge(merged,Ranges,by="row.names")
mergedlocation<-mergedlocation[order(mergedlocation$PValue),]
filter<-mergedlocation[(mergedlocation$PValue<0.01 & mergedlocation$logFC<(-0.3)),]
filterlocation<-filter[,11:15]

chr<-paste("chr",filterlocation$seqnames,sep="")
filterlocation$seqnames<-chr

x<-data.frame(data_@rowRanges)
chr <- paste("chr",x$seqnames,sep="")
x$seqnames<-chr
write.table(x,file = "alllocation.bed",row.names = F,col.names = F, quote = F,sep = "\t")


pheatmap(adj.counts_sortedbypvalue[1:100,],scale = "row",show_rownames = F,border_color = NA,
         main = c("sig diff regions"),file = "Pheatmap.pdf")

t.test(adj.counts_sortedbypvalue[1:250,3:4])
write.csv(adj.counts_sortedbypvalue,file = "countsbyPub.csv")

merged <- mergeWindows(rowRanges(data), tol=100L)
tabcom <- combineTests(merged$id, results$table)
tabcom1 <- cbind(tabcom,merged$region)

#########ggplot boxplot ##########
boxplot(adj.counts_sortedbypvalue[1:100,])
boxplotshort<-as.data.frame(t(adj.counts_sortedbypvalue[1:100,]))
boxplotshort$group<-c("B","B","B","A","A")
boxplotshort$group <- factor(boxplotshort$group,
                                 levels = c("A","B"),ordered = TRUE)
melt<- melt(boxplotshort, id.var = "group")

ggplot(data = melt, aes(x=group, y=value)) + geom_boxplot(aes(fill=group)) + 
  theme_light() + scale_fill_brewer(palette="Set1",direction=-1)+ xlab("") +
  ylab("Normalised CPM")+ 
  stat_compare_means(method = "t.test",aes(label = paste('p =', ..p.format..)))+ 
  ggtitle("Top 100 regions")

#########################################annotate with genomic region################################################################



filterublocation$chr<-paste ("chr",filterublocation$chr,sep = "")
y<-as(filterublocation,"GRanges")

filter27melocation<-merged27melocation[merged27melocation$PValue<0.01,11:15]
l<-as(filter27melocation,"GRanges")

f<-as(x,"GRanges")

filter4melocation<-merged4melocation[merged4melocation$PValue<0.01,11:15]
q<-as(filter4melocation,"GRanges")

refseq_mm10<-readTranscriptFeatures("/home/mkelly/research/mm10 genes (1).bed")

anot = annotateWithGeneParts(y,refseq_mm10)

plotTargetAnnotation(anot ,precedence=TRUE,
                     main="differential ub location", cex.legend = 0.6)

#annotate with h3k4me3 or h3k27me3 overlap

anotk27<-annotateWithFeature(l,y)



