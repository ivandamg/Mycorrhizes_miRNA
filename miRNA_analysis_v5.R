#####################################################
########################
#######
####### miRNA seq analysis with DEseq2
# V4
# ivan mateus
# sources:
#http://www.bioconductor.org/help/workflows/rnaseqGene/#deseq2-import-functions
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
#http://www.bioconductor.org/help/workflows/rnaseqGene/

# count reads in bash
#for i in $(ls *fastq.gz); do zcat $i | echo $((`wc -l`/4)); done
#for i in $(ls *fastq.gz); do cat $i | echo $((`wc -l`/4)); done

#######################################################
##libraries
#source("https://bioconductor.org/biocLite.R")

#install.packages(c("DESeq2","pheatmap","RColorBrewer","ggplot2"))
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library(pheatmap)
library("RColorBrewer")
library(ggplot2)



###############################
#??import data
setwd('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/mesc_ref_only_DE/') 

setwd('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/cassava_mapped_miRNA/') 
filesToProcess <- dir(pattern = "*res.csv$")  #files to process.
#import files
listOfFiles <- lapply(filesToProcess, function(x) read.csv(x, header=T,sep="\t"))
names(listOfFiles) <- gsub('miRNAs_expressed_all_samples_','',gsub('_res.csv$','',filesToProcess,perl=TRUE),perl=TRUE)
names(listOfFiles) 


###############################
#Summary results sequencing
summa<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/summary_seq_results_v1.txt",h=T)
head(summa)
summa2<-summa[summa$Discard=="No",]
summa2<-na.omit(summa2)
head(summa2)
barplot(summa2$microRNA_10reads,names.arg = summa2$Name,las=2)
boxplot(summa2$microRNA_10reads~interaction(summa2$AMF,summa2$PLANT),las=2)
boxplot(summa2$microRNA_10reads~summa2$Lane,las=2)
barplot(summa2$microRNA_10reads,names.arg = summa2$Name,las=2)
par(mfrow=c(1,3))
boxplot(summa2$microRNA_10reads[summa2$PLANT=="BRA337"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT=="BRA337"&summa2$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
boxplot(summa2$microRNA_10reads[summa2$PLANT=="COL2215"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT=="COL2215"&summa2$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
boxplot(summa2$microRNA_10reads[summa2$PLANT!="BRA337"&summa2$PLANT!="COL2215"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT!="BRA337"&summa2$PLANT!="COL2215"&summa2$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
kruskal.test(summa2$microRNA_10reads[summa2$PLANT=="BRA337"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT=="BRA337"&summa2$AMF!="CO_INOC"])
kruskal.test(summa2$microRNA_10reads[summa2$PLANT=="COL2215"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT=="COL2215"&summa2$AMF!="CO_INOC"])
kruskal.test(summa2$microRNA_10reads[summa2$PLANT!="BRA337"&summa2$PLANT!="COL2215"&summa2$AMF!="CO_INOC"]~summa2$AMF[summa2$PLANT!="BRA337"&summa2$PLANT!="COL2215"&summa2$AMF!="CO_INOC"])

summa2$Name
summa2$Sample
###############################
#Data pre-process
listOfFiles<-listOfFiles[summa2$Sample]
df_vf<-cbind.data.frame(listOfFiles[[1]][,c(1,3,2,6)],listOfFiles[[2]][,c(2,6)],listOfFiles[[3]][,c(2,6)],listOfFiles[[4]][,c(2,6)],
                 listOfFiles[[5]][,c(2,6)],listOfFiles[[6]][,c(2,6)],listOfFiles[[7]][,c(2,6)],listOfFiles[[8]][,c(2,6)],
                 listOfFiles[[9]][,c(2,6)],listOfFiles[[10]][,c(2,6)],listOfFiles[[11]][,c(2,6)],listOfFiles[[12]][,c(2,6)],
                 listOfFiles[[13]][,c(2,6)],listOfFiles[[14]][,c(2,6)],listOfFiles[[15]][,c(2,6)],listOfFiles[[16]][,c(2,6)],
                 listOfFiles[[17]][,c(2,6)],listOfFiles[[18]][,c(2,6)],listOfFiles[[19]][,c(2,6)],listOfFiles[[20]][,c(2,6)],
                 listOfFiles[[21]][,c(2,6)],listOfFiles[[22]][,c(2,6)],listOfFiles[[23]][,c(2,6)],listOfFiles[[24]][,c(2,6)],
                 listOfFiles[[25]][,c(2,6)],listOfFiles[[26]][,c(2,6)],listOfFiles[[27]][,c(2,6)],listOfFiles[[28]][,c(2,6)],
                 listOfFiles[[29]][,c(2,6)],listOfFiles[[30]][,c(2,6)],listOfFiles[[31]][,c(2,6)],listOfFiles[[32]][,c(2,6)],
                 listOfFiles[[33]][,c(2,6)],listOfFiles[[34]][,c(2,6)],listOfFiles[[35]][,c(2,6)],listOfFiles[[36]][,c(2,6)],
                 listOfFiles[[37]][,c(2,6)],listOfFiles[[38]][,c(2,6)],listOfFiles[[39]][,c(2,6)],listOfFiles[[40]][,c(2,6)],
                 listOfFiles[[41]][,c(2,6)],listOfFiles[[42]][,c(2,6)],listOfFiles[[43]][,c(2,6)],listOfFiles[[44]][,c(2,6)]
                 )
head(df_vf)
table(df_vf$X.miRNA)
colnames(df_vf)[seq(3,90,2)]<-as.character(summa2$Name)
names(listOfFiles)
summa2$Name
# count number reads in files
df_vft1<-df_vf[,grep("seq",colnames(df_vf),invert=T)]
df_vft2<-df_vft1[,c(-1,-2)]


table(df_vf$X.miRNA)[table(df_vf$X.miRNA)>1]

df_vf<-df_vf[!duplicated(df_vf),]
names(table(paste(df_vf$X.miRNA,df_vf$precursor))[table(paste(df_vf$X.miRNA,df_vf$precursor))>1])
# take out repeted name miRNA
df_vf<-df_vf[df_vf$X.miRNA!="mmu-miR-466i-5p",]
df_vf$X.miRNA[grep('miR171',df_vf$X.miRNA)]
###############################
#format data
# PLANT samples only

colnames(df_vf)[seq(3,90,2)]<-as.character(summa2$Name)



###############################
#format data
# PLANT samples only
plant_samples<-df_vf[,grep("seq.norm",colnames(df_vf),invert=T)]
plant_samples<-plant_samples[,grep("X.",colnames(plant_samples),invert=T)]
plant_samples<-plant_samples[,grep("precursor",colnames(plant_samples),invert=T)]

dim(plant_samples)
rownames(plant_samples)<-paste(df_vf$X.miRNA,df_vf$precursor) #
length(df_vf$X.miRNA)
#table(paste(df_vf$X.miRNA,df_vf$precursor))[table(paste(df_vf$X.miRNA,df_vf$precursor))>1]

table(rownames(plant_samples))[table(rownames(plant_samples))>1]

#filter out microRNA without reads
dim(plant_samples)
##################################################################################################################################



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
### Individual analysis per variety
##################################################################################################################################



###############
# BRA337
###############

#filter out microRNA without reads
plant_samplesF2<-plant_samples

# only work with var BRA337
plant_samplesF2<-plant_samplesF2[,grep("BRA337",colnames(plant_samplesF2),invert=F)]

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.2",
              "CM4574-7.B1.3")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]

# design
AMF_treat<-sapply(strsplit(colnames(plant_samplesF2),"\\."), "[[", 2)


###############################
# Design
colnames(plant_samplesF2)

colData<-  cbind.data.frame(AMF_treat,rep("single-end",length(colnames(plant_samplesF2))) )
colnames(colData)<-c("AMF","type")
rownames(colData)<-colnames(plant_samplesF2)
head(plant_samplesF2)
table(interaction(colData$AMF,colData$PLANT))
###############################
# DEseq object
ddsBRA337 <- DESeqDataSetFromMatrix(countData = plant_samplesF2,
                              colData =colData,
                              design = ~AMF)
ddsBRA337 <- estimateSizeFactors(ddsBRA337)
ddsBRA337@colData
###############################
#Filtering
ddsBRA337 <- ddsBRA337[ rowSums(counts(ddsBRA337)) > 20, ]
rld <- rlog(ddsBRA337, blind=FALSE)
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rownames(colData)#paste( ddsBRA337$AMF, ddsBRA337$PLANT, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(12)
pheatmap(sampleDistMatrix,col=colors)



pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/cmdscale_BRA337.pdf" ,width=10, height=10,useDingbats = F)
p10<-plotPCA(rld, intgroup = c("AMF"),returnData=T)
percentVar <- round(100 * attr(p10, "percentVar"))
ggplot(p10, aes(PC1, PC2, shape=AMF)) + geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16)) 

dev.off()
################################
#differential expression
ddsBRA337 <- DESeq(ddsBRA337)

# results  MOCK  vs B1
res <-results(ddsBRA337,cooksCutoff=F,contrast = c("AMF","MOCK","B1"))
res@listData$padj[res@listData$padj<0.2]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_MOCK_B1.pdf" ,width=10, height=10,useDingbats = F)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

BRA337_M_B1<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(BRA337_M_B1)<-res@rownames
write.table(BRA337_M_B1,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_MOCK_B1.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)


#####################################
# results  MOCK  vs DAOM
res <-results(ddsBRA337,cooksCutoff=F,contrast = c("AMF","MOCK","DAOM197198"))
res@listData$padj[res@listData$padj<0.1]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_MOCK_DAOM.pdf" ,width=10, height=10,useDingbats = F)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

BRA337_M_CAN<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(BRA337_M_CAN)<-res@rownames
write.table(BRA337_M_CAN,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_MOCK_DAOM.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)


#####################################
# results  B1 vs DAOM197198
res <-results(ddsBRA337,cooksCutoff=F,contrast = c("AMF","B1","DAOM197198"))
res@listData$padj[res@listData$padj<0.1]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_B1_DAOM.pdf" ,width=10, height=10,useDingbats = F)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

BRA337_B1_CAN<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(BRA337_B1_CAN)<-res@rownames
write.table(BRA337_B1_CAN,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_BRA337_B1_DAOM.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)


###############
# CM4574-7
###############

#filter out microRNA without reads
plant_samplesF2<-plant_samples

# only work with var CM4574-7
plant_samplesF2<-plant_samplesF2[,grep("CM4574-7",colnames(plant_samplesF2),invert=F)]

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.1",
              "CM4574-7.B1.2")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]

# design
AMF_treat<-sapply(strsplit(colnames(plant_samplesF2),"\\."), "[[", 2)


###############################
# Design
colnames(plant_samplesF2)

colData<-  cbind.data.frame(AMF_treat,rep("single-end",length(colnames(plant_samplesF2))) )
colnames(colData)<-c("AMF","type")
rownames(colData)<-colnames(plant_samplesF2)
head(plant_samplesF2)
table(interaction(colData$AMF,colData$PLANT))
###############################
# DEseq object
ddsCM4574 <- DESeqDataSetFromMatrix(countData = plant_samplesF2,
                              colData =colData,
                              design = ~AMF)
ddsCM4574 <- estimateSizeFactors(ddsCM4574)
ddsCM4574@colData
###############################
#Filtering
ddsCM4574 <- ddsCM4574[ rowSums(counts(ddsCM4574)) > 20, ]
rld <- rlog(ddsCM4574, blind=FALSE)
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rownames(colData)#paste( ddsCM4574$AMF, ddsCM4574$PLANT, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(12)
pheatmap(sampleDistMatrix,col=colors)


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/cmdscale_CM4574-7.pdf" ,width=10, height=10,useDingbats = F)
p10<-plotPCA(rld, intgroup = c("AMF"),returnData=T)
percentVar <- round(100 * attr(p10, "percentVar"))
ggplot(p10, aes(PC1, PC2, shape=AMF)) + geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
theme(plot.title = element_text(size = 20, face = "bold"),
      text = element_text(size = 18),
      axis.title = element_text(face="bold"),
      axis.text.x=element_text(size = 16)) 

dev.off()
################################
#differential expression
ddsCM4574 <- DESeq(ddsCM4574)

# results  MOCK  vs B1
res <-results(ddsCM4574,cooksCutoff=F,contrast = c("AMF","MOCK","B1"))
res@listData$padj[res@listData$padj<0.5]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)
res@listData$pvalue[grep("miR171",res@rownames)]
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_MOCK_B1.pdf" ,width=10, height=10,useDingbats = F)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

CM4574_M_B1<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(CM4574_M_B1)<-res@rownames
write.table(CM4574_M_B1,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_MOCK_B1.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)


#####################################
# results  MOCK  vs DAOM
res <-results(ddsCM4574,cooksCutoff=F,contrast = c("AMF","MOCK","DAOM197198"))
res@listData$padj[res@listData$padj<0.1]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_MOCK_DAOM.pdf" ,width=10, height=10,useDingbats = F)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

CM4574_M_CAN<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(CM4574_M_CAN)<-res@rownames
write.table(CM4574_M_CAN,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_MOCK_DAOM.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)


#####################################
# results  B1 vs DAOM197198
res <-results(ddsCM4574,cooksCutoff=F,contrast = c("AMF","B1","DAOM197198"))
res@listData$padj[res@listData$padj<0.1]
topGene <- rownames(res)[res@listData$padj<0.1]
topGene <-na.omit(topGene)

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_B1_DAOM.pdf" ,width=10, height=10,useDingbats = F,)

plotMA(res, main="DESeq2 MOCK vs B1", ylim=c(-2,2),las=1)

for ( i in 1:length(topGene)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=topGene[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(topGene[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  tryCatch(print(p10),error=function(e) 1 )
}

dev.off()

CM4574_B1_CAN<-cbind.data.frame(res@listData$baseMean,res@listData$log2FoldChange,res@listData$lfcSE,
                             res@listData$stat,res@listData$pvalue,res@listData$padj)
rownames(CM4574_B1_CAN)<-res@rownames
write.table(CM4574_B1_CAN,"~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_CM4574-7_B1_DAOM.txt",quote=FALSE,sep='\t',
            col.names = T, row.names = T)

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
# do venn diagram to see per variety how many isolates shared. then do this for tree vars.
Venn_BRA337_M_B1<-na.omit(rownames(BRA337_M_B1)[BRA337_M_B1$`res@listData$padj`<0.1])
Venn_BRA337_M_CAN<-na.omit(rownames(BRA337_M_CAN)[BRA337_M_CAN$`res@listData$padj`<0.1])
Venn_COL2215_M_B1<-na.omit(rownames(COL2215_M_B1)[COL2215_M_B1$`res@listData$padj`<0.1])
Venn_COL2215_M_CAN<-na.omit(rownames(COL2215_M_CAN)[COL2215_M_CAN$`res@listData$padj`<0.1])
Venn_CM4574_M_B1<-na.omit(rownames(CM4574_M_B1)[CM4574_M_B1$`res@listData$padj`<0.1])
Venn_CM4574_M_CAN<-na.omit(rownames(CM4574_M_CAN)[CM4574_M_CAN$`res@listData$padj`<0.1])

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")
geneLists<-list()
geneLists[[1]]<-as.character(Venn_BRA337_M_B1)
geneLists[[2]]<-as.character(Venn_BRA337_M_CAN)
geneLists[[3]]<-as.character(Venn_COL2215_M_B1)
geneLists[[4]]<-as.character(Venn_COL2215_M_CAN)
geneLists[[5]]<-as.character(Venn_CM4574_M_B1)
geneLists[[6]]<-as.character(Venn_CM4574_M_CAN)
# We can rename our list vectors
names(geneLists) <- c("BRA337_M_B1", "BRA337_M_CAN","COL2215_M_B1", "COL2215_M_CAN",
                      "CM4574_M_B1", "CM4574_M_CAN")


# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_all_BRA337.pdf" ,width=10, height=10,useDingbats = F,)
venn.plot <- venn.diagram(geneLists[1:2] , NULL, fill=c("grey33", "grey75"), cat.cex=4,height = 2,width = 2,
                          alpha=c(0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                          category.names=c("BRA337_M_B1", "BRA337_M_DAOM197198"), main="miRNA DE in BRA337")
grid.draw(venn.plot)
dev.off()
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_all_COL2215.pdf" ,width=10, height=10,useDingbats = F,)

venn.plot <- venn.diagram(geneLists[3:4]  ,NULL, fill=c("grey33", "grey75"), cat.cex=4,height = 2,width = 2,
                          alpha=c(0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                          category.names=c("COL2215_M_B1", "COL2215_M_DAOM197198"), main="miRNA DE in COL2215")
grid.draw(venn.plot)
dev.off()
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_all_CM4574.pdf" ,width=10, height=10,useDingbats = F,)

venn.plot <- venn.diagram(geneLists[5:6]  , NULL, fill=c("grey33", "grey75"), cat.cex=4,height = 2,width = 2,
                          alpha=c(0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                          category.names=c("CM4574-7_M_B1", "CM4574-7_M_DAOM197198"), main="miRNA DE in CM4574-7")
grid.draw(venn.plot)
dev.off()
# To get the list of gene present in each Venn compartment we can use the gplots package
require("gplots")
a <- venn(geneLists, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(a)
# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")
# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head) 

### ACross varieties
geneListsV<-list()
geneListsV[[1]]<-unique(c(as.character(Venn_BRA337_M_B1),as.character(Venn_BRA337_M_CAN)))
geneListsV[[2]]<-unique(c(as.character(Venn_COL2215_M_B1),as.character(Venn_COL2215_M_CAN)))
geneListsV[[3]]<-unique(c(as.character(Venn_CM4574_M_B1),as.character(Venn_CM4574_M_CAN)))

geneListsV[[1]]<-as.character(Venn_BRA337_M_B1)
geneListsV[[2]]<-as.character(Venn_CM4574_M_B1)
geneListsV[[3]]<-as.character(Venn_BRA337_M_CAN)
geneListsV[[4]]<-as.character(Venn_CM4574_M_CAN)


# We can rename our list vectors
names(geneListsV) <- c("BRA337","COL2215","CM4574")
names(geneListsV) <- c("BRA337_MOCK_B1","CM4574_MOCK_B1",
                       "BRA337_MOCK_DAOM","CM4574_MOCK_DAOM")



pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_all_cultivars.pdf" ,width=10, height=10,useDingbats = F,)
venn.plotV <- venn.diagram(geneListsV , NULL, fill=c("grey33","grey75","grey33","grey75"), cat.cex=4,height = 2,width = 2,
                           alpha=c(0.5,0.5,0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                           category.names=c("BRA337_M_B1","CM4574-7_M_B1","BRA337_M_DAOM","CM4574-7_M_DAOM"), main="miRNA DE shared in cultivars")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram

grid.draw(venn.plotV)
dev.off()


# only miRNA members across varieties
geneListsV<-list()
geneListsV[[1]]<-unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_BRA337_M_B1)))," "),
              function (x) x[2])),"-"), "[[", 2))
geneListsV[[2]]<-unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_CM4574_M_B1)))," "),
                                     function (x) x[2])),"-"), "[[", 2))
geneListsV[[3]]<-unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_BRA337_M_CAN)))," "),
                                                      function (x) x[2])),"-"), "[[", 2))
geneListsV[[4]]<-unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_CM4574_M_CAN)))," "),
                                                      function (x) x[2])),"-"), "[[", 2))

names(geneListsV) <- c("BRA337_MOCK_B1","CM4574_MOCK_B1",
                       "BRA337_MOCK_DAOM","CM4574_MOCK_DAOM")

geneListsV
# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_members_cultivars.pdf" ,width=10, height=10,useDingbats = F,)
venn.plotV <- venn.diagram(geneListsV , NULL, fill=c("grey33","grey75","grey33","grey75"), cat.cex=4,height = 2,width = 2,
                           alpha=c(0.5,0.5,0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                           category.names=c("BRA337_M_B1","CM4574-7_M_B1","BRA337_M_DAOM","CM4574-7_M_DAOM"), main="miRNA DE shared in cultivars")

grid.draw(venn.plotV)
dev.off()
# only miRNA families across varieties
geneListsV<-list()
geneListsV[[1]]<-unique(gsub("[a-z]","",
                             unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_BRA337_M_B1)))," "),
                                                                  function (x) x[2])),"-"), "[[", 2))
))

geneListsV[[2]]<-unique(gsub("[a-z]","",
                               unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_CM4574_M_B1)))," "),
                                                                    function (x) x[2])),"-"), "[[", 2))
))
geneListsV[[3]]<-unique(gsub("[a-z]","",
                             unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_BRA337_M_CAN)))," "),
                                                                  function (x) x[2])),"-"), "[[", 2))
))

geneListsV[[4]]<-unique(gsub("[a-z]","",
                             unique(sapply(strsplit(unlist(lapply(strsplit(unique(c(as.character(Venn_CM4574_M_CAN)))," "),
                                                                  function (x) x[2])),"-"), "[[", 2))
))
names(geneListsV) <- c("BRA337_MOCK_B1","CM4574_MOCK_B1",
                       "BRA337_MOCK_DAOM","CM4574_MOCK_DAOM")

geneListsV[1:4]

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_families_cultivars.pdf" ,width=10, height=10,useDingbats = F,)
venn.plotV <- venn.diagram(geneListsV , NULL, fill=c("grey33","grey75","grey33","grey75"), cat.cex=4,height = 2,width = 2,
                           alpha=c(0.5,0.5,0.5,0.5), cex = 4, cat.fontface=4, main.cex= 4,sub.cex=4,
                           category.names=c("BRA337_M_B1","CM4574-7_M_B1","BRA337_M_DAOM","CM4574-7_M_DAOM"), main="miRNA DE shared in cultivars")

grid.draw(venn.plotV)
dev.off()

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
###### Super heat map. ploting different miRNA x samples + plant DW
#devtools::install_github("rlbarter/superheat")
library(superheat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)
library(RColorBrewer)
#https://rlbarter.github.io/superheat-examples/fMRI/


# do clustering

plant_samplesF2<-plant_samples

# only work with var CM4574-7
plant_samplesF2<-plant_samplesF2[,grep("BRA337",colnames(plant_samplesF2),invert=F)]

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.1",
              "CM4574-7.B1.2")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]

GENE2SUPERHEAT<-plant_samplesF2[grep( "miR171", rownames(plant_samplesF2)),]
GENE2SUPERHEAT<-GENE2SUPERHEAT[ rowSums(GENE2SUPERHEAT) > 20, ]

samples.clusters <-sapply(strsplit(colnames(plant_samplesF2[1:5,]),"\\."), "[[", 2)
Ext<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/External_data_miRNA_samples_v1.txt",h=T)

External_data<-Ext[Ext$PLANT=="BRA337",]
External_data<-External_data[match(colnames(plant_samplesF2), External_data$R_NAMES),]

png("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miR171_DW_BRA337.png" ,width=1000, height=1000,)

superheat(GENE2SUPERHEAT, 
          #X.text =as.matrix(GENE2SUPERHEAT),
          yt =External_data$DW_total,
          yt.axis.name = "Plant total\n dry weight \n (gms)",
          X.text.angle=0,
          heat.pal = brewer.pal(6, "GnBu"),
         # heat.pal.values=c(0,0.1,0.2,0.3,0.5,1),
          yt.num.ticks=5,
          yt.axis.size=20,
          yt.axis.name.size=20,
          yt.obs.col = rep("slategray4", 9),
          #yt.point.alpha = 0.6,
          #yt.axis.name.size = 10,
          #yt.plot.size = 0.7,
          yt.point.size = 2,
          left.label.text.size=6,
          bottom.label.text.size=6,
          left.label.size=0.6,
          #bottom.label.size=0.5,
          yt.plot.type="bar",
          #bottom.label.text.angle=45,
          #left.label.text.angle=45,
          #membership.rows = miRNA.clusters,
          membership.cols = samples.clusters,
          clustering.method ="hierarchical",
          row.dendrogram=TRUE,
          
          #left.label = "none",
          #bottom.label = "none",
          grid.hline.col = "white",
          grid.vline.col = "white",
          #grid.hline.size = 2,
          #grid.vline.size = 2,
          #row.title = "Validation images (120)",
          #column.title = "Voxels (1,294)",
          #title = "(a)"
          )
dev.off()


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

### mirNA correlation to pheno
plant_samplesF2<-plant_samples

# only work with var CM4574-7

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.1",
              "CM4574-7.B1.2")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]

GENE2SUPERHEAT<-plant_samplesF2#[grep( "miR171", rownames(plant_samplesF2)),]
GENE2SUPERHEAT<-GENE2SUPERHEAT[ rowSums(GENE2SUPERHEAT) > 20, ]
head(GENE2SUPERHEAT)

Ext<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/External_data_miRNA_samples_v1.txt",h=T)
External_data<-Ext
External_data<-External_data[match(colnames(plant_samplesF2), External_data$R_NAMES),]
par(mfrow=c(1,1))
cor(t(GENE2SUPERHEAT[grep("COL2215",colnames(GENE2SUPERHEAT))]),External_data$DW_total[External_data$PLANT=="COL2215"])
cor(t(GENE2SUPERHEAT[grep("BRA337",colnames(GENE2SUPERHEAT))]),External_data$DW_total[External_data$PLANT=="BRA337"])
cor(t(GENE2SUPERHEAT[grep("CM4574",colnames(GENE2SUPERHEAT))]),External_data$DW_total[External_data$PLANT=="CM4574"])

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
# resume data in varieties .

# MOCK vs CAN
dim(COL2215_M_CAN[COL2215_M_CAN[,6]<0.1,])
dim(BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,])
dim(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,])

table(rbind(COL2215_M_CAN[COL2215_M_CAN[,6]<0.1,],
        BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,],
        CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,]))
# MOCK vs B1
dim(COL2215_M_B1[COL2215_M_B1[,6]<0.1,])
dim(BRA337_M_B1[BRA337_M_B1[,6]<0.1,])
dim(CM4574_M_B1[CM4574_M_B1[,6]<0.1,])

# CAN vs B1
dim(COL2215_B1_CAN[COL2215_B1_CAN[,6]<0.1,])
dim(BRA337_B1_CAN[BRA337_B1_CAN[,6]<0.1,])
dim(CM4574_B1_CAN[CM4574_B1_CAN[,6]<0.1,])

##############################################################
##############################################################

table(c(all_sign_miRNA,col))
col<-unique(c(rownames(COL2215_M_CAN[COL2215_M_CAN[,6]<0.1,]),rownames(COL2215_M_B1[COL2215_M_B1[,6]<0.1,])))

all_sign_miRNA<-unique(c(rownames(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,]),rownames(CM4574_M_B1[CM4574_M_B1[,6]<0.1,]),
               rownames(na.omit(BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,])),rownames(na.omit(BRA337_M_B1[BRA337_M_B1[,6]<0.1,]))))

all_sign_miRNA<-sapply(strsplit(all_sign_miRNA," "), "[[", 1)

sapply(strsplit(rownames(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,])," "), "[[", 1)
sapply(strsplit(rownames(CM4574_M_B1[CM4574_M_B1[,6]<0.1,])," "), "[[", 1)

sapply(strsplit(rownames(na.omit(BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,]))," "), "[[", 1)
sapply(strsplit(rownames(na.omit(BRA337_M_B1[BRA337_M_B1[,6]<0.1,]))," "), "[[", 1)

###########################################################################
###########################################################################
#### different families that are significant
table(sapply(strsplit(all_sign_miRNA,"-"), "[[", 2))

miRfamilies<-c("miR156","miR159","miR160","miR167","miR168",
  "miR171","miR172","miR396","miR397","miR398","miR408","miR6445")

atha_significant<-all_sign_miRNA[grep("ath",all_sign_miRNA)]

###########################################################################

# extract names from fasta xargs faidx -d "|" My_seq.fasta < My_Ids.txt
all_sign_miRNA<-unique(c(rownames(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,]),rownames(CM4574_M_B1[CM4574_M_B1[,6]<0.1,]),
                         rownames(na.omit(BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,])),rownames(na.omit(BRA337_M_B1[BRA337_M_B1[,6]<0.1,]))))
all_sign_miRNA
# barplot 3 varieties, 3 treatments

library(RColorBrewer)


geneXX <-"vvi-miR160b vvi-MIR160b"
geneCountsB <- plotCounts(ddsBRA337, gene=geneXX , intgroup=c("AMF"), returnData=TRUE)
geneCountsC <- plotCounts(ddsCM4574, gene=geneXX , intgroup=c("AMF"), returnData=TRUE)
geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
colnames(geneCounts)[3]<- "PLANT"

pdf(paste("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/",sapply(strsplit(geneXX," "), "[[", 1),"_2var.pdf",sep="") ,width=6, height=6,useDingbats = F)
p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Normalized counts") +
  scale_x_discrete(name = "Plant cultivar") +
  ggtitle(sapply(strsplit(geneXX," "), "[[", 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16)) 
p10 = p10 + scale_fill_grey(start = 0, end = .9)
p10
dev.off()
################################################################
geneXX <-all_sign_miRNA[grep("miR408",all_sign_miRNA)]

sapply(strsplit(geneXX," "), "[[", 1)
for ( i in 1:length(geneXX)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=geneXX[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=geneXX[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(geneXX[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) +
    scale_fill_brewer(palette = "Spectral")
  tryCatch(print(p10),error=function(e) 1 )
  #tryCatch(ggsave("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA171_noF.pdf" ,width=10, height=10,useDingbats = F,p10),error=function(e) 1 )
}
################################################################
################################################################
################################################################
# barplot 3 varieties, 3 treatments, all mir171
geneXX<-

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA171_noF.pdf" ,width=10, height=10,useDingbats = F)
for ( i in 1:length(geneXX)) {
  # barplot 3 varieties, 3 treatments, all mir171
  geneCountsB <- plotCounts(ddsBRA337, gene=geneXX[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCountsC <- plotCounts(ddsCM4574, gene=geneXX[i] , intgroup=c("AMF"), returnData=TRUE)
  geneCounts<-cbind(rbind(geneCountsA,geneCountsB,geneCountsC),rep(c("BRA337","CM4574_7"),each=9))
  colnames(geneCounts)[3]<- "PLANT"
  p10 <- ggplot(geneCounts, aes(x = PLANT, y = count, fill = AMF)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(geneXX[i]) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) +
    scale_fill_brewer(palette = "Spectral")
  tryCatch(print(p10),error=function(e) 1 )
  #tryCatch(ggsave("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA171_noF.pdf" ,width=10, height=10,useDingbats = F,p10),error=function(e) 1 )
}
dev.off()
#################################################################################################################################
#################################################################################################################################
################################################################################################################################

#################################################################################################################################
############ make phylogenetic three for miRNA x
################################################################################################################################

# procedure in bash, web, cluster to make phylogenetic tree.
# 1) select miRNA in highfidelity mature file. grep -A1 'miRXXX' mature_hig... high_conf_mature.fasta > miRXX_HF.fasta 
# 2) format fasta in textwrangler, delete "--" lines
# 3) align them in  http://www.ebi.ac.uk/Tools/msa/muscle/ using muscle aligner, fasta output
# 4) rename file miRXXX_HF_align.fasta  and format name with only mirna notation
# 5) convert to .nex http://users-birc.au.dk/biopv/php/fabox/fasta2mrbayes.php unix output
# 6) then use mr bayes to calculate the tree. in cluster
# module add Phylogeny/mrbayes/3.2.6;
# mb
# execute miRNAXXX_HF_align.nex 
# lset nst=2
# mcmc ngen=1000000 samplefreq=25 printfreq=100 diagnfreq=1000 temp=0.25
# sump
# sumt Conformat=simple Outputname=SimpleTree_miRXXX_HF_.nex
# 7) import alligned file into R
#install.packages(c("ape","phangorn","seqinr"))
library(ape)
library(phangorn)
library(seqinr)
# import mr bayes phylogeny and check coinoculation
#install.packages("phangorn")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
library(devtools)

# plot mr bayes with different color to significant miRNA


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_phylogeny/miR408_tree.pdf" ,width=13, height=13,useDingbats = F) 
mrbayes.tree<-read.nexus("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_phylogeny/SimpleTree_miR156_HF.nex.con.tre") 
color<-gsub("_","-",mrbayes.tree$con_50_majrule$tip.label) %in% all_sign_miRNA
color[color=="FALSE"]<-"grey"
color[color=="TRUE"]<-"black"
p <- character(length(mrbayes.tree[[1]]$node.label))
#The following three lines define your labeling scheme. 
p[mrbayes.tree[[1]]$node.label >= 0.95] <- "black"
p[mrbayes.tree[[1]]$node.label < 0.95 & mrbayes.tree[[1]]$node.label >= 0.75] <- "gray"
p[mrbayes.tree[[1]]$node.label < 0.75] <- "white"

plot(mrbayes.tree[[1]],tip.color = color,cex=1.4,type="fan")
nodelabels(pch=21, cex = 2.7, bg = p)
dev.off()

############################
##### info about all significant miRNA sequences

## tree
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_phylogeny/miRNA408_tree_v2.pdf" ,width=13, height=13,useDingbats = F) 
mrbayes.tree<-read.nexus("/Users/ivadako/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_miRNA/miRNA_phylogeny/SimpleTree_miR408_sign.nex.con.tre") 
color<-gsub("_","-",mrbayes.tree$con_50_majrule$tip.label) %in% all_sign_miRNA
color[color=="FALSE"]<-"grey"
color[color=="TRUE"]<-"black"
p <- character(length(mrbayes.tree[[1]]$node.label))
#The following three lines define your labeling scheme. 
p[mrbayes.tree[[1]]$node.label >= 0.95] <- "black"
p[mrbayes.tree[[1]]$node.label < 0.95 & mrbayes.tree[[1]]$node.label >= 0.7] <- "gray"
p[mrbayes.tree[[1]]$node.label < 0.7] <- "white"

plot(mrbayes.tree[[1]],tip.color = color,cex=1.4,type="fan")
nodelabels(pch=21, cex = 2.7, bg = p)
dev.off()

##  regroupe similar miRNA 
all_sign_miRNA



#################################################################################################################################
#################################################################################################################################
################################################################################################################################

###############################
#External data info
Ext<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/External_data_miRNA_samples_v1.txt",h=T)
head(Ext)
options(scipen=999)
Ext<-Ext[grep("COL2215",Ext$Name,invert=T),]
head(Ext)
barplot(Ext$SEQ_nb_mirna,names.arg = Ext$Name,las=2)
boxplot(Ext$microRNA_10reads~interaction(Ext$AMF,Ext$PLANT),las=2)
boxplot(Ext$microRNA_10reads~Ext$Lane,las=2)
barplot(Ext$microRNA_10reads,names.arg = Ext$Name,las=2)
par(mfrow=c(1,1))
boxplot(Ext$microRNA_10reads[Ext$PLANT=="BRA337"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT=="BRA337"&Ext$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
boxplot(Ext$microRNA_10reads[Ext$PLANT=="COL2215"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT=="COL2215"&Ext$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
boxplot(Ext$microRNA_10reads[Ext$PLANT!="BRA337"&Ext$PLANT!="COL2215"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT!="BRA337"&Ext$PLANT!="COL2215"&Ext$AMF!="CO_INOC"],las=2,ylim=c(0,1000))
kruskal.test(Ext$microRNA_10reads[Ext$PLANT=="BRA337"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT=="BRA337"&Ext$AMF!="CO_INOC"])
kruskal.test(Ext$microRNA_10reads[Ext$PLANT=="COL2215"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT=="COL2215"&Ext$AMF!="CO_INOC"])
kruskal.test(Ext$microRNA_10reads[Ext$PLANT!="BRA337"&Ext$PLANT!="COL2215"&Ext$AMF!="CO_INOC"]~Ext$AMF[Ext$PLANT!="BRA337"&Ext$PLANT!="COL2215"&Ext$AMF!="CO_INOC"])


# COLONIZATION per treatment
boxplot(Ext$COL_perc~interaction(Ext$TREAT,Ext$PLANT),las=2)
boxplot(Ext$COL_perc[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"],las=2,ylim=c(0,1))
kruskal.test(Ext$COL_perc[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"])
boxplot(Ext$COL_perc[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"],las=2,ylim=c(0,1))
kruskal.test(Ext$COL_perc[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"])
boxplot(Ext$COL_perc[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"],las=2,ylim=c(0,1))
kruskal.test(Ext$COL_perc[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"])


# TOTAL DW per treatment
  boxplot(Ext$DW_total~interaction(Ext$TREAT,Ext$PLANT),las=2)
  boxplot(Ext$DW_total[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"],las=2,ylim=c(0,20))
  kruskal.test(Ext$DW_total[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"])
  boxplot(Ext$DW_total[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"],las=2,ylim=c(0,20))
  kruskal.test(Ext$DW_total[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"])
  boxplot(Ext$DW_total[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"],las=2,ylim=c(0,15))
  kruskal.test(Ext$DW_total[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"])
  
  # nb miRNA per treatment
boxplot(Ext$SEQ_nb_mirna~interaction(Ext$TREAT,Ext$PLANT),las=2,ylim=c(0,380))
boxplot(Ext$SEQ_nb_mirna[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_nb_mirna[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"])
boxplot(Ext$SEQ_nb_mirna[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_nb_mirna[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"])
boxplot(Ext$SEQ_nb_mirna[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_nb_mirna[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"])

# mean read per miRNA per treatment
boxplot(Ext$SEQ_Mean_read_miRNA~interaction(Ext$TREAT,Ext$PLANT),las=2,ylim=c(0,850))
boxplot(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="COL2215"]~Ext$TREAT[Ext$PLANT=="COL2215"])
boxplot(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="BRA337"]~Ext$TREAT[Ext$PLANT=="BRA337"])
boxplot(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"],las=2,ylim=c(0,380))
kruskal.test(Ext$SEQ_Mean_read_miRNA[Ext$PLANT=="CM4574"]~Ext$TREAT[Ext$PLANT=="CM4574"])



# COLONIZATION PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_COL_perc.pdf" ,width=10, height=10,useDingbats = F)
COL_perc_plot<-Ext[,colnames(Ext) %in% c("COL_perc","PLANT","TREAT")]
  p10 <- ggplot(COL_perc_plot, aes(x = PLANT, y = COL_perc, fill = TREAT)) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Colonization (% of colonized roots)") +
    scale_x_discrete(name = "Plant variety") +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16))
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
  print(p10)
dev.off()

# DW TOTAL PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_DW_total.pdf" ,width=10, height=10,useDingbats = F)
DW_total_plot<-Ext[,colnames(Ext) %in% c("DW_total","PLANT","TREAT")]
p10 <- ggplot(DW_total_plot, aes(x = PLANT, y = DW_total, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Plant total dry weigth (gms)") +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16))
p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()

# DW AG PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_DW_AG.pdf" ,width=10, height=10,useDingbats = F)
DW_AG_plot<-Ext[,colnames(Ext) %in% c("DW_AG","PLANT","TREAT")]
p10 <- ggplot(DW_AG_plot, aes(x = PLANT, y = DW_AG, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Plant AG dry weigth (gms)") +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16))
p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()

# DW UG PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_DW_UG.pdf" ,width=10, height=10,useDingbats = F)
DW_UG_plot<-Ext[,colnames(Ext) %in% c("DW_UG","PLANT","TREAT")]
p10 <- ggplot(DW_UG_plot, aes(x = PLANT, y = DW_UG, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Plant UG dry weigth (gms)") +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16))
p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()

# Height_125dai PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_Height_125dai.pdf" ,width=10, height=10,useDingbats = F)
Height_125dai_plot<-Ext[,colnames(Ext) %in% c("Height_125dai","PLANT","TREAT")]
p10 <- ggplot(Height_125dai_plot, aes(x = PLANT, y = Height_125dai, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Plant Height 125dai") +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16))
p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()


# SEQ_nb_mirna PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_SEQ_nb_mirna.pdf" ,width=10, height=10,useDingbats = F)
SEQ_nb_mirna_plot<-Ext[,colnames(Ext) %in% c("SEQ_nb_mirna","PLANT","TREAT")]
p10 <- ggplot(SEQ_nb_mirna_plot, aes(x = PLANT, y = SEQ_nb_mirna, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Nb. miRNA",limits=c(0,400)) +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()

# SEQ_Mean_read_miRNA PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_samples_SEQ_Mean_read_miRNA.pdf" ,width=10, height=10,useDingbats = F)
SEQ_Mean_read_miRNA_plot<-Ext[,colnames(Ext) %in% c("SEQ_Mean_read_miRNA","PLANT","TREAT")]
p10 <- ggplot(SEQ_Mean_read_miRNA_plot, aes(x = PLANT, y = SEQ_Mean_read_miRNA, fill = TREAT)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Mean cov. miRNA",limits=c(0,400)) +
  scale_x_discrete(name = "Plant variety") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16))
p10 = p10 + scale_fill_grey(start = 0, end = .9)
print(p10)
dev.off()

# SEQ_Mean_read_miRNA PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_collapsedfa_reads.pdf" ,width=10, height=10,useDingbats = F)
SEQ_Mean_read_miRNA_plot<-Ext[,colnames(Ext) %in% c("collapsed.fa","SEQ_nb_mirna","PLANT","TREAT")]
p10 <- ggplot(SEQ_Mean_read_miRNA_plot, aes(y = SEQ_nb_mirna, x = collapsed.fa,color=PLANT)) +
  geom_point(shape=16,size=6) +
  scale_y_continuous(name = " Nb. miRNA ",limits=c(0,400)) +
  scale_x_continuous(name = "Reads in file",limits=c(0,25311120)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16)) 
p10 = p10 + scale_color_grey(start = 0, end = .9)
print(p10)
dev.off()

# SEQ_Mean_read_miRNA PLOT
pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/miRNA_collapsedfa_reads2.pdf" ,width=10, height=10,useDingbats = F)
SEQ_Mean_read_miRNA_plot<-Ext[,colnames(Ext) %in% c("collapsed.fa","SEQ_nb_mirna","PLANT","TREAT")]
p10 <- ggplot(SEQ_Mean_read_miRNA_plot, aes(y = SEQ_nb_mirna, x = collapsed.fa,color=TREAT)) +
  geom_point(shape=16,size=6) +
  scale_y_continuous(name = " Nb. miRNA ",limits=c(0,400)) +
  scale_x_continuous(name = "Reads in file",limits=c(0,25311120)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 16)) 
  p10 = p10 + scale_color_grey(start = 0, end = .9)
print(p10)
dev.off()

