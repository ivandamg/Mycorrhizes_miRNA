########################################################################################
#' Import transcript-level abundances and estimated counts for gene-level analysis packages
########################################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
#install.packages("pheatmap")
#biocLite("topGO")
#biocLite("limma")
#biocLite("edgeR")
library("gplots")
library("RColorBrewer")
library(Rgraphviz)
library(topGO)
library(edgeR)
library(limma)
library(heatmap2)
library(dynamicTreeCut)
library(WGCNA)
library(pheatmap)
options("scipen"=3,"digits"=2)
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/WGCNA_1.49.tgz",type = "source", repos = NULL, lib=.Library) 

########################################################################################
# ANALYSIS AT GENE LEVEL 
# input data: Kallisto pseudo-alignment. 

tximport <- function(files,
                     type=c("none","kallisto","salmon","sailfish","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     reader=read.delim,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer,
                     collatedFiles,
                     ignoreTxVersion=FALSE) {
  
  type <- match.arg(type, c("none","kallisto","salmon","sailfish","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment='#') 
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    importer <- reader
  }
  
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      
      #####################################################################
      # some temporary code for detecting older fishes
      if ((i == 1) &
          (type %in% c("salmon","sailfish")) &
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length" 
        # because the comment lines have the same comment character
        # as the header, need to name the column names
        importer <- function(x) {
          tmp <- reader(x, comment="#")
          names(tmp) <- c("Name","Length","TPM","NumReads")
          tmp
        }
        # re-read the first file
        raw <- as.data.frame(importer(files[i]))
      }
      #####################################################################
      
      # does the table contain gene association or was an external tx2gene table provided?
      if (is.null(tx2gene) & !txOut) {
        # e.g. Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    message("")
    
    txi <- list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx,
                countsFromAbundance="no")
    
    # if the user requested just the transcript-level data:
    if (txOut) {
      return(txi)
    }
    
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
    # e.g. RSEM already has gene-level summaries
    # just combine the gene-level summaries across files
  } else {
    # stating the obvious:
    if (txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
    
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}

# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize tx-level matrices to gene-level
#' @export
summarizeToGene <- function(txi,
                            tx2gene,
                            ignoreTxVersion=FALSE,
                            countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
) {
  
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- txi$abundance
  countsMatTx <- txi$counts
  lengthMatTx <- txi$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {
    colnames(tx2gene) <- c("tx","gene")
    if (ignoreTxVersion) {
      txId <- sapply(strsplit(as.character(txId), "\\."), "[[", 1)
    }
    tx2gene$gene <- factor(tx2gene$gene)
    tx2gene$tx <- factor(tx2gene$tx)
    # remove transcripts (and genes) not in the abundances
    tx2gene <- tx2gene[tx2gene$tx %in% txId,]
    tx2gene$gene <- droplevels(tx2gene$gene)
    ntxmissing <- sum(!txId %in% tx2gene$tx)
    if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)
    sub.idx <- txId %in% tx2gene$tx
    abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
    countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
    lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
    txId <- txId[sub.idx]
    geneId <- tx2gene$gene[match(txId, tx2gene$tx)]
  }
  
  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
  message("summarizing counts")
  countsMat <- fastby(countsMatTx, geneId, colSums)
  message("summarizing length")
  
  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
  lengthMat <- weightedLength / abundanceMat   
  
  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)
  
  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))
  
  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)
  
  if (countsFromAbundance != "no") {
    countsSum <- colSums(countsMat)
    if (countsFromAbundance == "lengthScaledTPM") {
      newCounts <- abundanceMat * rowMeans(lengthMat)
    } else {
      newCounts <- abundanceMat
    }
    newSum <- colSums(newCounts)
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  if (ncol(m) > 1) {
    t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
  } else {
    matrix(sapply(idx, function(i) fun(m[i,,drop=FALSE])),
           dimnames=list(levels(f), colnames(m)))
  }
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
        # otherwise use the geometric mean of the lengths from the other samples
        idx <- is.nan(lengthMat[i,])
        lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
      }
    }
  }
  lengthMat
}

#####################################################################################################################################################
#####################################################################################################################################################
# transfrom transcripts to genes & import data
########################################################################################################################################
#Mesculenta
setwd('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/mRNA_on_samples/') 

filesToProcess <- dir(pattern = "*_Cassava_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_Cassava_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/mRNA_on_samples/V4_10_Cassava_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')
txi_YUCA <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_YUCA<-txi_YUCA[[2]]



########################################################################################################################################
########################################################################################################################################
# CASSAVA
########################################################################################################################################
########################################################################################################################################
# DATA SELECTION
########################################################################################################################################
# rename samples
summa<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/summary_seq_results_v1.txt",h=T)
head(summa)
summa2<-summa[summa$Discard=="No",]
summa2<-na.omit(summa2)


# Selection all but co-inoculation
colnames(FOUR_VARS_YUCA)
colnames(df_vf)[seq(3,90,2)]<-as.character(summa2$Name)

to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.2",
              "CM4574-7.B1.3")


FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("_3"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("_7"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("_11"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("_15"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("_19"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("CANB1"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_16b"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8"), colnames(FOUR_VARS_YUCA),invert =T)]

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
# DGE
## PLOT GENE EXPRESSION OF TARGET GENES

VAR<-c("CM4574","COL2215","COL2215","COL2215","COL2215","COL2215","COL2215","COL2215","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574")
TREAT<-c("DAOM197198","B1","MOCK","DAOM197198","B1","MOCK","MOCK","DAOM197198","DAOM197198","B1","MOCK","V5_14","V5_16","V5_17","B1","V5_20","MOCK","B1","MOCK","DAOM197198","B1","MOCK","V6_13","V6_14","MOCK","DAOM197198","B1","MOCK","DAOM197198")
REP<-c("2","1","1","1","2","3","5","3","2","1","1","V5_14","V5_16","V5_17","2","V5_20","4","3","5","4","1","1","V6_13","V6_14","2","3","4","4","4")
NAME<-c("CM4574.DAOM197198.2","COL2215.B1.1","COL2215.MOCK.1","COL2215.DAOM197198.1","COL2215.B1.2","COL2215.MOCK.3","COL2215.MOCK.5","COL2215.DAOM197198.3","BRA337.DAOM197198.2","BRA337.B1.1","BRA337.MOCK.1","BRA337.V5_14.V5_14","BRA337.V5_16.V5_16","BRA337.V5_17.V5_17","BRA337.B1.2","BRA337.V5_20.V5_20","BRA337.MOCK.4","BRA337.B1.3","BRA337.MOCK.5","BRA337.DAOM197198.4","CM4574.B1.1","CM4574.MOCK.1","CM4574.V6_13.V6_13","CM4574.V6_14.V6_14","CM4574.MOCK.2","CM4574.DAOM197198.3","CM4574.B1.4","CM4574.MOCK.4","CM4574.DAOM197198.4")

TREAT_YUCA<-factor(TREAT,levels=c("MOCK","DAOM197198","B1"))
colnames(FOUR_VARS_YUCA)<-NAME
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("COL2215.MOCK.5"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("COL2215"), colnames(FOUR_VARS_YUCA),invert =T)]

REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_YUCA),"\\."),function (x) x[[3]]))
VAR<-unlist(lapply(strsplit(colnames(FOUR_VARS_YUCA),"\\."),function (x) x[[1]]))
TREAT<-unlist(lapply(strsplit(colnames(FOUR_VARS_YUCA),"\\."),function (x) x[[2]]))
TREAT_YUCA<-factor(TREAT,levels=c("MOCK","DAOM197198","B1"))
Treatsss <- factor(paste(VAR,TREAT_YUCA,sep="."))
design <- model.matrix(~0+Treatsss)
colnames(design) <- c(levels(Treatsss))

# DATA FILTERING AND NORMALIZATION
DGE_YUCA <- DGEList(FOUR_VARS_YUCA)
#FILTERING  
keep_Y <- rowSums(DGE_YUCA$counts>100) >= 3
DGE_YUCA<- DGE_YUCA[keep_Y,]
#NORMALIZATION
DGE_YUCA <- calcNormFactors(DGE_YUCA)
DGE_YUCA_N <- voom(DGE_YUCA, design,plot=TRUE)

####################################################################################
fit_YUCA <- lmFit(DGE_YUCA_N,design)
cm <- makeContrasts(BRA337_B1.MOCK= BRA337.B1-BRA337.MOCK,
                    BRA337_DAOM197198.MOCK = BRA337.DAOM197198-BRA337.MOCK,
                    BRA337_DAOM197198.B1 = BRA337.DAOM197198-BRA337.B1,
                    CM4574_B1.MOCK= CM4574.B1-CM4574.MOCK,
                    CM4574_DAOM197198.MOCK = CM4574.DAOM197198-CM4574.MOCK,
                    CM4574_DAOM197198.B1 = CM4574.DAOM197198-CM4574.B1,
                    levels=design)

fit_YUCA<- eBayes(contrasts.fit(fit_YUCA, cm))

sign_YUCA_V1<-topTable(fit_YUCA,number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V1<-sign_YUCA_V1[sign_YUCA_V1$adj.P.Val<0.05,]





######################################################################################

#PREDICTED TARGET
target<-read.table('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v1_mirna/prediction_miRNA/psRNATargetJob-mesc305_vf.txt',h=T,sep = "\t") 
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

xtable(target[1:10,1:5])
head(target)

table(target$miRNA)

genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in%
                                  unique(target$Target), ] 
miRNA_origin<-target[target$Target %in% rownames(genes_YUCA_compe2),]

#colnames(genes_YUCA_compe2)<-NAME

sampleDists <- dist( t( DGE_YUCA_N$E ) )
sampleDistMatrix <- as.matrix( sampleDists )
pheatmap(sampleDistMatrix)

mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, VAR,TREAT_YUCA)
ggplot(mds, aes(X1,X2,shape=TREAT_YUCA)) + geom_point(size=6) +
  coord_fixed()


#genes_ggp<-cbind.data.frame(t(genes_YUCA_compe2),sapply(strsplit(colnames(genes_YUCA_compe2),"\\."), "[[", 1),sapply(strsplit(colnames(genes_YUCA_compe2),"\\."), "[[", 2))
#colnames(genes_ggp)[c(length(colnames(genes_ggp))-1,length(colnames(genes_ggp)))]<-c("VAR","TREAT")


#miRNA_origin$miRNA[miRNA_origin$Target %in% colnames(genes_ggp)[i]][1]
#miRNA_origin$At_ortholog[miRNA_origin$Target %in% colnames(genes_ggp)[i]][1]

mRNA_YUCA<-cbind.data.frame(t(genes_YUCA_compe2),VAR,TREAT_YUCA)
colnames(mRNA_YUCA)[c(length(colnames(mRNA_YUCA))-1,length(colnames(mRNA_YUCA)))]<-c("VAR","TREAT")
rownames(mRNA_YUCA)
#genes_ggp<-genes_ggp[grep("V",rownames(genes_ggp),invert=T),]

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/miRpredicted_target_mesculenta.pdf",width=10, height=10,useDingbats = F)

for (i in 1:(dim(mRNA_YUCA)[2]-2)) {
  p10 <- ggplot(mRNA_YUCA, aes_string(x = "VAR", y = colnames(mRNA_YUCA)[i], fill = "TREAT_YUCA")) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(paste(colnames(mRNA_YUCA)[i],miRNA_origin$miRNA[miRNA_origin$Target %in% colnames(mRNA_YUCA)[i]][1])) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) +
    scale_fill_brewer(palette = "Spectral")
  print(p10)
}
dev.off()


########################################################################################
########################################################################################
########################################################################################
### miRNA data
# do clustering

plant_samplesF2<-plant_samples

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.1",
              "CM4574-7.B1.2")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]
#reorder samples
colnames(plant_samplesF2)<-gsub("-","",colnames(plant_samplesF2))
colnames(plant_samplesF2) 

miRNA<-plant_samplesF2[ rowSums(plant_samplesF2) > 20, ]
miRNA<-t(miRNA)
#miRNA<-miRNA[grep( "miR171", rownames(miRNA)),]

samples.clusters <-sapply(strsplit(colnames(miRNA[1:5,]),"\\."), "[[", 2)

External_data<-plant_samplesF2[match(colnames(plant_samplesF2), rownames(genes_ggp))]


target[target$Target=="Manes.11G146500",]

########################################################################################################################################
########################################################################################################################################
## miRNA and mRNA on same plot
# Facet according to the cut variable


#same samples in data.frame
rownames(miRNA)<-gsub("47\\.","4\\.",rownames(miRNA))
colnames(miRNA)<-unlist(lapply(strsplit(colnames(miRNA)," "),function (x) x[[1]]))
colnames(miRNA)<-gsub("-","\\.",colnames(miRNA))

rownames(mRNA_YUCA)
colnames(mRNA_YUCA)


both_miRNA_mRNA<-merge(miRNA, mRNA_YUCA, by="row.names", all=TRUE)
rownames(both_miRNA_mRNA)<-both_miRNA_mRNA$Row.names
both_miRNA_mRNA$VAR<-unlist(lapply(strsplit(rownames(both_miRNA_mRNA),"\\."),function (x) x[[1]]))
both_miRNA_mRNA$TREAT_YUCA<-unlist(lapply(strsplit(rownames(both_miRNA_mRNA),"\\."),function (x) x[[2]]))

# filter by var
both_miRNA_mRNA<-both_miRNA_mRNA[grep("COL2215",rownames(both_miRNA_mRNA),invert=T),]
#both_miRNA_mRNA<-both_miRNA_mRNA[grep("B1",rownames(both_miRNA_mRNA),invert=T),]

# filter targets 
target_equiv<-target[,c(1,4)][!duplicated(target[,c(1,4)]),]
target_equiv<-target_equiv[target_equiv$miRNA_Acc. %in% all_sign_miRNA,]

#target_equiv<-target_equiv[target_equiv$miRNA_Acc. %in%  unique(c(sapply(strsplit(rownames(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,])," "), "[[", 1),
#                                                                  sapply(strsplit(rownames(CM4574_M_B1[CM4574_M_B1[,6]<0.1,])," "), "[[", 1)))
                          
                           
 #                          ,]

# different var filters
unique(c(sapply(strsplit(rownames(na.omit(BRA337_M_CAN[BRA337_M_CAN[,6]<0.1,]))," "), "[[", 1),
sapply(strsplit(rownames(na.omit(BRA337_M_B1[BRA337_M_B1[,6]<0.1,]))," "), "[[", 1)))

unique(c(sapply(strsplit(rownames(CM4574_M_CAN[CM4574_M_CAN[,6]<0.1,])," "), "[[", 1),
sapply(strsplit(rownames(CM4574_M_B1[CM4574_M_B1[,6]<0.1,])," "), "[[", 1)))

unique(c(sapply(strsplit(rownames(COL2215_M_CAN[COL2215_M_CAN[,6]<0.1,])," "), "[[", 1),
sapply(strsplit(rownames(COL2215_M_B1[COL2215_M_B1[,6]<0.1,])," "), "[[", 1)))


target_equiv[,1]<-gsub("-","\\.",target_equiv[,1])


colnames(both_miRNA_mRNA)

both_miRNA_mRNA

#install.packages("gridExtra")
library (gridExtra)


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_ALL.pdf",width=8, height=4,useDingbats = F)
cor_miRNA_mRNA<-list()
for (gene in 1:dim(target_equiv)[1]){
#plot
plot_miRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = target_equiv[gene,1], fill = "TREAT_YUCA")) +
  geom_boxplot(alpha=0.7)  + theme_bw()
plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
tryCatch(plot_mRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = as.character(target_equiv[gene,2]), fill = "TREAT_YUCA")) +
  geom_boxplot(alpha=0.7) + theme_bw() ,error=function(e) plot(1,1,main="NO mRNA sequenced"))
plot_mRNA = plot_mRNA + scale_fill_grey(start = 0, end = .9)

tryCatch(grid.arrange(plot_miRNA, plot_mRNA, nrow=1, ncol=2),error=function(e) 1)

#correlation
corre<-tryCatch(cbind(both_miRNA_mRNA[,grep(target_equiv[gene,2],colnames(both_miRNA_mRNA))],
             both_miRNA_mRNA[,grep( target_equiv[gene,1],colnames(both_miRNA_mRNA))]),error=function(e) NA)
corre<-tryCatch(corre[complete.cases(corre),],error=function(e) NA)
cor_miRNA_mRNA[[gene]]<-tryCatch(cor(corre[,1],corre[,2]),error=function(e) NA)
names(cor_miRNA_mRNA)[[gene]]<-paste(target_equiv[gene,1],as.character(target_equiv[gene,2]))

}
dev.off()

hist(unlist(cor_miRNA_mRNA),breaks=30,xlim = c(-1,1))
abline(v = -.8,col="red")
abline(v = .8,col="red")
names(cor_miRNA_mRNA[cor_miRNA_mRNA<(-0.6)])
cor_miRNA_mRNA[quantile(cor_miRNA_mRNA,.95)]

cor_miRNA_mRNA2<-unlist(cor_miRNA_mRNA)[!is.na(unlist(cor_miRNA_mRNA))]
quantile(cor_miRNA_mRNA2,c(.02,.98))

cor_miRNA_mRNA2[cor_miRNA_mRNA2<(-0.35)]

tosee<-na.omit(names(cor_miRNA_mRNA[cor_miRNA_mRNA<(-0.6)]))
unlist(lapply(strsplit(tosee," "),function (x) x[[1]]))
tosee
for i in 1:length(unlist(lapply(strsplit(tosee," "),function (x) x[[1]]))) {
plot_miRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = unlist(lapply(strsplit(tosee," "),function (x) x[[1]]))[i], fill = "TREAT_YUCA")) +
  geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
tryCatch(plot_mRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = unlist(lapply(strsplit(tosee," "),function (x) x[[2]]))[i], fill = "TREAT_YUCA")) +
           geom_boxplot(alpha=0.7) + theme_bw() ,error=function(e) plot(1,1,main="NO mRNA sequenced"))
plot_mRNA = plot_mRNA + scale_fill_grey(start = 0, end = .9)

tryCatch(grid.arrange(plot_miRNA, plot_mRNA, nrow=1, ncol=2),error=function(e) 1)
}

########################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#################################################################################################################################################
Ext
DW_total_plot<-Ext[,colnames(Ext) %in% c("DW_total","PLANT","TREAT")]

rownames(DW_total_plot)<-gsub("4-7.","4.",Ext$R_NAMES)

DW_total_plot<-DW_total_plot[match(rownames(both_miRNA_mRNA), rownames(DW_total_plot)),]

rownames(both_miRNA_mRNA)

both_miRNA_mRNA_DW<-cbind.data.frame(both_miRNA_mRNA,DW_total_plot$DW_total)

cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
      )
both_miRNA_mRNA_DW[,grep( target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))]

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_miRNA_DW.pdf",width=10, height=10,useDingbats = F)
cor_miRNA_DW<-list()
for (gene in 1:dim(target_equiv)[1]){
  #plot
  plot_miRNA<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = target_equiv[gene,1], fill = "TREAT_YUCA")) +
    geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
 plot_DW<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = "DW_total_plot$DW_total", fill = "TREAT_YUCA")) +
             geom_boxplot(alpha=0.7) + theme_bw() 
  plot_DW = plot_DW + scale_fill_grey(start = 0, end = .9)
  
grid.arrange(plot_miRNA, plot_DW, nrow=2, ncol=1)
  
  #correlation

  cor_miRNA_DW[[gene]]<-cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
  )
  
  names(cor_miRNA_DW)[[gene]]<-paste(target_equiv[gene,1],"DW")
  
}
dev.off()

hist(unlist(cor_miRNA_DW),breaks=30,xlim = c(-1,1))
abline(v = -.8,col="red")
abline(v = .8,col="red")
cor_miRNA_DW[cor_miRNA_DW<(-0.3)]
cor_miRNA_DW[quantile(cor_miRNA_DW,.95)]

cor_miRNA_DW2<-unlist(cor_miRNA_DW)[!is.na(unlist(cor_miRNA_DW))]
quantile(cor_miRNA_DW2,c(.02,.98))

cor_miRNA_DW2[cor_miRNA_DW2<(-0.35)]



#######################################################################################################################################


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_miRNA_DW.pdf",width=10, height=10,useDingbats = F)
cor_mRNA_DW<-list()
for (gene in 1:dim(target_equiv)[1]){
  #plot
  plot_miRNA<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = target_equiv[gene,2], fill = "TREAT_YUCA")) +
    geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
  plot_DW<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = "DW_total_plot$DW_total", fill = "TREAT_YUCA")) +
    geom_boxplot(alpha=0.7) + theme_bw() 
  plot_DW = plot_DW + scale_fill_grey(start = 0, end = .9)
  
  grid.arrange(plot_miRNA, plot_DW, nrow=2, ncol=1)
  
  #correlation
  
  cor_mRNA_DW[[gene]]<-cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,2],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
  )
  
  names(cor_mRNA_DW)[[gene]]<-paste(target_equiv[gene,2],"DW")
  
}
dev.off()

########################################################################################################################################












########################################################################################################################################
########################################################################################################################################
# AMF
########################################################################################################################################
########################################################################################################################################
# DATA SELECTION
########################################################################################################################################
# rename samples
summa<-read.table("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/summary_seq_results_v1.txt",h=T)
head(summa)
summa2<-summa[summa$Discard=="No",]
summa2<-na.omit(summa2)


# Selection all but co-inoculation
colnames(FOUR_VARS_AMF)
colnames(df_vf)[seq(3,90,2)]<-as.character(summa2$Name)

to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.2",
              "CM4574-7.B1.3")


FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("_3"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("_7"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("_11"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("_15"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("_19"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CANB1"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_16b"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V1"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8"), colnames(FOUR_VARS_AMF),invert =T)]

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
# DGE
## PLOT GENE EXPRESSION OF TARGET GENES

VAR<-c("CM4574","COL2215","COL2215","COL2215","COL2215","COL2215","COL2215","COL2215","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","BRA337","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574","CM4574")
TREAT<-c("DAOM197198","B1","MOCK","DAOM197198","B1","MOCK","MOCK","DAOM197198","DAOM197198","B1","MOCK","V5_14","V5_16","V5_17","B1","V5_20","MOCK","B1","MOCK","DAOM197198","B1","MOCK","V6_13","V6_14","MOCK","DAOM197198","B1","MOCK","DAOM197198")
REP<-c("2","1","1","1","2","3","5","3","2","1","1","V5_14","V5_16","V5_17","2","V5_20","4","3","5","4","1","1","V6_13","V6_14","2","3","4","4","4")
NAME<-c("CM4574.DAOM197198.2","COL2215.B1.1","COL2215.MOCK.1","COL2215.DAOM197198.1","COL2215.B1.2","COL2215.MOCK.3","COL2215.MOCK.5","COL2215.DAOM197198.3","BRA337.DAOM197198.2","BRA337.B1.1","BRA337.MOCK.1","BRA337.V5_14.V5_14","BRA337.V5_16.V5_16","BRA337.V5_17.V5_17","BRA337.B1.2","BRA337.V5_20.V5_20","BRA337.MOCK.4","BRA337.B1.3","BRA337.MOCK.5","BRA337.DAOM197198.4","CM4574.B1.1","CM4574.MOCK.1","CM4574.V6_13.V6_13","CM4574.V6_14.V6_14","CM4574.MOCK.2","CM4574.DAOM197198.3","CM4574.B1.4","CM4574.MOCK.4","CM4574.DAOM197198.4")

TREAT_AMF<-factor(TREAT,levels=c("MOCK","DAOM197198","B1"))
colnames(FOUR_VARS_AMF)<-NAME
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("COL2215.MOCK.5"), colnames(FOUR_VARS_AMF),invert =T)]

REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"\\."),function (x) x[[3]]))
VAR<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"\\."),function (x) x[[1]]))
TREAT<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"\\."),function (x) x[[2]]))
TREAT_AMF<-factor(TREAT,levels=c("MOCK","DAOM197198","B1"))
Treatsss <- factor(paste(VAR,TREAT_AMF,sep="."))
design <- model.matrix(~0+Treatsss)
colnames(design) <- c(levels(Treatsss))

# DATA FILTERING AND NORMALIZATION
DGE_AMF <- DGEList(FOUR_VARS_AMF)
#FILTERING  
keep_Y <- rowSums(DGE_AMF$counts>100) >= 3
DGE_AMF<- DGE_AMF[keep_Y,]
#NORMALIZATION
DGE_AMF <- calcNormFactors(DGE_AMF)
DGE_AMF_N <- voom(DGE_AMF, design,plot=TRUE)

#PREDICTED TARGET
target<-read.table('~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/psRNATargetJob-mesc305_vf.txt',h=T,sep = "\t") 
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

xtable(target[1:10,1:5])
head(target)

table(target$miRNA)

genes_AMF_compe2<-DGE_AMF$counts[rownames(DGE_AMF$counts) %in%
                                     unique(target$Target), ] 
miRNA_origin<-target[target$Target %in% rownames(genes_AMF_compe2),]

#colnames(genes_AMF_compe2)<-NAME

sampleDists <- dist( t( DGE_AMF_N$E ) )
sampleDistMatrix <- as.matrix( sampleDists )
pheatmap(sampleDistMatrix)

mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, VAR,TREAT_AMF)
ggplot(mds, aes(X1,X2,shape=TREAT_AMF)) + geom_point(size=6) +
  coord_fixed()


#genes_ggp<-cbind.data.frame(t(genes_AMF_compe2),sapply(strsplit(colnames(genes_AMF_compe2),"\\."), "[[", 1),sapply(strsplit(colnames(genes_AMF_compe2),"\\."), "[[", 2))
#colnames(genes_ggp)[c(length(colnames(genes_ggp))-1,length(colnames(genes_ggp)))]<-c("VAR","TREAT")


#miRNA_origin$miRNA[miRNA_origin$Target %in% colnames(genes_ggp)[i]][1]
#miRNA_origin$At_ortholog[miRNA_origin$Target %in% colnames(genes_ggp)[i]][1]

mRNA_AMF<-cbind.data.frame(t(genes_AMF_compe2),VAR,TREAT_AMF)
colnames(mRNA_AMF)[c(length(colnames(mRNA_AMF))-1,length(colnames(mRNA_AMF)))]<-c("VAR","TREAT")
rownames(mRNA_AMF)
#genes_ggp<-genes_ggp[grep("V",rownames(genes_ggp),invert=T),]

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/miRpredicted_target_mesculenta.pdf",width=10, height=10,useDingbats = F)

for (i in 1:(dim(mRNA_AMF)[2]-2)) {
  p10 <- ggplot(mRNA_AMF, aes_string(x = "VAR", y = colnames(mRNA_AMF)[i], fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7) +
    scale_y_continuous(name = "Normalized counts") +
    scale_x_discrete(name = "Plant variety") +
    ggtitle(paste(colnames(mRNA_AMF)[i],miRNA_origin$miRNA[miRNA_origin$Target %in% colnames(mRNA_AMF)[i]][1])) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16)) +
    scale_fill_brewer(palette = "Spectral")
  print(p10)
}
dev.off()


########################################################################################
########################################################################################
########################################################################################
### miRNA data
# do clustering

plant_samplesF2<-plant_samples

# filter out bad samples 
to_discard<-c("BRA337.MOCK.3","BRA337.CO_INOC.2",
              "BRA337.B1.4","BRA337.CO_INOC.5",
              "COL2215.MOCK.2","COL2215.MOCK.4",
              "CM4574-7.MOCK.3","CM4574-7.DAOM197198.1",
              "CM4574-7.B1.2")
plant_samplesF2<-plant_samplesF2[,!colnames(plant_samplesF2) %in% to_discard ]
plant_samplesF2<-plant_samplesF2[,grep("CO_INOC",colnames(plant_samplesF2),invert=T)  ]
#reorder samples
colnames(plant_samplesF2)<-gsub("-","",colnames(plant_samplesF2))
colnames(plant_samplesF2) 

miRNA<-plant_samplesF2[ rowSums(plant_samplesF2) > 20, ]
miRNA<-t(miRNA)
#miRNA<-miRNA[grep( "miR171", rownames(miRNA)),]

samples.clusters <-sapply(strsplit(colnames(miRNA[1:5,]),"\\."), "[[", 2)

External_data<-plant_samplesF2[match(colnames(plant_samplesF2), rownames(genes_ggp))]


target[target$Target=="Manes.11G146500",]

########################################################################################################################################
########################################################################################################################################
## miRNA and mRNA on same plot
# Facet according to the cut variable


target_equiv<-target[,c(1,4)][!duplicated(target[,c(1,4)]),]
target_equiv[,1]<-gsub("-","\\.",target_equiv[,1])
#same samples in data.frame
rownames(miRNA)<-gsub("47\\.","4\\.",rownames(miRNA))
colnames(miRNA)<-unlist(lapply(strsplit(colnames(miRNA)," "),function (x) x[[1]]))
colnames(miRNA)<-gsub("-","\\.",colnames(miRNA))

rownames(mRNA_AMF)
colnames(mRNA_AMF)


both_miRNA_mRNA<-merge(miRNA, mRNA_AMF, by="row.names", all=TRUE)
rownames(both_miRNA_mRNA)<-both_miRNA_mRNA$Row.names
both_miRNA_mRNA$VAR<-unlist(lapply(strsplit(rownames(both_miRNA_mRNA),"\\."),function (x) x[[1]]))
both_miRNA_mRNA$TREAT_AMF<-unlist(lapply(strsplit(rownames(both_miRNA_mRNA),"\\."),function (x) x[[2]]))

# filter by var
both_miRNA_mRNA<-both_miRNA_mRNA[both_miRNA_mRNA$VAR=="BRA337",]

colnames(both_miRNA_mRNA)

#install.packages("gridExtra")
library (gridExtra)


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_miRNA_mRNA.pdf",width=10, height=10,useDingbats = F)
cor_miRNA_mRNA<-list()
for (gene in 1:dim(target_equiv)[1]){
  #plot
  plot_miRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = target_equiv[gene,1], fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
  tryCatch(plot_mRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = as.character(target_equiv[gene,2]), fill = "TREAT_AMF")) +
             geom_boxplot(alpha=0.7) + theme_bw() ,error=function(e) plot(1,1,main="NO mRNA sequenced"))
  plot_mRNA = plot_mRNA + scale_fill_grey(start = 0, end = .9)
  
  tryCatch(grid.arrange(plot_miRNA, plot_mRNA, nrow=2, ncol=1),error=function(e) 1)
  
  #correlation
  corre<-tryCatch(cbind(both_miRNA_mRNA[,grep(target_equiv[gene,2],colnames(both_miRNA_mRNA))],
                        both_miRNA_mRNA[,grep( target_equiv[gene,1],colnames(both_miRNA_mRNA))]),error=function(e) NA)
  corre<-tryCatch(corre[complete.cases(corre),],error=function(e) NA)
  cor_miRNA_mRNA[[gene]]<-tryCatch(cor(corre[,1],corre[,2]),error=function(e) NA)
  names(cor_miRNA_mRNA)[[gene]]<-paste(target_equiv[gene,1],as.character(target_equiv[gene,2]))
  
}
dev.off()

hist(unlist(cor_miRNA_mRNA),breaks=30,xlim = c(-1,1))
abline(v = -.8,col="red")
abline(v = .8,col="red")
cor_miRNA_mRNA[cor_miRNA_mRNA<(-0.4)]
cor_miRNA_mRNA[quantile(cor_miRNA_mRNA,.95)]

cor_miRNA_mRNA2<-unlist(cor_miRNA_mRNA)[!is.na(unlist(cor_miRNA_mRNA))]
quantile(cor_miRNA_mRNA2,c(.02,.98))

cor_miRNA_mRNA2[cor_miRNA_mRNA2<(-0.35)]

plot_miRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = "osa.miR160f.5p", fill = "TREAT_AMF")) +
  geom_boxplot(alpha=0.7)  + theme_bw()
plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
tryCatch(plot_mRNA<-ggplot(both_miRNA_mRNA, aes_string(x = "VAR", y = "Manes.07G099400", fill = "TREAT_AMF")) +
           geom_boxplot(alpha=0.7) + theme_bw() ,error=function(e) plot(1,1,main="NO mRNA sequenced"))
plot_mRNA = plot_mRNA + scale_fill_grey(start = 0, end = .9)

tryCatch(grid.arrange(plot_miRNA, plot_mRNA, nrow=2, ncol=1),error=function(e) 1)


########################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#################################################################################################################################################
Ext
DW_total_plot<-Ext[,colnames(Ext) %in% c("DW_total","PLANT","TREAT")]

rownames(DW_total_plot)<-gsub("4-7.","4.",Ext$R_NAMES)

DW_total_plot<-DW_total_plot[match(rownames(both_miRNA_mRNA), rownames(DW_total_plot)),]

rownames(both_miRNA_mRNA)

both_miRNA_mRNA_DW<-cbind.data.frame(both_miRNA_mRNA,DW_total_plot$DW_total)

cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
)
both_miRNA_mRNA_DW[,grep( target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))]

pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_miRNA_DW.pdf",width=10, height=10,useDingbats = F)
cor_miRNA_DW<-list()
for (gene in 1:dim(target_equiv)[1]){
  #plot
  plot_miRNA<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = target_equiv[gene,1], fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
  plot_DW<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = "DW_total_plot$DW_total", fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7) + theme_bw() 
  plot_DW = plot_DW + scale_fill_grey(start = 0, end = .9)
  
  grid.arrange(plot_miRNA, plot_DW, nrow=2, ncol=1)
  
  #correlation
  
  cor_miRNA_DW[[gene]]<-cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,1],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
  )
  
  names(cor_miRNA_DW)[[gene]]<-paste(target_equiv[gene,1],"DW")
  
}
dev.off()

hist(unlist(cor_miRNA_DW),breaks=30,xlim = c(-1,1))
abline(v = -.8,col="red")
abline(v = .8,col="red")
cor_miRNA_DW[cor_miRNA_DW<(-0.3)]
cor_miRNA_DW[quantile(cor_miRNA_DW,.95)]

cor_miRNA_DW2<-unlist(cor_miRNA_DW)[!is.na(unlist(cor_miRNA_DW))]
quantile(cor_miRNA_DW2,c(.02,.98))

cor_miRNA_DW2[cor_miRNA_DW2<(-0.35)]



#######################################################################################################################################


pdf("~/Google Drive/Doctorat_shared_unil/miRNA-seq/Result_miR/v2_mirna/prediction_miRNA/biplot_miRNA_DW.pdf",width=10, height=10,useDingbats = F)
cor_mRNA_DW<-list()
for (gene in 1:dim(target_equiv)[1]){
  #plot
  plot_miRNA<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = target_equiv[gene,2], fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7)  + theme_bw()
  plot_miRNA = plot_miRNA + scale_fill_grey(start = 0, end = .9)
  plot_DW<-ggplot(both_miRNA_mRNA_DW, aes_string(x = "VAR", y = "DW_total_plot$DW_total", fill = "TREAT_AMF")) +
    geom_boxplot(alpha=0.7) + theme_bw() 
  plot_DW = plot_DW + scale_fill_grey(start = 0, end = .9)
  
  grid.arrange(plot_miRNA, plot_DW, nrow=2, ncol=1)
  
  #correlation
  
  cor_mRNA_DW[[gene]]<-cor(both_miRNA_mRNA_DW[,grep(target_equiv[gene,2],colnames(both_miRNA_mRNA_DW))],DW_total_plot$DW_total
  )
  
  names(cor_mRNA_DW)[[gene]]<-paste(target_equiv[gene,2],"DW")
  
}
dev.off()

########################################################################################################################################

