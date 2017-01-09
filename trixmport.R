# process data from kallisto runs and calculate differentially expressed genes

# Deseq2 can't process transcript level data. So I have used R package called tximport to get gene level estimated counts

# tximport needs a transcript conversion file
#run biomart to convert transcript IDs to hgnc (gene ids)

library(biomaRt)
listEnsembl()
ensembl=useMart("ensembl")
listDatasets(ensembl)
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(mart)
attributes= listAttributes(ensembl)
file = read.table("/Users/swati/Desktop/transcripts", header=F)
A=getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), values=file$V1, mart=mart, filters='ensembl_transcript_id',uniqueRows = TRUE)
write.table(A, file="/Users/swati/Desktop/Transcriptogene.txt", sep="\t")

# use tximport
library(tximportData)
setwd("/Volumes/sblab/users/swatik/Kevin-RNA-seq-data/3kallisto/Bootstrapped-kallisto/")
dir <- "/Volumes/sblab/users/swatik/Kevin-RNA-seq-data/3kallisto/Bootstrapped-kallisto"
sample <- read.table("samples", header = T)
files <- file.path(dir, sample$samplesrun, "abundance.tsv")
txt2gene <- read.table("Txt2gene.txt", heade=T)
names(files) <- paste0("sample", 151:158)
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = txt2gene)
write.table(txi.kallisto, file="txi.kallisto.txt", sep="\t")


# get prepared for deseq2
library(DESeq2)
sampleTable <- data.frame(condition = factor(rep(c("NT-DMSO", "NT-BMN","GPBP1-DMSO", "GPBP1-BMN"), each = 2)))
rownames(sampleTable) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)


# run deseq2 
dds <- DESeq(dds)
res <- results(dds)

##################  Calculate fold change per condition #############################################
#fold change GPBP1-DMSO vs NT-DMSO 

NT.DMSO_GPBP1.DMSO <- results (dds,contrast=c("condition","GPBP1-DMSO", "NT-DMSO"))
write.table(NT.DMSO_GPBP1.DMSO, file="GPBP1.DMSO.NT.DMSO", sep="\t") 
 
#GPBP1.BMN_NT.BMN

GPBP1.BMN_NT.BMN <- results (dds,contrast=c("condition","GPBP1-BMN", "NT-BMN")) 
write.table(GPBP1.BMN_NT.BMN, file="GPBP1.BMN_NT.BMN", sep="\t") 

#GPBP1.BMN_GPBP1.DMSO 
GPBP1.BMN_GPBP1.DMSO <- results (dds,contrast=c("condition","GPBP1-BMN", "GPBP1-DMSO"))
write.table(GPBP1.BMN_GPBP1.DMSO, file="GPBP1.BMN_GPBP1.DMSO", sep="\t") 

#GPBP1.BMN_GPBP1.DMSO 
GPBP1.BMN_NT.BMN <- results (dds,contrast=c("condition","GPBP1-BMN", "NT-BMN"))
write.table(GPBP1.BMN_NT.BMN , file="GPBP1.BMN_NT.BMN", sep="\t") 
 
##################################################################################################### 
# base mean for all the samples
#The base mean is the mean of normalized counts of all samples, normalizing
#for sequencing depth.

#It does not take into account gene length. The base mean is used in DESeq2
#only for estimating the dispersion of a gene (it is used to estimate the
#fitted dispersion). For this task, the range of counts for a gene is
#relevant but not the gene's length (or other technical factors influencing
#the count, like sequence content).

rowMeans(counts(dds, normalized=TRUE))

# to get normalized count
dds <- estimateSizeFactors(dds)
counts(dds, normalized=TRUE)
