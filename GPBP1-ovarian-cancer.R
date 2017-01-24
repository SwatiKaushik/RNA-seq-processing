#!/usr/bin/env Rscript
# Take set of genelist and compare the expression of genes with CNV of another selected gene (eg. here - GPBP1)
# Compare the expression in gene loss vs gene diploid samples and identifies the genes with significant difference in expression


setwd("/Users/swati/Desktop")

rm(list=ls())

substrRight <- function(x, n){
    substr(x, 1, n)
}

gene.exp <- read.table("unc.edu_OV_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor" , header=T)
m.gene.exp=sapply(strsplit(rownames(gene.exp), "[|]"), `[`, 1)
rownames(gene.exp) <- m.gene.exp

genename = read.table("genename.txt", header=T)
#genename =c("A1CF", "A4GNT", "A4GALT")
subset.gene.exp <- gene.exp[grep(paste(genename$genes, collapse="|"), rownames(gene.exp)),]
subset.gene.exp.log <- log(subset.gene.exp)
write.table(subset.gene.exp.log, file= "genenames.exp.log", sep="\t")

subset.gene.exp.log  <- read.table("genenames.exp.log", header=T)
subset.gene.exp.log[subset.gene.exp.log == -Inf] <- 0
t.subset.gene.exp.log = scale(t(subset.gene.exp.log))
attr(t.subset.gene.exp.log, "scaled:scale") =NULL
attr(t.subset.gene.exp.log, "scaled:center") =NULL
r.t.subset.gene.exp.log = round(t.subset.gene.exp.log, digits = 3)
t.gene.exp = t(r.t.subset.gene.exp.log)
reduced.col.t.gene.exp <- substrRight(colnames(t.gene.exp), 15)
colnames(t.gene.exp) <- reduced.col.t.gene.exp
write.table(t.gene.exp, file= "genename.log.Zscored.txt", sep="\t")

cnv.file <- read.table("OV.gistic.all_thresholded.by_genes", header=T, sep="\t")
gpbp1.cnv.data <- cnv.file[grep("^GPBP1$",cnv.file$Gene.Symbol),]
reduced.col <- substrRight(colnames(gpbp1.cnv.data), 15)
colnames(gpbp1.cnv.data) <- reduced.col
write.table(gpbp1.cnv.data, file="GPBP1.cnv.ovarian.txt", sep="\t",row.names = FALSE)

df1 <- t.gene.exp[,intersect(colnames(gpbp1.cnv.data), colnames(t.gene.exp))]
write.table(df1, file="intersect.GPBP1.cnv.txt", sep="\t")
df2 <- gpbp1.cnv.data[,intersect(colnames(gpbp1.cnv.data), colnames(t.gene.exp))]
write.table(df2, file="intersect.gene.exp.txt", sep="\t")

GPBP1.cnv.interesting.gene.exp <- rbind(df2,df1)
write.table(GPBP1.cnv.interesting.gene.exp, file="GPBP1.cnv.interesting.gene.exp.txt", sep="\t")
trans.GPBP1.cnv.interesting.gene.exp <- t(GPBP1.cnv.interesting.gene.exp)
write.table(trans.GPBP1.cnv.interesting.gene.exp, file="GPBP1.cnv.interesting.gene.exp.T", sep="\t")

trans.file = read.table("GPBP1.cnv.interesting.gene.exp.T", header=T)
sorted.matrix <- trans.file[order(trans.file$X6359),]
#sorted.matrix <- t(apply(t(GPBP1.cnv.interesting.gene.exp),1,sort))
#t.sorted.matrix =t(sorted.matrix)
write.table(sorted.matrix, file="GPBP1.cnv.interesting.gene.exp.sorted.txt", sep="\t")
ttest <- t(apply(sorted.matrix[-1],2,function(x) unlist({t.test(x[1:17],x[150:247])[c("p.value","estimate")]})))
write.table(ttest, file="-2-0.txt", sep="\t")
ttest2 <- t(apply(sorted.matrix[-1],2,function(x) unlist({t.test(x[1:149],x[150:247])[c("p.value","estimate")]})))
write.table(ttest2, file="-2-1-0.txt", sep="\t")

#t.sorted.matrix <- read.table("GPBP1.cnv.interesting.gene.exp.sorted.txt", header=T, sep="\t")

#file = read.table("boxplot-input.txt", header=T, fill= TRUE, sep="\t")
#boxplot(file$TOP3A[1:149],file$TOP3A[150:247],file$RPA1[1:149], file$RPA1[150:247],file$BRCA1[1:149], file$BRCA1[150:247], file$XRCC3[1:149], file$XRCC3[150:247], file$XRCC2[1:149], file$XRCC2[150:247], file$RAD51AP1[1:149], file$RAD51AP1[150:247], file$TOP3B[1:149], file$TOP3B[150:247], file$BRCA2[1:149], file$BRCA2[150:247], notch = TRUE, col=c("mediumaquamarine", "slategray3"), lwd=1.5,at =c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17, 19,20, 22,23),xaxt='n')
#at1 <- seq(2, 24, 3)
# axis(side =1, at1, labels = F)
# legend("topright", legend=c("GPBP1 loss", "GPBP1 diploid"),fill=c("mediumaquamarine", "slategray3"),box.lwd = 0,box.col = "white",bg = "white")
