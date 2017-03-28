#Calculates network activity of the genes from the given expression data and generates a heatmap
#Swati Kaushik March 28 17

rm(list=ls())

library("gplots")
library("biomaRt")
setwd("/Users/swati/Desktop/Engelman-RNA-seq/")

args <- commandArgs(TRUE)

####
#use biomart
####

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)

restgenelist <- read.table(args[1], header=TRUE)
data <- getBM(attributes = c("ensembl_gene_id",'hgnc_symbol'), values=restgenelist$A, mart=mart, filters="ensembl_gene_id",uniqueRows = TRUE)
write.table(data, file="ensmblttohgnc.txt", sep="\t")
merged.data <- merge(data, restgenelist, by.x ="ensembl_gene_id" , by.y= "A")
write.table(merged.data, file ="merged.data", sep="\t")


#idx <- sapply(gene.names, function(x) {
#    which(subset.matrix$hgnc_symbol == x)
#})

####
##get names of the genes which are part of network activity
####

na.file = read.table(args[2], header=T) #FINAL-map
genes <- rownames(na.file)

# last name is the network acivity calculation so remove the last element of the vector	
gene.names <- genes[1:(length(genes)-1)]

####
## open the gene expression data 
####

merged.data <- read.table("merged.data" , header=T, sep="\t")

#to make exact match
add.pre <- paste("^", gene.names, sep="")
add.suff <- paste(add.pre,"$", sep="")

#grep the genes and rows
subset.matrix <- merged.data[grep(paste(add.suff , collapse="|"), merged.data$hgnc_symbol), ]
#sorted.subset.matrix <- subset.matrix[idx,]
sorted.subset.matrix.c <- subset.matrix[match(gene.names, subset.matrix$hgnc_symbol),]
sorted.subset.matrix <- sorted.subset.matrix.c[,order(colnames(sorted.subset.matrix.c))]
write.table(sorted.subset.matrix, file = "sorted.subset.matrix.txt", sep="\t")

####
## Z score calculation on the subsetted matrix
####

#replace the rownames with actual hgnc symbols 
rownames(sorted.subset.matrix) <- sorted.subset.matrix$hgnc_symbol
#take log of CPM values
sorted.subset.matrix.log <- log(sorted.subset.matrix[,3:ncol(sorted.subset.matrix)])
write.table(sorted.subset.matrix.log, file= "sorted.subset.matrix.log", sep="\t")

#calculate Zscore
sorted.subset.matrix.log.T <- scale(t(sorted.subset.matrix.log))
attr(sorted.subset.matrix.log.T, "scaled:scale") =NULL
attr(sorted.subset.matrix.log.T, "scaled:center") =NULL
sorted.subset.matrix.log.T.t <- t(sorted.subset.matrix.log.T)
write.table(sorted.subset.matrix.log.T.t, file ="sorted.subset.matrix.zscore", sep="\t")

####
# calculate network activity
####

output <- vector()
#take colnames in a vector
p <- as.vector(colnames(sorted.subset.matrix.log.T.t))

for (i in 1:ncol(sorted.subset.matrix.log.T.t)){

	sum.cols <- sum(sorted.subset.matrix.log.T.t[1:14,i])
	sub.cols <- sum(sorted.subset.matrix.log.T.t[15:nrow(sorted.subset.matrix.log.T.t),i])
	network.activity <- sum.cols - sub.cols
	no.data <- paste(p[i], network.activity,collapse='	')
	output <- rbind(output, data.frame(p[i],network.activity))
}

rownames(output) <- output$p.i.
write.table(output[,2], file = "NA.out", sep="\t", row.names = rownames(output))

####
# calculate z score of NA
####

na.vector <- t(output[,2])
na.vectpr.t <- scale(t(na.vector))
attr(na.vectpr.t, "scaled:scale") =NULL
attr(na.vectpr.t, "scaled:center") =NULL
na.vector <- t(na.vectpr.t)

####
#create all NA component matrix
####

new.matrix <- rbind(sorted.subset.matrix.log.T.t, na.vector)
names.matrix <- rownames(sorted.subset.matrix)
appended.names.matrix <- append(names.matrix, "EGFR-NA")
rownames(new.matrix) <- appended.names.matrix
write.table(new.matrix, file = "New.entire.matrix.out", sep="\t")

####
#Take average of replicates
####

new.matrix.t <- t(new.matrix)
names.row <- rownames(new.matrix.t)
rownames(new.matrix.t) <- gsub("\\..*","",names.row)
averaged.new.matrix <- aggregate(new.matrix.t[, 1:ncol(new.matrix.t)], list(rownames(new.matrix.t)), mean)
averaged.new.matrix.t <- t(averaged.new.matrix)
colnames(averaged.new.matrix.t) <- averaged.new.matrix.t[1, ]
averaged.new.matrix.t <- averaged.new.matrix.t[-1, ] 
#averaged.new.matrix.t.nq <- noquote(averaged.new.matrix.t) 
write.table(averaged.new.matrix.t, file = "averaged.new.matrix.txt" ,sep="\t")

#averaged.new.matrix.t.nq [averaged.new.matrix.t.nq  == " "] <- ""

####
#plot heatmap
####
#file = read.table("averaged.new.matrix.txt", header=T)
pdf(file="heatmap-averaged.pdf")
#col.corr <- colorRampPalette(c("white","red"),(10))
heatmap.2 (as.matrix(file),
           Rowv = TRUE,
           Colv = TRUE,
           distfun = dist,
           hclustfun = hclust,
           dendrogram = "both",
           scale = "none",
           na.rm = TRUE,
           col = bluered(256),
           trace = "none",
           cexRow = 0.6, 
		   cexCol = 0.7)
		   
dev.off()	

####

