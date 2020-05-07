#author: Colin Tang
#
#RNA-seq pipeline for analyzing salmon aligned quant.sf files.
#Will perform normalizations and differential expression analysis.
#outputs abundance (TPM) file, count file, LFC files indicated in the conditions input,
#write GSEA files (.rnk), make volcano plots for each LFC file, make PCA graph for all input data
#
#run in the folder with all the .sf files and conditions file input with the format: (Sample, Condition, OutputFileName)
#sample names will be the same as the filenames without the .sf extension.
#LFC analyses are separated by unique outputfilenames.

library("ggbiplot")
library("rgl")
library("tximport")
library("biomaRt")
library("ggplot2")
library("ggrepel")
library("dplyr")

RNAseq.Pipeline<-function(conditionsfile) {
  #processes quant.sf files to counts
  RNApipeline.salmon.process()
  rnaData <- read.table("salmongene_counts.txt", sep = "\t", header = T)
  
  #Checks for duplicated genes in the counts file and deletes the second entry if found.
  duplicatedGenes <- which(duplicated(rnaData[,1]) == T)
  if (length(duplicatedGenes) > 0) {
    print("Duplicate ENSGs. Deleted entries:")
    print(rnaData[duplicatedGenes,1:2])
    rnaData <- rnaData[-duplicatedGenes,]
  }
  rownames(rnaData) <- rnaData[,1]   #makes rownames the unique ensembl gene ID
  geneKey <- rnaData[,1:2]
  
  #reorganizes data, makes pca groups table, and makes PCA graph
  PCAdata <- read.table("salmongene_abundance.txt", sep = "\t", header = T)
  PCAdata <- PCAdata[,-1:-2]
  RNAsamples <- colnames(PCAdata)
  labelsdf <- t(as.data.frame(strsplit(RNAsamples, "_")))
  PCAsamples <- paste(labelsdf[,3], labelsdf[,4], labelsdf[,5], sep = "_")
  colnames(PCAdata) <- PCAsamples
  PCAconditions <- data.frame(Samples = PCAsamples, 
                              Condition = factor(paste(labelsdf[,2], labelsdf[,3], sep = "_")))
  RNApipeline.pca(PCAdata, PCAconditions)
  
  conditionsdf <- read.table(conditionsfile, sep = "\t", header = T)
  #determines number of analyses for deseq loop
  listAnalysis <- levels(conditionsdf$OutputFileName)
  
  #performs all the differential analysis specified in the conditions file
  for (i in listAnalysis) {
    #gets RNAdata columns for samples in the ith analysis
    treatedSamples <- conditionsdf$Sample[which(conditionsdf$OutputFileName == i & conditionsdf$Condition == "treated")]
    treatedCols <- match(treatedSamples, colnames(rnaData))
    untreatedSamples <- conditionsdf$Sample[which(conditionsdf$OutputFileName == i & conditionsdf$Condition == "untreated")]
    untreatedCols <- match(untreatedSamples, colnames(rnaData))
    #creates a RNA data subset df for each analysis
    tempRNAdata <- as.matrix(subset(rnaData, select = c(treatedCols, untreatedCols)))
    conditions <- conditionsdf[which(conditionsdf$OutputFileName == i),1:2]
    conditions <- conditions[match(colnames(tempRNAdata), conditions$Sample),]
    colnames(conditions) <- c("sample", "condition")
    
    #performs differential expression
    RNAlfc <- RNApipeline.DESeq2(tempRNAdata, conditions, geneKey)
    #writes LFC and gsea files
    write.table(RNAlfc, paste(i, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(subset(RNAlfc, select = c(2,4)), paste(i, "_GSEA.rnk", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F)
    #makes volcano plot
    RNApipeline.volcano.plot(RNAlfc, alphavalue = 0.05, topgenes = 20, graphname = i)
  }
}


#processes .sf files to counts and abundance files
#
RNApipeline.salmon.process<-function() {
  dir <- system.file("extdata", package = "tximportData")
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
  
  quantFiles <- list.files(pattern = ".sf")
  txi <- tximport(quantFiles, type = "salmon", tx2gene = tx2gene)
  
  sampleNames <- gsub(".sf", "", quantFiles)
  abundancedf <- data.frame(rownames(txi$abundance), txi$abundance)
  countsdf <- data.frame(rownames(txi$counts), txi$counts)
  colnames(abundancedf) <- c("gene_id", sampleNames)
  colnames(countsdf) <- c("gene_id", sampleNames)
  
  abundancedf$ensembl_id <- gsub('\\..+$', '', abundancedf$gene_id)
  countsdf$ensembl_id <- gsub('\\..+$', '', countsdf$gene_id)
  #gets symbols from biomart
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  symbol <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), filters = "ensembl_gene_id", values = abundancedf$ensembl_id, mart = mart)
  
  abundancedf <- merge(x = symbol, y = abundancedf, by.x="ensembl_gene_id", by.y="ensembl_id")
  countsdf <- merge(x = symbol, y = countsdf, by.x="ensembl_gene_id", by.y="ensembl_id")
  
  #removes empty gene names, and ensembl id column
  abundancedf <- abundancedf[-which(abundancedf[,2] == ""),]
  countsdf <- countsdf[-which(countsdf[,2] == ""),]
  abundancedf$ensembl_gene_id <- NULL
  countsdf$ensembl_gene_id <- NULL
  abundancedf <- subset(abundancedf, select = c(2,1,3:length(abundancedf)))
  countsdf <- subset(countsdf, select = c(2,1,3:length(countsdf)))
  countsdf[,3:length(countsdf)] <- round(countsdf[,3:length(countsdf)])   #rounds counts to nearest integer
  
  colnames(abundancedf)[1:2] <- c("GeneID", "GeneSymbol")
  colnames(countsdf)[1:2] <- c("GeneID", "GeneSymbol")
  write.table(abundancedf, "salmongene_abundance.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(countsdf, "salmongene_counts.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}


#input selected RNA seq data matrix and conditions dataframe with no extraneous data
#will run differential expression using the conditions dataframe and will return the LFC dataframe
#
RNApipeline.DESeq2<-function(tempRNAdata, conditions, geneKey) {
  dds <- DESeqDataSetFromMatrix(countData = tempRNAdata, colData = conditions, design = ~ condition)
  
  #prefiltering removes low count features of < 10 reads for a row of a feature
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #factor leveling
  dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
  #differential expression calculation
  dds <- DESeq(dds)
  res <- results(dds)
  res <- results(dds, contrast=c("condition","treated","untreated"))
  summary(res)
  resOrdered <- res[order(res$pvalue),]
  resOrdered <- data.frame(GeneID = rownames(resOrdered), GeneSymbol = geneKey[match(rownames(resOrdered), geneKey[,1]),2], resOrdered)
  return(resOrdered)
}


#makes volcanoplot
RNApipeline.volcano.plot<-function(datadf, alphavalue = 0.05, topgenes = 20, graphname = "Volcanoplot") {
  datadf <- na.omit(datadf)
  input <- mutate(datadf, sig = ifelse(datadf$padj < alphavalue, paste("FDR<",alphavalue, sep = ""), "Not Sig")) #Will have different colors depending on significance
  v <- ggplot(input, aes(log2FoldChange, -log10(pvalue))) #volcanoplot with log2Foldchange versus pvalue
  v + geom_point(aes(col = sig)) + #add points colored by significance
    scale_color_manual(values = c("red","black")) + 
    ggtitle(paste(graphname, "Volcano Plot")) + #e.g. 'Volcanoplot DESeq2'
    geom_text_repel(data = head(input, 20), aes(label = GeneSymbol)) #adding text for the top 20 genes
  ggsave(paste(graphname, "_volcano.pdf", sep = ""), width = 8,height = 6) #In case you want to easily save to disk
}

#expects data with no gene names or IDs as columns
RNApipeline.pca<-function(pcaData, pcaGroups) {
  pcaData <- na.omit(pcaData)
  screenmatrix <- log(pcaData)
  screenmatrix <- screenmatrix[is.finite(rowSums(screenmatrix)),]
  screenmatrix <- t(screenmatrix)
  
  pca.input <- screenmatrix[,apply(screenmatrix,2,var,na.rm=TRUE) != 0]
  pca.output <- prcomp(pca.input, center = T, scale. = T) # computes variance
  ggbiplot(pca.output, cex.lab=3, scale = .75, choices = c(1,2), ellipse = F,
           groups = pcaGroups[,2], labels.size = 2,
           var.axes = 'F', label = pcaGroups[,1]) + 
    geom_point(aes(color = pcaGroups[,2]), size=7, pch=19, alpha=0.5) +
    theme(title = element_text(size = 20),axis.title = element_text(size=20),axis.text = element_text(size=20)) +
    ggtitle('Unsupervised PCA', subtitle = "Using Genes from Panel")
  ggsave("PCAnalysis.pdf", width = 12, height = 8)
}


