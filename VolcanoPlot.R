library("ggplot2")
library("ggrepel")
library("dplyr")

# volcano plot function used to make custom volcano plots and label genes of interest.

volcano.plot<-function(rnaseqdata, alphavalue = 0.05, lfcLim = 1.5, top = 10, genesetFile = NULL, graphname = "volcano",
                       xlims = c(-2.5,2.5), ylims = c(0,20), width = 6, height = 6, pointsize = 1.5,
                       boxpadding = 0.5, labelRepel = 10, labelSize = 5, axisLabelsize = 12, axisTitle = 14) {
  datadf <- read.table(rnaseqdata, sep = "\t", header = T)
  datadf$GeneSymbol <- as.character(datadf$GeneSymbol)
  datadf$padj <- p.adjust(datadf$pvalue, method = "BH")
  datadf <- na.omit(datadf)
  
  # manages gene labels
  input <- mutate(datadf, sig = ifelse(datadf$padj < alphavalue & abs(datadf$log2FoldChange) >= lfcLim, "Sig", "Not_Sig")) #Will have different colors depending on significance
  input$label <- ""
  if (is.null(genesetFile) == F) {
    geneset <- read.table(genesetFile, stringsAsFactors = F)
    geneset$rownum <- match(geneset[,1], datadf$GeneSymbol)
    input[geneset[,2],length(input)] <- as.character(geneset[,1])
  } else if (top > 0 && is.null(genesetFile) == T) {
    input$label[1:top] <- input$GeneSymbol[1:top]
  }
  # changes outliers to triangles
  input$shape <- ifelse(-log10(input$padj) > ylims[2] | abs(input$log2FoldChange) > xlims[2], "triangle", "circle")
  input$padj[-log10(input$padj) > ylims[2]] <- 10^-(ylims[2])
  input$log2FoldChange[input$log2FoldChange > xlims[2]] <- xlims[2]
  input$log2FoldChange[input$log2FoldChange < -xlims[2]] <- -xlims[2]
  # begins plotting
  v <- ggplot(input, aes(log2FoldChange, -log10(padj))) #volcanoplot with log2Foldchange versus pvalue
  v + geom_point(aes(col = sig, shape=shape), size = pointsize) + #add points colored by significance
    scale_color_manual(values = c("grey50","red")) + xlim(xlims) + ylim(ylims) +
    xlab("log2(fold change)") + ylab("-log10(FDR)") + theme_bw() + 
    theme(legend.position = "none", axis.text = element_text(size=axisLabelsize), axis.title = element_text(size=axisTitle,face="bold")) + #ggtitle("Volcanoplot")
    geom_text_repel(label = input$label, box.padding = boxpadding, size = labelSize, force = labelRepel, segment.size = 0.3, 
                    min.segment.length = 0.1, segment.alpha = 0.5, na.rm = TRUE, ylim = c(ylims[1],ylims[2])) #adding text for the top 20 genes
  ggsave(paste(graphname, "_Volcanoplot.pdf", sep = ""), width = width, height = height)
}
