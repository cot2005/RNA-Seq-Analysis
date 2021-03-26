# RNA-Seq-Analysis
This bulk RNA-seq analysis uses Salmon aligned data to perform differential expression analysis (Using DESeq2), make volcano plots, and output GSEA prerank files(.rnk). The pipeline requires a conditions file which specifies the sample names (they should match the names of the .sf files), which group they belong to and the analysis (output file name) they will be a part of. Different output file names will perform separate differential expression analyses. The analysis is performed as treated/untreated. The genesetFile input allows for customization of the labeling of the volcano plots. The genesetFile must be a list of genes to be labeled in the volcano plots.

The conditions table is a tab-delimited file with the format below:
```
Sample	Condition	OutputFileName
Sample1_replicate1	untreated	Sample1vsSample2
Sample1_replicate2	untreated	Sample1vsSample2
Sample2_replicate1	treated Sample1vsSample2
Sample2_replicate2	treated	Sample1vsSample2
Sample1_replicate1	untreated	Sample1vsSample3
Sample1_replicate2	untreated	Sample1vsSample3
Sample3_replicate1	treated	Sample1vsSample3
Sample3_replicate2	treated	Sample1vsSample3
(etc)
```

The pipeline can be run on pre-analyzed count data by inputing the name of the count file as the RNAcountFile. Pre-analyzed count data must have the following format:
```
GeneID	GeneSymbol	Sample1	Sample2	Sample3 (etc)
ENSG_ID1	GENE1	1000	1100	1200
ENSG_ID2	GENE2	1 2 3
(etc)
```

# Volcano Plot
This plotting function is used to create volcano plots from the output log fold change file from the RNA-seq pipeline. The plots and labels can be customized through the inputs of this function. Plot will automatically color significant points as red and convert genes outside of the graph limits to triangles.
```
Usage:
rnaseqdata = file name containing the differential expression data
alphavalue = threshold for FDR significance
lfcLim = threshold for log fold change significance
top = 10
genesetFile = optional input for gene labels. Input file name containing the genes to be labeled.
graphname = graph name prefix
xlims = input to set x-axis limits, defaults to +/- 2.5
ylims = input to set y-axis limits, defaults to 0,20
width = width for output plot size, defaults to 6
height = height for output plot size, defaults to 6
boxpadding = box padding for labels by ggrepel, defaults to 0.5
labelRepel = force for labels by ggrepel, defaults to 10
labelSize = font size for labels, defaults to 5
axisLabelsize = axis font size, defaults to 16
axisTitle = axis title font size, defaults to 18
```
