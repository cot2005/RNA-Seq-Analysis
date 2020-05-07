# RNA-Seq-Analysis
This bulk RNA-seq analysis uses Salmon aligned data to perform differential expression (Using DESeq2), make volcano plots, and output GSEA prerank files(.rnk). The pipeline requires a conditions file which specifies the sample names (they should match the names of the .sf files), which group they belong to and the analysis (output file name) they will be a part of. Different output file names will perform separate differential expression analyses. The analysis is performed as treated/untreated.

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
