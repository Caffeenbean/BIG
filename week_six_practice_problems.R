
# 1. Find a data set on this page: http://www.bioconductor.org/packages/release/data/experiment/


# 2. Install and load the data set
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("curatedAdipoRNA")

# 3. Convert the data to a DESeq dataset object with an appropriate design formula

library(curatedAdipoRNA)

suppressMessages( library("DESeq2") ) #Need these libraries
suppressMessages( library("ggplot2") )
data("adipo_counts") #load data 
se = adipo_counts    #place all rna seq data from NKI into variable se

pThr        <- 0.01   # for adj. p-value threshold (<= 0.01)
logFCThr    <- 1      # for log2 fold-change (at least 2 fold)
baseMeanThr <- 20     # for average expression level (at least 20 counts).
cpmThr      <- 1      # for copy-per-million (at least 1 cpm). 

head( assay( se ) )  # examine all data component
rowData( se )
colData( se )

dds <- DESeqDataSet(se, design = ~ stage)#how h3k4me3 changes in different stages of cancer


# 4. Run DESeq

dds <- DESeq( dds )

# 5. Get the results 

res <- results( dds )


head(res)

# 6. If not already done, convert the gene symbols to gene IDs
BiocManager::install("org.Mm.eg.db") #Mouse Gene Annotations
library(org.Mm.eg.db) 
library( AnnotationDbi ) 
library( dplyr ) 

keytypes(org.Mm.eg.db)

rownames(res)


anno <- anno <- AnnotationDbi::select(org.Mm.eg.db, rownames( res ), 
                                      columns=c("SYMBOL", "ENTREZID", "SYMBOL", "GENENAME"), 
                                      keytype="SYMBOL")

# adding ENSEMBL gene ID as a column in significant differentially expression gene table.
res = cbind( SYMBOL = rownames( res), res )
wholeTable <- left_join( as.data.frame( res ), anno )

# let's take a look...
head( wholeTable ) 




# 7. Use EnhancedVolcano to plot the results

library( EnhancedVolcano )
EnhancedVolcano( as.data.frame(wholeTable), lab = wholeTable$SYMBOL, 
                 x = 'log2FoldChange', y = 'padj',
                 xlim = c(-8, 8), title = 'H3K4me3 change over time',
                 pCutoff = 0.01, FCcutoff = 2, pointSize = 2.0, 
                 labSize = 3.0 )


# 8. Save the results to a .csv file 

# When finished, send your R script to Henry Miller (millerh1@uthscsa.edu)
# Make sure to explain what comparison you made with DESeq and why.

# For example -- to test the effect of dex on gene expression in the "airway" dataset, I would start with...
library(airway)
dds <- DESeqDataSet(airway, design = ~dex)
# And then continue the analysis from there...




