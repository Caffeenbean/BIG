setwd("C:/Users/s_aio/Desktop/BIG/Week 8")

#### Week 8: Biological interpretation ####

# Download any missing packages:
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("recount")
BiocManager::install("apeglm")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("clusterProfiler")
install.packages("msigdbr")

# Load libraries
library(recount)
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(clusterProfiler)

## Effects of EWS-FLI1 KD in Ewing sarcoma ##

# Biological question: What is the effect of knocking
# down EWS-FLI1 in Ewing sarcoma? How does it change
# cellular behavior?

# The dataset we use here is downloaded from recount2
# Find a project of interest (note: it's hit or miss, the database is not limitless)
project_info <- abstract_search(query = "Ewing sarcoma")

# View the available projects and select one
View(project_info)
selected_study <- "SRP015989"

# Download the study data
download_study(selected_study)

# Load the data
load("SRP015989/rse_gene.Rdata")

class(rse_gene)

# Examine the read count matrix
cts <- assay(rse_gene)
View(cts)

# Examine the sample metadata to see if it's ready for DESeq2 analysis
cd <- as.data.frame(colData(rse_gene))
View(cd) #we don't know which are scrambled and which are controls: 
#go to sra and look up experiment type: 

# Fix the colData to give a column with the appropriate groups
rse_gene$condition <- c(rep("shCTR", 3), rep("shEF1", 4)) #using repeat command

#so now you have which ones are controls and which are knockdown

# Do DESeq2 analysis
# -- Make the dds
dds <- DESeqDataSet(rse_gene, design = ~condition)#we want to compare shCTR to shEF1
# -- Get rlog
rld <- rlog(dds)#tranform rlog and make sure that your data variants look ok (they are nicely split between biological groups)
# -- plot PCA
plotPCA(rld) #this means that the technical aspect such as the way you ran your experiment isn't nearly as important as your biological component
# -- analyze
dds <- DESeq(dds)
# -- get results
res <- results(dds, contrast = c("condition", "shEF1", "shCTR")) #contrast will you give you the fold difference
# -- plotMA
plotMA(res)#how does your fold change relate to the degree that each gene was expressed
# LFC shrink-a more accurate way to calculate fold change-where it finds where the fold change was exaggerated and shirnks them to a more appropriate place
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)
# -- plotMA with resNorm
plotMA(resNorm)
# Make a DF
resdf <- as.data.frame(resNorm)
View(resdf)
# -- convert ENSG to gene symbol
str(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))
View(ens2sym)
# -- wrangle the resdf and join with ens2sym map
resdf <- resdf %>%
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-rowname) %>%
  inner_join(y = ens2sym, by = "GENEID")
View(resdf)
# -- volcano plot
EnhancedVolcano(resdf, lab = resdf$SYMBOL, pCutoff = 1e-100,
                FCcutoff = 3,
                x = "log2FoldChange", y = "padj")
# -- Get the over-expressed genes
resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange > 2) %>%
  write_csv(file = "over_expressed_genes.csv")
# Get the under-expressed genes
resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange < -2) %>%
  write_csv(file = "under_expressed_genes.csv")
 #https://maayanlab.cloud/Enrichr/  go to this site and copy and paste all the upregulated or downregulated genes
#this will give you which of all of these are transcription factors, or long non coding rna etc 
#are drugs that knock down ewsfli1 a good idea? based on enrichr maybe not because they seem to cause metastasis (pathway section)

## Hands-on activity #1 ##

# Find a study in recount and analyze it!
# You may be asked to share the results...

##########################


## DIY over-representation analysis ##


project_info <- abstract_search(query = "Ewing sarcoma")
project_info
selected_study2 <- "SRP028344"
download_study(selected_study2)
load("SRP028344/rse_gene.Rdata")

class(rse_gene)

cts2 <- assay(rse_gene)
View(cts2)

cd2 <- as.data.frame(colData(rse_gene))
View(cd2)

rse_gene$condition <- c(rep("shCTR", 3), rep("shERG", 4))
dds2 <- DESeqDataSet(rse_gene, design = ~condition)


# Get the over-expressed genes as a vector
over_expressed_genes <- resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange > 2) %>%
  pull(SYMBOL)

# Get the gene sets and wrangle
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Run over-representation analysis
egmt <- enricher(gene = over_expressed_genes,
                 TERM2GENE = gene_sets)
edf <- as.data.frame(egmt)
View(edf)

# Plot results with clusterProfiler
dotplot(egmt)
barplot(egmt)

## Hands-on activity #2 ##

# Continue your analysis from hands-on activity #1 by
# using enricher() to find the "curated" gene sets ("C2")

##########################

## GSEA ##

# Discrete cutoff -- or continuous distribution?
EnhancedVolcano(resdf, lab = resdf$SYMBOL,
                pCutoff = 1e-2,
                FCcutoff = 2,
                x = "log2FoldChange", y = "padj")

# Adding a score for GSEA
resdf2 <- resdf %>%
  arrange(padj) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Deal with inf
resdf2 <- resdf2 %>%
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Remove NAs and order by GSEA
resdf2 <- resdf2 %>%
  filter(! is.na(gsea_metric)) %>%
  arrange(desc(gsea_metric))
View(resdf2)

# GSEA value histogram
hist(resdf2$gsea_metric, breaks = 100)

# Get the ranked GSEA vector
ranks <- resdf2 %>%
  select(SYMBOL, gsea_metric) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  deframe()

# Run GSEA
gseares <- GSEA(geneList = ranks, 
                TERM2GENE = gene_sets)
gsearesdf <- as.data.frame(gseares)
View(gsearesdf)

# Make summary plots
dotplot(gseares)

# Make GSEA plot for "GO_CELL_MATRIX_ADHESION"
gseaplot(gseares, geneSetID = "GO_CELL_MATRIX_ADHESION",
         title = "GO_CELL_MATRIX_ADHESION")

# Make GSEA plot for "GO_RIBOSOMAL_SUBUNIT"
gseaplot(gseares, geneSetID = "GO_RIBOSOMAL_SUBUNIT",
         title = "GO_RIBOSOMAL_SUBUNIT")

# Make GSEA plot for top and bottom results
# -- Get top 4 over-expressed pathways
top_pathways <- gsearesdf %>%
  top_n(n = 4, wt = NES) %>%
  pull(ID)
# -- Make gseaplot for each and return as list
top_pathway_plots <- lapply(top_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})
# -- Arrange with labels as a multi-panel plot
top_pathway_plot <- ggarrange(plotlist = top_pathway_plots,
                              ncol = 2, nrow = 2, labels = "AUTO")
# -- Save it
ggsave(top_pathway_plot, filename = "top_GSEA_up.png",
       height = 11, width = 18)
# -- Repeat steps with top 4 under-expressed pathways
bottom_pathways <- gsearesdf %>%
  top_n(n = 4, wt = -NES) %>%
  pull(ID)
bottom_pathway_plots <- lapply(bottom_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})
bottom_pathway_plot <- ggarrange(plotlist = bottom_pathway_plots,
                                 ncol = 2, nrow = 2, labels = "AUTO")
ggsave(bottom_pathway_plot, filename = "top_GSEA_down.png",
       height = 11, width = 18)


## Hands-on activity #3 ##

# Continue your analysis from hands-on activity #2 by
# using GSEA() to find the over/under-expressed "curated" gene sets ("C2")

##########################



## Homework ##

# Continue to practice this...
# Can you address a biological question from your research?

##############





