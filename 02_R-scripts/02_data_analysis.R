# RNAseq_DV_ILCs Project
# Part 2 : data analysis
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This script cover data visualization in Figure 3

# ------------------------------------------
# Load required libraries and set up basic parameters
# ------------------------------------------

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(stringr)

# set up basic parameter
wk_dir = "../ILCs_DV_analysis/github/"
setwd(wk_dir)

save_dir = paste0(wk_dir,"rds/") # Directory for processed/intermediate rds files to be saved
data_dir = paste0(wk_dir,"data/") # Directory fot downloaded processed data

# ------------------------------------------
# Load data
# ------------------------------------------

countData <- read.delim(paste0(data_dir,"count_mtx.txt"))
head(countData)

tpmData <- read.delim(paste0(data_dir,"tpm_mtx.txt"))
head(tpmData)

metadata <- read.delim(paste0(data_dir,"metadata.txt"))
metadata

# ------------------------------------------
# Figure 3B
# ------------------------------------------

# First, we list genes of interested
gene_plot <- c("PTPRC","LTB","IL7R","KLRB1","TBX21","PTGDR2","GATA3","AHR","KIT",
               "CD3E","CD4","CD8A","MS4A1","CD79B","NCAM1")

col_orders <- c("DF1_F","DF2_F","DF3_F","DF1_C","DF2_C","DF3_C",
                "DHF1_F","DHF2_F","DHF3_F","DHF1_C","DHF2_C","DHF3_C")

# filter expression matrix to only genes of interested
emat <- tpmData %>%
  dplyr::select(-ENSG) %>%
  filter(gene %in% gene_plot) %>%
  filter(rowSums(.[,-1]) > 0) %>%
  tibble::column_to_rownames(var = "gene")
emat
emat <- log1p(emat)

# heatmap annotation
ha = HeatmapAnnotation(
  severity = metadata$severity,
  timepoint = metadata$timepoint,
  col = list(
    severity = c("DHF"="red","DF"="orange"),
    timepoint = c("Febrile"="#B2182B","Conv"="#92C5DE")
  ),
  simple_anno_size = unit(0.3, "cm")
)

# Heatmap plot using ComplexHeatmap
Heatmap(emat,
        name = "log1p TPM\n(non z-score)",
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "white", lwd = 0.5),
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_side = "top",
        column_order = col_orders,
        row_order = gene_plot,
        top_annotation = ha
)

# ------------------------------------------
# Figure 3D
# ------------------------------------------

# calculate pearson correlation
corm = cor(tpmData[,3:14], method = "pearson")

# heatmap annotation
ha = HeatmapAnnotation(
  severity = metadata$severity,
  timepoint = metadata$timepoint,
  col = list(
    severity = c("DHF"="red","DF"="orange"),
    timepoint = c("Febrile"="#B2182B","Conv"="#92C5DE")
  ),
  simple_anno_size = unit(0.3, "cm")
)

# Heatmap plot using ComplexHeatmap
Heatmap(corm,
        name = "correlation",
        column_title = "Pearson Correlation",
        rect_gp = gpar(col = "white", lwd = 0.5),
        col = circlize::colorRamp2(c(0.89, 1), c("white", "red")),
        show_row_dend = F,
        top_annotation = ha)

# ------------------------------------------
# Figure 3E
# ------------------------------------------

# we first filtered out contaminated genes
# Since IG genes have nothing to do with ILCs
# We hypothesized that during cell isolation and preparation,
# damaged cells released their internal mRNA content thus causing background/ambient RNA contamination

# get IG genes list then filtered out from expression matrix
IGH_genes <- countData$gene %>% stringr::str_subset(pattern = "^IGH")
countData <- countData[!(countData$gene %in% IGH_genes),]

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countData[,3:14],
                              colData = metadata,
                              design = ~ patient + timepoint)

# Clean up non/lowly expressed genes
keep <- rowSums(counts(dds)==0) <= 6
dds <- dds[keep,]

#vsd normalization
vsd <- vst(dds)

# set some params here
intgroup = "condition" #match with metadata col
ntop = 200 #number of top variable genes used
arrow_multiplier = 20 #loading size on plot

# calculate the variance for each gene
rv <- rowVars(assay(vsd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

intgroup.df <- as.data.frame(colData(vsd)[,intgroup, drop=FALSE])

# add the intgroup factors together to create a new grouping factor
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(vsd)[[intgroup]]
}

# merge data for the plot
df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(vsd))
col_list <- list("Febrile_DF" = "orange", "Febrile_DHF" = "red", "Conv_DF" = "lightblue", "Conv_DHF" = "dodgerblue")

# get gene loadings for PC1 and PC2
loadings <- as.data.frame(pca$rotation)
loadings <- loadings %>% tibble::rownames_to_column("ENSG")
loadings <- left_join(loadings, countData %>% select(ENSG, gene), by = "ENSG")
loadings_genes <- c(loadings %>% select(gene, PC1) %>% arrange(-PC1) %>% top_n(5, PC1) %>% .$gene %>% as.character(),
                    loadings %>% select(gene, PC2) %>% arrange(-PC2) %>% top_n(5, PC2) %>% .$gene %>% as.character(),
                    loadings %>% select(gene, PC1) %>% arrange(-PC1) %>% top_n(5, -PC1) %>% .$gene %>% as.character(),
                    loadings %>% select(gene, PC2) %>% arrange(-PC2) %>% top_n(5, -PC2) %>% .$gene %>% as.character())
loading_plot <- loadings %>% filter(gene %in% loadings_genes)

# Plot PCA with loadings
ggplot(data=df, aes_string(x="PC1", y="PC2", color="group")) +
  geom_segment(data = loading_plot,
               aes(x = 0, y = 0,
                   xend = (PC1*arrow_multiplier),
                   yend = (PC2*arrow_multiplier)), arrow = arrow(length = unit(1/2, "picas")),
               color = "lightgrey") +
  annotate("text",
           x = (loading_plot$PC1*arrow_multiplier*1.1),
           y = (loading_plot$PC2*arrow_multiplier*1.1),
           label = loading_plot$gene,
           color = "grey30") +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  ggtitle("PCA") +
  scale_color_manual(values = col_list) +
  theme_classic()

# ------------------------------------------
# End of script
# ------------------------------------------
