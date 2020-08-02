# RNAseq_DV_ILCs Project
# Part 3 : DEseq2
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This script cover DESeq2 analysis and cover figure 4A,B

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
library(ggrepel)
library(eulerr)

# set up basic parameter
wk_dir = "../ILCs_DV_analysis/github/"
setwd(wk_dir)

save_dir = paste0(wk_dir,"rds/") # Directory for processed/intermediate rds files to be saved
data_dir = paste0(wk_dir,"data/") # Directory fot downloaded processed data

n_workers = 4 #no. of worker for paralellization

# set up significant threshold for differentially expressed genes
threshold_log2fc = 1 #Log2 fold change
threshold_padj = 0.01 #adjusted p value

# ------------------------------------------
# Load data
# ------------------------------------------

countData <- read.delim(paste0(data_dir,"count_mtx.txt"))
head(countData)

tpmData <- read.delim(paste0(data_dir,"tpm_mtx.txt"))
head(tpmData)

metadata <- read.delim(paste0(data_dir,"metadata.txt"))
metadata

# Subset data to only DF samples
countData_DF <- countData %>% select(ENSG,gene,starts_with("DF"))
metadata_DF <- metadata %>% filter(severity == "DF")

# Subset data to only DHF samples
countData_DHF <- countData %>% select(ENSG,gene,starts_with("DHF"))
metadata_DHF <- metadata %>% filter(severity == "DHF")

# ------------------------------------------
# custom fx
# ------------------------------------------

# This fx process result from DESeq2 into dataframe and add gene symbol to the dataframe
res_to_dataframe <- function(dds.res = res_DF, tpm.mtx = tpmData){
  res.name <- dds.res@elementMetadata@listData[["description"]][2] %>% substring(25)
  res.df <- as.data.frame(dds.res) %>% tibble::rownames_to_column("ENSG")
  ensg_hgnc <- tpm.mtx %>% filter(ENSG %in% res.df$ENSG) %>% dplyr::select(ENSG, gene)
  res.df <- merge(res.df,ensg_hgnc)
  res.df$resname <- res_name
  rownames(res.df) <- res.df$ENSG
  return(res.df)
}

# This fx take processed dataframe from fx res_to_dataframe and flag significant DE genes according to cutoff setting
resdf_prep_volc <- function(res.df, threshold.log2fc = threshold_log2fc, threshold.padj = threshold_log2fc){
  res.tmp <- as.data.frame(res.df)
  res.tmp$log10p_val <- -log10(res.tmp$pvalue)
  res.tmp$log10p_val_adj <- -log10(res.tmp$padj)
  res.tmp$threshold = as.factor(abs(res.tmp$log2FoldChange) > threshold.log2fc & res.tmp$padj < threshold.padj)
  return(res.tmp)
}

# This fx take dataframe output from resdf_prep_volc to make a volcano plot
# genes to label can assign via gene_plot parameters
volc_plot <- function(res.tmp, gene_plot = c("CD69","KLRB1"), xlim = c(-25,25)){
  p <- ggplot(res.tmp, aes(x=log2FoldChange, y=log10p_val_adj, label = gene)) +
    geom_point(aes(colour = threshold, text = gene )) +
    geom_text_repel(
      data = subset(res.tmp, res.tmp$gene %in% gene_plot)
    ) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
    coord_cartesian(xlim = xlim) +
    theme_light() +
    ggtitle(res.tmp[1,]$resname,
            paste0("LFC cutoff = ",threshold_log2fc,", p adj cutoff = ",threshold_padj)
    )
  # p <- ggplotly(p)
  return(p)
}

# ------------------------------------------
# DESeq2 DF
# ------------------------------------------

# create DESeq2 obj
dds_DF <- DESeqDataSetFromMatrix(countData = countData_DF[,-(1:2)],
                              colData = metadata_DF,
                              design = ~ patient + timepoint)

# Clean up data
keep <- rowSums(counts(dds_DF)==0) <= 3
dds_DF <- dds_DF[keep,]

# Perform DESeq2 analysis
dds_DF <- DESeq(dds_DF, parallel = TRUE, BPPARAM=MulticoreParam(workers = n_workers))

# get result and sort by padj
res_DF <- results(dds_DF, name="timepoint_Febrile_vs_Conv", alpha=0.05)
res_DF <- res_DF[order(res_DF$padj),]

# Perform apeglm shrinkage estimation
res_DF_apeglm <- lfcShrink(dds_DF, coef="timepoint_Febrile_vs_Conv", type="apeglm", parallel = TRUE)
res_DF_apeglm <- res_DF_apeglm[order(res_DF_apeglm$padj),]

# get gene symbol
res_DF_df <- res_to_dataframe(dds.res = res_DF, tpm.mtx = tpmData)

#volc plot
res_DF_df <- resdf_prep_volc(res_DF_df, threshold.log2fc = threshold_log2fc, threshold.padj = threshold_padj)
gene_to_plot <- res_DF_df %>% filter(threshold == TRUE) %>% arrange(-log10p_val_adj) %>% .$gene %>% as.character()
volc_plot(res_DF_df, gene_to_plot[1:50])
# volc_plot(res_DF_df, "CD69")

# get gene symbol
res_DF_apeglm_df <- res_to_dataframe(dds.res = res_DF_apeglm, tpm.mtx = tpmData)

# volc plot
res_DF_apeglm_df <- resdf_prep_volc(res_DF_apeglm_df, threshold.log2fc = threshold_log2fc, threshold.padj = threshold_padj)
gene_to_plot <- res_DF_apeglm_df %>% filter(threshold == TRUE) %>% arrange(-log10p_val_adj) %>% .$gene %>% as.character()
volc_plot(res_DF_apeglm_df, gene_to_plot[1:50])
# volc_plot(res_DF_apeglm_df, "CD69")

# ------------------------------------------
# DESeq2 DHF
# ------------------------------------------

dds_DHF <- DESeqDataSetFromMatrix(countData = countData_DHF[,-(1:2)],
                                 colData = metadata_DHF,
                                 design = ~ patient + timepoint)

keep <- rowSums(counts(dds_DHF)==0) <= 3
dds_DHF <- dds_DHF[keep,]
dds_DHF

dds_DHF <- DESeq(dds_DHF, parallel = TRUE, BPPARAM=MulticoreParam(workers = n_workers))
resultsNames(dds_DHF)

res_DHF <- results(dds_DHF, name="timepoint_Febrile_vs_Conv", alpha=0.05)
res_DHF <- res_DHF[order(res_DHF$padj),]

res_DHF_apeglm <- lfcShrink(dds_DHF, coef="timepoint_Febrile_vs_Conv", type="apeglm", parallel = TRUE)
res_DHF_apeglm <- res_DHF_apeglm[order(res_DHF_apeglm$padj),]

res_DHF_df <- res_to_dataframe(dds.res = res_DHF, tpm.mtx = tpmData)

res_DHF_df <- resdf_prep_volc(res_DHF_df, threshold.log2fc = threshold_log2fc, threshold.padj = threshold_padj)
gene_to_plot <- res_DHF_df %>% filter(threshold == TRUE) %>% arrange(-log10p_val_adj) %>% .$gene %>% as.character()
volc_plot(res_DHF_df, gene_to_plot[1:50])
# volc_plot(res_DHF_df, "CD69")
# res_DHF_df %>% filter(gene == "CD69")

res_DHF_apeglm_df <- res_to_dataframe(dds.res = res_DHF_apeglm, tpm.mtx = tpmData)

res_DHF_apeglm_df <- resdf_prep_volc(res_DHF_apeglm_df, threshold.log2fc = threshold_log2fc, threshold.padj = threshold_padj)
gene_to_plot <- res_DHF_apeglm_df %>% filter(threshold == TRUE) %>% arrange(-log10p_val_adj) %>% .$gene %>% as.character()
volc_plot(res_DHF_apeglm_df, gene_to_plot[1:50])
# volc_plot(res_DHF_apeglm_df, "CD69")
# res_DHF_apeglm_df %>% filter(gene == "CD69")

# ------------------------------------------
# LFC plot Figure 4A
# ------------------------------------------

# get significant DE genes in either DF or DHF
sig_DHF <- res_DHF_df %>% filter(threshold == TRUE) %>% dplyr::select(log2FoldChange, gene)
colnames(sig_DHF) <- c("DHF_LFC","gene")

sig_DF <- res_DF_df %>% filter(threshold == TRUE) %>% dplyr::select(log2FoldChange, gene)
colnames(sig_DF) <- c("DF_LFC","gene")

sig_LFC <- inner_join(sig_DHF,sig_DF)

# plot only significant DE genes
ggplot(sig_LFC, aes(x = DHF_LFC, y =DF_LFC, label = gene)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(size = 3, colour = "#F8766D") + geom_text_repel() +
  labs(title = "DEG Fold Change: DHF vs DF (only both True)") +
  theme_minimal()

# get all DE genes in either DF or DHF
DHF_all <- res_DHF_df %>%
  dplyr::select(log2FoldChange, threshold, log10p_val_adj, gene, ENSG)
colnames(DHF_all) <- c("DHF_LFC","DHF_Threshold","DHF_log10padj","genes", "ENSG")

DF_all <- res_DF_df %>%
  dplyr::select(log2FoldChange, threshold, log10p_val_adj, gene, ENSG)
colnames(DF_all) <- c("DF_LFC","DF_Threshold","DF_log10padj","genes", "ENSG")

# merge table and triage significance in either or both conditions
all_LFC <- inner_join(DHF_all,DF_all)
all_LFC <- all_LFC %>%
  mutate(new_threshold = case_when(
    DHF_Threshold == TRUE & DF_Threshold == TRUE ~ "Both True",
    DHF_Threshold == TRUE & DF_Threshold == FALSE ~ "DHF True",
    DHF_Threshold == FALSE & DF_Threshold == TRUE ~ "DF True",
    TRUE ~ "False"
  ))

all_LFC$new_threshold <- all_LFC$new_threshold %>% as.factor()
gene_label <- sig_LFC$gene %>% as.character()

# Plot Figure 4A
ggplot(all_LFC, aes(x = DHF_LFC, y =DF_LFC, label = genes)) +
  geom_point(aes(alpha = new_threshold, colour = new_threshold, size = new_threshold)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_rug(alpha = 0.1, size = 0.1) +
  scale_color_manual(values=c("#F8766D","#7CAE00","#00BFC4","lightgrey"),
                     name = "Significance",
                     labels = c("DHF & DF", "DF", "DHF","None")) +
  scale_alpha_manual(values=c(1,1,1,0.4),
                     name = "Significance",
                     labels = c("DHF & DF", "DF", "DHF","None")) +
  scale_size_manual(values=c(4,2,2,1),
                    name = "Significance",
                    labels = c("DHF & DF", "DF", "DHF","None")) +
  labs(title = "DEG Fold Change: DHF vs DF",
       subtitle = paste0("Log2FC threshold > ",threshold_log2fc,", padj threshold < ",threshold_padj)) +
  geom_text_repel(
    data = subset(all_LFC, all_LFC$genes %in% gene_label ),
    size = 4.5
  ) +
  theme_classic() +
  xlab("DHF Log2 Fold Change") +
  ylab("DF Log2 Fold Change") +
  theme(axis.title = element_text(size = rel(1.6), face = "bold"),
        axis.text = element_text(size = rel(1.6)),
        axis.line = element_line(size = rel(1.6)),
        legend.title = element_text(size = rel(1.6), face = "bold"),
        legend.text = element_text(size = rel(1.6))
  )

# ------------------------------------------
# Venn diagram Fig 4B
# ------------------------------------------

# filter only significant DE genes
Ac_DHF <- res_DHF_df %>% filter(log2FoldChange >0) %>% filter(threshold == TRUE)
Ac_DF <- res_DF_df %>% filter(log2FoldChange >0) %>% filter(threshold == TRUE)

# plot Venn diagram
combined <- list()
combined[["Febrile DF"]] <- Ac_DF[["gene"]] %>% as.character()
combined[["Febrile DHF"]] <- Ac_DHF[["gene"]] %>% as.character()
fit <- eulerr::euler(combined)
plot(fit, quantities = TRUE)

# ------------------------------------------
# save files
# ------------------------------------------

saveRDS(res_DF, file = paste0(save_dir,"res_DF.rds"))
saveRDS(res_DF_df, file = paste0(save_dir,"res_DF_df.rds"))
saveRDS(res_DF_apeglm, file = paste0(save_dir,"res_DF_apeglm.rds"))

saveRDS(res_DHF, file = paste0(save_dir,"res_DHF.rds"))
saveRDS(res_DHF_df, file = paste0(save_dir,"res_DHF_df.rds"))
saveRDS(res_DHF_apeglm, file = paste0(save_dir,"res_DHF_apeglm"))

# ------------------------------------------
# End of script
# ------------------------------------------
