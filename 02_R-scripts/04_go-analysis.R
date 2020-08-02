# RNAseq_DV_ILCs Project
# Part 4 : GO term analysis
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This script cover GO term enrichment analysis using gProfiler
# and cover code for plotting Figure 4C-F

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
library(gprofiler2)

# set up basic parameter
wk_dir = "../ILCs_DV_analysis/github/"
setwd(wk_dir)

save_dir = paste0(wk_dir,"rds/") # Directory for processed/intermediate rds files to be saved
data_dir = paste0(wk_dir,"data/") # Directory fot downloaded processed data

# ------------------------------------------
# Load data
# ------------------------------------------

tpmData <- read.delim(paste0(data_dir,"tpm_mtx.txt"))
head(tpmData)

metadata <- read.delim(paste0(data_dir,"metadata.txt"))
metadata

# results from previous part (03_deseq2.R)
res_DF_df <- readRDS(file = paste0(save_dir,"res_DF_df.rds"))
res_DHF_df <- readRDS(file = paste0(save_dir,"res_DHF_df.rds"))

# ------------------------------------------
# custom_gpr_plot
# ------------------------------------------

# this fx make barplots from gprofiler results
gostres_plot <- function(gostres, n.top = 10) {
  gostres_df <- gostres[["result"]] %>% dplyr::select(term_id,term_name,p_value,source,query,term_size,intersection_size) %>%
    mutate(log_p_value = -log10(p_value)) %>%
    mutate(signif = ifelse(p_value < 0.05, "p < 0.05","p >= 0.05")) %>%
    mutate(term_name = stringr::str_wrap(term_name, 80)) %>%
    mutate(int_sec_txt = paste0(intersection_size,"/",term_size)) %>%
    filter(source == "GO:BP") %>%
    filter(p_value < 1) %>%
    top_n(n.top, log_p_value)
  p <- ggplot(gostres_df, aes(x = reorder(term_name, log_p_value), y = log_p_value)) +
    geom_bar(stat = "identity", aes(fill = signif)) +
    geom_text(aes(label=int_sec_txt, y=0.05, hjust = 0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = list("p < 0.05" = "salmon","p >= 0.05" = "grey60")) +
    geom_hline(yintercept = c(1.3,2), linetype='dashed') +
    xlab(element_blank()) + ylab("-log10 p value") +
    theme_classic() +
    ggtitle(gostres_df$query, subtitle = "(intersection_size/term_size)") +
    coord_flip()
  print(p)
}

# ------------------------------------------
# Figure 4C,E
# ------------------------------------------

# make a list of DE genes in febrile phase of DF or DHF
DE_gene_list <- list(
  Febrile_DHF = res_DHF_df %>% filter(threshold == TRUE, log2FoldChange > 0) %>% .$gene %>% as.character(),
  Febrile_DF = res_DF_df %>% filter(threshold == TRUE, log2FoldChange > 0) %>% .$gene %>% as.character()
)

# Perform gProfiler analysis for DHF
gostres_DHF <- gost(DE_gene_list["Febrile_DHF"], organism = "hsapiens", ordered_query = TRUE,
                    multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
                    measure_underrepresentation = FALSE, evcodes = TRUE,
                    user_threshold = 0.05, custom_bg = NULL,
                    numeric_ns = "", sources = c("GO"))

# Print results
gostres_DHF[["result"]] %>%
  filter(source == "GO:BP") %>%
  dplyr::select(term_id, term_name, term_size, query_size, intersection_size, intersection, p_value, significant) %>%
  filter(significant == TRUE) %>%
  print(row.names = FALSE)

# Plot result (Figure 4E)
gostres_plot(gostres_DHF, n.top = 9)

# Perform gProfiler analysis for DF
gostres_DF <- gost(DE_gene_list["Febrile_DF"], organism = "hsapiens", ordered_query = TRUE,
                    multi_query = FALSE, significant = FALSE, exclude_iea = TRUE,
                    measure_underrepresentation = FALSE, evcodes = TRUE,
                    user_threshold = 0.05, custom_bg = NULL,
                    numeric_ns = "", sources = c("GO"))

# Print results
gostres_DF[["result"]] %>%
  filter(source == "GO:BP") %>%
  dplyr::select(term_id, term_name, term_size, query_size, intersection_size, intersection, p_value, significant) %>%
  filter(significant == TRUE) %>%
  print(row.names = FALSE)

# Plot result (Figure 4C)
gostres_plot(gostres_DF, n.top = 10)

# ------------------------------------------
# custom_gpr_plot 2
# ------------------------------------------

# this function take genes in enriched GO term and compare their expression fold-change between febrile and conv phase
# This fx can handle 1 GO term at a time by taking filtered result from gProfiler (describe below)
# row_to_plot parameter indicate which row (GO term) should this fx will plot
geneset_plot_l2fc <- function(emat, terms_sig, row_to_plot = 1){
  term_name <- paste0(terms_sig[row_to_plot,"term_id"]," ",terms_sig[row_to_plot,"term_name"])
  subtitle <- paste0("term_size = ",terms_sig[row_to_plot,"term_size"],
                     "   intersection_size = ",terms_sig[row_to_plot,"intersection_size"])
  gene_list <- terms_sig[row_to_plot,"intersection"] %>% str_split(pattern = ",") %>% unlist()

  # filter selected gene from emat and transpose sample into row
  df <- emat[emat$gene %in% gene_list, 2:14]

  #We add pseudo count here
  df[,2:13] <- df[,2:13]+0.1

  df$DHF1 <- log2(df$DHF1_F/df$DHF1_C)
  df$DHF2 <- log2(df$DHF2_F/df$DHF2_C)
  df$DHF3 <- log2(df$DHF3_F/df$DHF3_C)
  df$DF1 <- log2(df$DF1_F/df$DF1_C)
  df$DF2 <- log2(df$DF2_F/df$DF2_C)
  df$DF3 <- log2(df$DF3_F/df$DF3_C)

  df <- df[,c(1,14:19)]

  df <- df %>% pivot_longer(-gene, names_to = "patient", values_to = "L2FC")

  colData <- metadata %>% select(patient,severity) %>% distinct()
  df <- left_join(df, colData, by = "patient")
  df$gene_pt <- paste(df$gene,df$patient,sep = "_")

  p <- ggplot(df, aes(x = severity, y = L2FC, fill = severity)) +
    geom_violin(scale = "width", alpha = 0.6, width = 0.7) +
    geom_point(position = position_jitter(width = 0.2, seed = 123, height = 0), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    scale_fill_manual(values = list("DF" = "orange", "DHF" = "red")) +
    ggtitle(term_name, subtitle = subtitle) +
    theme_classic() +
    ylab("log2FC") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1.6)),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(size = rel(1.8), face ="bold"),
          axis.title.x = element_blank()
    )
  return(p)
}

# Filter only significant term from gProfileR result
terms_sig_DF <- gostres_DF[["result"]] %>%
  filter(source == "GO:BP") %>%
  dplyr::select(term_id, term_name, term_size, query_size, intersection_size, intersection, p_value, significant) %>%
  filter(significant == TRUE)
terms_sig_DF

# Plot those terms (Figure 4D)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DF, row_to_plot = 1)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DF, row_to_plot = 2)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DF, row_to_plot = 3)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DF, row_to_plot = 4)

# Filter only significant term from gProfileR result
terms_sig_DHF <- gostres_DHF[["result"]] %>%
  filter(source == "GO:BP") %>%
  dplyr::select(term_id, term_name, term_size, query_size, intersection_size, intersection, p_value, significant) %>%
  filter(significant == TRUE)
terms_sig_DHF

# Plot those terms (Figure 4F)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 1)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 2)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 3)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 4)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 5)
geneset_plot_l2fc(emat = tpmData, terms_sig = terms_sig_DHF, row_to_plot = 6)

# ------------------------------------------
# End of script
# ------------------------------------------
