## Description

This repository contains several scripts necessary to perform RNA-sequencing analysis in this publication  *"Innate Lymphoid Cells Activation and Transcriptional Changes in Human Dengue Infection"*.

Primary purpose of this repository is to provide additional details for method of analyses performed in the publication.

## Data
Data used in this publication was deposited under ancession number GSEXXXXXXX
fastq files were deposited on SRA.
Processed files (TPM, Read counts) of all samples are also available

## Raw data processing
This part demonstrate how fastq files were processed
After FastQC evaluation, we trimmed Nextera adapter using Trimmomatic.
Read mapping was done using HISAT2
TPM was obtained by StingTie.
Read counts were obtain by HTseq-count
Processed data of all samples were deposited on GEO (see above)

## R scripts
R scripts for data analysis in Fig 3 and 4 in the publication
* Part 02 demonstrate data visualization in Figure 3
* Part 03 perform DESeq2 analysis and data visualization in Figure 4A,B
* Part 04 focus on GO term enrichment analysis with gProfileR and data disualization in Figure 4C-F

## Packages used in this analysis

[Trimmomatic (0.36)](http://www.usadellab.org/cms/?page=trimmomatic)

[HISAT2 (2.1.0)](http://daehwankimlab.github.io/hisat2/)

[StringTie (1.3.5)](https://ccb.jhu.edu/software/stringtie/)

[HTseq-count (0.6.1p1)](https://htseq.readthedocs.io/)

R 3.6.0
* [ComplexHeatmap (2.2.0)](https://github.com/jokergoo/ComplexHeatmap)
* [dplyr (0.8.5)](https://dplyr.tidyverse.org/)
* [stringr (1.4.0)](https://stringr.tidyverse.org/)
* [ggplot2 (3.3.0)](https://ggplot2.tidyverse.org/)
* [ggrepel (0.8.2)](https://github.com/slowkow/ggrepel)
* [tidyr (1.0.2)](https://tidyr.tidyverse.org/)
* [tibble (3.0.1)](https://tibble.tidyverse.org/)
* [gprofiler2 (0.1.8)](https://biit.cs.ut.ee/gprofiler/page/r)

Data analyses were performed on Ubuntu 16.04.6 LTS
