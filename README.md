# ucsd_stat_genom

**TITLE**: Identifying the role of cell-type specific genetic regulatory variants in Multiple Sclerosis

**PROJECT DESCRIPTION**: Multiple sclerosis (MS) is a common autoimmune disease where the immune system attacks myelin, the protective sheath around nerve fibers in the central nervous system. This damage disrupts communication between neurons, leading to symptoms like weakness, numbness, and fatigue. With approximately 2.8 million people affected worldwide, MS is influenced by both genetic and environmental factors. Genome-wide association studies (GWAS) have identified over 200 genetic loci linked to MS, including variations in regulatory regions that influence gene expression.

In this project, we focus on single nucleotide polymorphisms (SNPs) located in enhancer regions, which regulate gene activity and can cause cell-type-specific changes in gene expression. Using GWAS data and computational tools, we aim to identify the genes and cell types affected by these enhancer SNPs to better understand their role in MS development. This research could uncover new insights into disease mechanisms and guide the development of potential therapies.

**DATA**: Our project utilizes three key datasets to investigate the genetic and regulatory factors involved in multiple sclerosis (MS). First, we use GWAS summary statistics from the study "Analysis of immune-related loci identifies 48 new susceptibility variants for multiple sclerosis," which analyzed genotyping data from 14,802 MS patients and 26,703 controls of European ancestry across 161,311 autosomal variants. 

Second, we incorporate H3K4me1 ChIP-seq data from a human-derived brain organoid (ENCODE Accession: [ENCSR844OQI](https://www.encodeproject.org/experiments/ENCSR844OQI/)). This dataset provides genome-wide regulatory annotations, including active and poised enhancer regions marked by H3K4me1. Brain organoids were chosen because they represent a range of cell types found in the adult brain (e.g., neurons, microglia) and offer a practical alternative to sampling adult brain tissue directly. By intersecting GWAS-identified SNPs with these enhancer regions, we prioritize variants with potential regulatory roles in MS-relevant cell types.

Lastly, we use genotyping and gene expression data from the 1000 Genomes Project for 358 individuals of European ancestry to perform eQTL analysis. This dataset includes genotyping information for 155,756 SNPs and expression data for 23,722 genes. While the eQTL cohort is different from the GWAS cohort, both datasets consist of large European populations, minimizing any potential bias. By integrating these datasets, we connect SNPs, enhancer activity, and gene expression to uncover their potential roles in MS.

**HOW TO USE OUR CODE?**
The code for the eQTL analysis can be found under /eqtl/FinalProject_eQTL.R
The following libraries are required to run eQTL analysis:
- data.table
- plink2R
- dplyr
- ggplot2
- biomaRt
- writexl

The code for IMPACT can be found under /impact/project.R
The following libraries are required to run IMPACT:
- data.table
- dplyr
- plink2R
- GenomicRanges
- Gviz
- rtracklayer
