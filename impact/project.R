# Specifying root directory
knitr::opts_knit$set(root.dir = '/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project')

# Loading libraries
library(data.table)
library(plink2R)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(Sushi)
library("devtools")



# Loading GWAS summary statistics
gwas_file <- "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv"
gwas_data <- fread(gwas_file, sep = "\t", header = TRUE)

#Exploring data
head(gwas_data)
str(gwas_data)

#Getting significant SNPs
significant_snps <- gwas_data %>%
  filter(p < 5e-8)

write.csv(significant_snps, "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/significant_snps.csv", row.names = FALSE)

# Create GRanges object from SNP data
snps_gr <- GRanges(
  seqnames = as.factor(significant_snps$chrom),  # Chromosome
  ranges = IRanges(start = significant_snps$pos, end = significant_snps$pos), # SNP Position
  mcols = significant_snps                    # Metadata (e.g., rsID, p-value, effect size)
)

# Loading enhancer region data
enhancer_data <- fread("/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/ENCFF026GBE.bed")
# match the format od chromosome column to SNP data
enhancer_data$V1 <- gsub("^chr", "", enhancer_data$V1)
colnames(enhancer_data) <- c("chrom", "start", "end", "V4", "V5", "V6", "V7", "V8", "V9", "V10")
enhancer_gr <- GRanges(
  seqnames = as.factor(enhancer_data$chrom),    
  ranges = IRanges(start = enhancer_data$start,
                   end = enhancer_data$end)   
)

# Find overlaps between SNPs and enhancers
overlaps <- findOverlaps(snps_gr, enhancer_gr)

# extract overlapping SNPs and annotate
overlapping_snps <- significant_snps[queryHits(overlaps), ]
write.csv(overlapping_snps, "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/enhancer_mapped_snps.csv", row.names = FALSE)

# Adding start and end positions of enhancers associated with SNPs to the dataframe
enhancer_start <- start(enhancer_gr[subjectHits(overlaps)])
enhancer_end <- end(enhancer_gr[subjectHits(overlaps)])

overlapping_snps$enhancer_start <- enhancer_start
overlapping_snps$enhancer_end <- enhancer_end


#In the other part of our project, we identified the genes associated with significant SNP's. I will use those genes to filter my SNPs.
#The reason why is because the Galaxy database contains bedgraph files seperately for each chromosome and transcription factor. Due to time and storage, I can't download and go through all these files.

#first getting the unique SNP's associated with a gene through eQTL
unique_eqtl_snps <- unique(eqtl_enhancer$SNP)

# Now I'm comparing the SNP's that I found that are in enhancer regions with these unique SNP's and intersecting them
snps_associated_w_genes<- overlapping_snps[overlapping_snps$rsid %in% unique_eqtl_snps, ]

cat(unique(snps_associated_w_genes$chrom)) #we only have SNP's from chromosomes 1,2,6,7,8,10,11,12,16,19 so using only those bedgraph files from Galaxy

# We chose the following transcription factors to look into: SOX2, CTCF, MEF2A, RXRA
# And the cell types: brain, blood, stem cell, eye

# a function to read and process the bedgraph data
process_bedgraph <- function(file_path, chrom) {
  data <- fread(file_path)
  data <- cbind(chrom, data)
  colnames(data) <- c("chrom", "start", "end", "value")
  return(data)
}

chrom_bounds <- snps_associated_w_genes %>%
  group_by(chrom) %>%
  summarise(
    min_pos = min(pos, na.rm = TRUE), # left-most position
    max_pos = max(pos, na.rm = TRUE)  # right-most position
  )

#Due to time limitations, we're only going to plot the top 3 chromosomes with the most SNPs --> 6,1,19
top_3_chromosomes <- snps_associated_w_genes %>%
  group_by(chrom) %>%
  summarise(snp_count = n()) %>%
  arrange(desc(snp_count)) %>%
  slice_head(n = 3)

#function for plotting
plot_bedgraph_by_tissue <- function(bedgraph_list, tissue_labels, chrom, chromstart, chromend, color_map, title) {
  # Subset tissue labels and colors for available data
  available_tissues <- names(bedgraph_list)
  available_colors <- color_map[available_tissues]
  
  # Plot each bedgraph
  for (i in seq_along(bedgraph_list)) {
    plotBedgraph(
      bedgraph_list[[i]],
      chrom = chrom,
      chromstart = chromstart,
      chromend = chromend,
      transparency = 0.50,
      overlay = TRUE,
      rescaleoverlay = FALSE,
      color = available_colors[i]
    )
    labelgenome(
      chrom = chrom,
      chromstart = chromstart,
      chromend = chromend,
      n = 4,
      scale = "Mb"
    )
    axis(side = 2, las = 2, tcl = 0.2)
  }
  
  # Add legend for available tissues
  legend(
    "topright",
    inset = 0.025,
    legend = tissue_labels[available_tissues],
    fill = available_colors,
    border = "black",
    text.font = 2,
    cex = 0.5
  )
  
  
  if (!is.null(title)) {
    title(main = title)
  }
}

# Colors for each tissue
color_map <- c(
  brain = "green",
  stem_cell = "purple",
  blood = "blue",
  eye = "orange"
)

# Tissue labels
tissue_labels <- c(
  brain = "Brain",
  stem_cell = "Stem Cell",
  blood = "Blood",
  eye = "Eye"
)

#LOADING THE DATA FOR CHR1
#BRAIN
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy736-[IMPACT_211_predictions_SOX2_chr1.bedgraph.gz",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy670-[IMPACT_529_predictions_CTCF_chr1.bedgraph.gz",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy538-[IMPACT_483_predictions_MEF2A_chr1.bedgraph.gz",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy494-[IMPACT_464_predictions_RXRA_chr1.bedgraph.gz"
)
brain_ch1_data <- lapply(file_paths, process_bedgraph, chrom = "chr1")
brain_ch1_sox2 <- brain_ch1_data[["SOX2"]]
brain_ch1_ctcf <- brain_ch1_data[["CTCF"]]
brain_ch1_mef2a <- brain_ch1_data[["MEF2A"]]
brain_ch1_rxra <- brain_ch1_data[["RXRA"]]

#BLOOD
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy4851-[IMPACT_391_predictions_CTCF_chr1.bedgraph.gz]",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy4037-[IMPACT_584_predictions_MEF2A_chr1.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy3575-[IMPACT_3_predictions_RXRA_chr1.bedgraph.gz]"
)
blood_ch1_data <- lapply(file_paths, process_bedgraph, chrom = "chr1")
blood_ch1_ctcf <- blood_ch1_data[["CTCF"]]
blood_ch1_mef2a <- blood_ch1_data[["MEF2A"]]
blood_ch1_rxra <- blood_ch1_data[["RXRA"]]

#STEM CELL
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy561-[IMPACT_370_predictions_SOX2_chr1.bedgraph.gz]",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy1419-[IMPACT_339_predictions_CTCF_chr1.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy1045-[IMPACT_618_predictions_RXRA_chr1.bedgraph.gz]"
)
sc_ch1_data <- lapply(file_paths, process_bedgraph, chrom = "chr1")
sc_ch1_sox2 <- sc_ch1_data[["SOX2"]]
sc_ch1_ctcf <- sc_ch1_data[["CTCF"]]
sc_ch1_rxra <- sc_ch1_data[["RXRA"]]

#EYE
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/Galaxy33-[IMPACT_327_predictions_CTCF_chr1.bedgraph.gz]"
)
eye_ch1_data <- lapply(file_paths, process_bedgraph, chrom = "chr1")
eye_ch1_ctcf <- eye_ch1_data[["CTCF"]]


#PLOTTING FOR CHROMOSOME 1
#first we need to find the left-most and right-most SNP's in chromosome 1
chrom1_bounds <- chrom_bounds %>% filter(chrom == 1)
lowest_chrom1 <- chrom1_bounds$min_pos
highest_chrom1 <- chrom1_bounds$max_pos

# Data for different tissues
bedgraph_list_ctcf <- list(
  brain = brain_ch1_ctcf,
  stem_cell = sc_ch1_ctcf,
  blood = blood_ch1_ctcf,
  eye = eye_ch1_ctcf
)

bedgraph_list_mef2a <- list(
  brain = brain_ch1_mef2a,
  blood = blood_ch1_mef2a
)

bedgraph_list_sox2 <- list(
  brain = brain_ch1_sox2,
  stem_cell = sc_ch1_sox2
)

bedgraph_list_rxra <- list(
  brain = brain_ch1_rxra,
  blood = blood_ch1_rxra,
  stem_cell = sc_ch1_rxra
)

chrom <- "chr1"
# Plot for CTCF
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_ctcf,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom1,
  chromend = highest_chrom1,
  color_map = color_map,
  title = "CTCF Transcription Factor"
)

# Plot for MEF2A
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_mef2a,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom1,
  chromend = highest_chrom1,
  color_map = color_map,
  title = "MEF2A Transcription Factor"
)

# Plot for SOX2
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_sox2,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom1,
  chromend = highest_chrom1,
  color_map = color_map,
  title = "SOX2 Transcription Factor"
)

#Plot for RXRA
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_rxra,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom1,
  chromend = highest_chrom1,
  color_map = color_map,
  title = "RXRA Transcription Factor"
)


#LOADING THE DATA FOR CHR6

#BRAIN
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy744-[IMPACT_211_predictions_SOX2_chr6.bedgraph.gz]",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy678-[IMPACT_529_predictions_CTCF_chr6.bedgraph.gz]",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy546-[IMPACT_483_predictions_MEF2A_chr6.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy502-[IMPACT_464_predictions_RXRA_chr6.bedgraph.gz]"
)
brain_ch6_data <- lapply(file_paths, process_bedgraph, chrom = "chr6")
brain_ch6_sox2 <- brain_ch6_data[["SOX2"]]
brain_ch6_ctcf <- brain_ch6_data[["CTCF"]]
brain_ch6_mef2a <- brain_ch6_data[["MEF2A"]]
brain_ch6_rxra <- brain_ch6_data[["RXRA"]]

#BLOOD
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy4859-[IMPACT_391_predictions_CTCF_chr6.bedgraph.gz]",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy4045-[IMPACT_584_predictions_MEF2A_chr6.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy3583-[IMPACT_3_predictions_RXRA_chr6.bedgraph.gz]"
)
blood_ch6_data <- lapply(file_paths, process_bedgraph, chrom = "chr6")
blood_ch6_ctcf <- blood_ch6_data[["CTCF"]]
blood_ch6_mef2a <- blood_ch6_data[["MEF2A"]]
blood_ch6_rxra <- blood_ch6_data[["RXRA"]]

#STEM CELL
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy1713-[IMPACT_369_predictions_SOX2_chr6.bedgraph.gz]",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy1427-[IMPACT_339_predictions_CTCF_chr6.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy1053-[IMPACT_618_predictions_RXRA_chr6.bedgraph.gz]"
)
sc_ch6_data <- lapply(file_paths, process_bedgraph, chrom = "chr6")
sc_ch6_sox2 <- sc_ch6_data[["SOX2"]]
sc_ch6_ctcf <- sc_ch6_data[["CTCF"]]
sc_ch6_rxra <- sc_ch6_data[["RXRA"]]

#EYE
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr6/Galaxy41-[IMPACT_327_predictions_CTCF_chr6.bedgraph.gz]"
)
eye_ch6_data <- lapply(file_paths, process_bedgraph, chrom = "chr6")
eye_ch6_ctcf <- eye_ch6_data[["CTCF"]]

#PLOTTING FOR CHROMOSOME 6
chrom6_bounds <- chrom_bounds %>% filter(chrom == 6)
lowest_chrom6 <- chrom6_bounds$min_pos
highest_chrom6 <- chrom6_bounds$max_pos

# Data for different tissues
bedgraph_list_ctcf <- list(
  brain = brain_ch6_ctcf,
  stem_cell = sc_ch6_ctcf,
  blood = blood_ch6_ctcf,
  eye = eye_ch6_ctcf
)

bedgraph_list_mef2a <- list(
  brain = brain_ch6_mef2a,
  blood = blood_ch6_mef2a
)

bedgraph_list_sox2 <- list(
  brain = brain_ch6_sox2,
  stem_cell = sc_ch6_sox2
)

bedgraph_list_rxra <- list(
  brain = brain_ch6_rxra,
  blood = blood_ch6_rxra,
  stem_cell = sc_ch6_rxra
)


chrom <- "chr6"
# Plot for CTCF
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_ctcf,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom6,
  chromend = highest_chrom6,
  color_map = color_map,
  title = "CTCF Transcription Factor"
)

# Plot for MEF2A
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_mef2a,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom6,
  chromend = highest_chrom6,
  color_map = color_map,
  title = "MEF2A Transcription Factor"
)

# Plot for SOX2
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_sox2,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom6,
  chromend = highest_chrom6,
  color_map = color_map,
  title = "SOX2 Transcription Factor"
)

#Plot for RXRA
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_rxra,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom6,
  chromend = highest_chrom6,
  color_map = color_map,
  title = "RXRA Transcription Factor"
)

#LOADING THE DATA FOR CHR19
#BRAIN
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy735-[IMPACT_211_predictions_SOX2_chr19.bedgraph.gz]",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy669-[IMPACT_529_predictions_CTCF_chr19.bedgraph.gz]",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy537-[IMPACT_483_predictions_MEF2A_chr19.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy493-[IMPACT_464_predictions_RXRA_chr19.bedgraph.gz]"
)
brain_ch19_data <- lapply(file_paths, process_bedgraph, chrom = "chr19")
brain_ch19_sox2 <- brain_ch19_data[["SOX2"]]
brain_ch19_ctcf <- brain_ch19_data[["CTCF"]]
brain_ch19_mef2a <- brain_ch19_data[["MEF2A"]]
brain_ch19_rxra <- brain_ch19_data[["RXRA"]]

#BLOOD
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy4850-[IMPACT_391_predictions_CTCF_chr19.bedgraph.gz]",
  "MEF2A" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy4036-[IMPACT_584_predictions_MEF2A_chr19.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy3574-[IMPACT_3_predictions_RXRA_chr19.bedgraph.gz]"
)
blood_ch19_data <- lapply(file_paths, process_bedgraph, chrom = "chr19")
blood_ch19_ctcf <- blood_ch19_data[["CTCF"]]
blood_ch19_mef2a <- blood_ch19_data[["MEF2A"]]
blood_ch19_rxra <- blood_ch19_data[["RXRA"]]

#STEM CELL
file_paths <- list(
  "SOX2" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy1704-[IMPACT_369_predictions_SOX2_chr19.bedgraph.gz]",
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy1418-[IMPACT_339_predictions_CTCF_chr19.bedgraph.gz]",
  "RXRA" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy1044-[IMPACT_618_predictions_RXRA_chr19.bedgraph.gz]"
)
sc_ch19_data <- lapply(file_paths, process_bedgraph, chrom = "chr19")
sc_ch19_sox2 <- sc_ch19_data[["SOX2"]]
sc_ch19_ctcf <- sc_ch19_data[["CTCF"]]
sc_ch19_rxra <- sc_ch19_data[["RXRA"]]


#EYE
file_paths <- list(
  "CTCF" = "/Users/zeynepkorkmaz/Desktop/UCSD/Statistical_Genomics/project/chr19/Galaxy32-[IMPACT_327_predictions_CTCF_chr19.bedgraph.gz]"
)
eye_ch19_data <- lapply(file_paths, process_bedgraph, chrom = "chr19")
eye_ch19_ctcf <- eye_ch19_data[["CTCF"]]


chrom19_bounds <- chrom_bounds %>% filter(chrom == 19)
lowest_chrom19 <- chrom19_bounds$min_pos
highest_chrom19 <- chrom19_bounds$max_pos

# Data for different tissues
bedgraph_list_ctcf <- list(
  brain = brain_ch19_ctcf,
  stem_cell = sc_ch19_ctcf,
  blood = blood_ch19_ctcf,
  eye = eye_ch19_ctcf
)

bedgraph_list_mef2a <- list(
  brain = brain_ch19_mef2a,
  blood = blood_ch19_mef2a
)

bedgraph_list_sox2 <- list(
  brain = brain_ch19_sox2,
  stem_cell = sc_ch19_sox2
)

bedgraph_list_rxra <- list(
  brain = brain_ch19_rxra,
  blood = blood_ch19_rxra,
  stem_cell = sc_ch19_rxra
)


chrom <- "chr19"
# Plot for CTCF
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_ctcf,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom19,
  chromend = highest_chrom19,
  color_map = color_map,
  title = "CTCF Transcription Factor"
)

# Plot for MEF2A
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_mef2a,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom19,
  chromend = highest_chrom19,
  color_map = color_map,
  title = "MEF2A Transcription Factor"
)

# Plot for SOX2
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_sox2,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom19,
  chromend = highest_chrom19,
  color_map = color_map,
  title = "SOX2 Transcription Factor"
)

#Plot for RXRA
plot_bedgraph_by_tissue(
  bedgraph_list = bedgraph_list_rxra,
  tissue_labels = tissue_labels,
  chrom = chrom,
  chromstart = lowest_chrom19,
  chromend = highest_chrom19,
  color_map = color_map,
  title = "RXRA Transcription Factor"
)
