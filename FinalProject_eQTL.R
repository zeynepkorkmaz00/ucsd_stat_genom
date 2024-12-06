library(data.table)
library(plink2R) 
library(dplyr)
library(ggplot2)

setwd("~/OneDrive - UC San Diego/DSC291/ucsd_stat_genom")

# load GWAS sumstats for MS
gwas <- fread("imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz") # 155756 SNPs

# load expression data  
gene_exp <- fread("GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", header = T)

# find significant SNPs
gwas_sig_8 <- gwas[gwas$p < 5e-8,] # 3782
gwas_sig_6 <- gwas[gwas$p < 5e-6,] # 5137

# construct Manhattan plot to visualize significant SNPs
# change p values of 0 to something very small (ie 1e-350)
gwas$p[gwas$p == 0] <- 1e-350
gwas$log10_p_values <- -log(x = gwas$p, base = 10) # add a column with -log10 p-value
gwas$color <- ifelse(gwas$chrom %% 2 == 0, "blue", "black") # add a column with the color to use

x_labs <- paste("Chr",unique(gwas$chrom)) # labels for x axis
# calculate position of label by counting # SNPs in each chrom, finding cumsum, and adding midpoint
snp_per_chrom <- gwas %>% group_by(chrom) %>% summarise(count = n()) 
pos <- c(0,cumsum(snp_per_chrom$count)[1:21]) + 1 + snp_per_chrom$count/2

# plotting Manhattan plot
p1 <- ggplot(gwas, aes(x = 1:nrow(gwas), y = log10_p_values, color = color)) +
  geom_point() + 
  geom_hline(yintercept = -log(x = 5e-8, base = 10), color = "red", linetype = "solid") +  # Red horizontal line
  scale_color_identity() +
  labs(x = "Chromosome", y = "-log10(P-value)", title = "Manhattan Plot for Multiple Sclerosis GWAS") + 
  scale_x_continuous(breaks = pos, labels = x_labs) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file = "ManhattanPlot.pdf", width = 4*4, height = 2*4)
p1
dev.off()

# many SNPs are in MHC region on chrom 6: can separate out the SNPs in the MHC region vs those not 
# taken from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
MHC_min <- 28477797
MHC_max <- 33448354

gwas_sig_MHC <- gwas_sig_8[(gwas_sig_8$chrom == 6) & (gwas_sig_8$pos >= MHC_min) & (gwas_sig_8$pos <= MHC_max),] # 2642
gwas_sig_nonMHC <- gwas_sig_8[!(gwas_sig_8$rsid %in% gwas_sig_MHC$rsid),] # 1140

# can find the genes that are within 500 kB of each of these SNPs before doing eQTL (takes too long to go through every potential gene) 
genes_list <- c()
for (ind in 1:nrow(gwas_sig_8)) {
  row <- gwas_sig_8[ind,]
  start <- row$pos - 500000
  end <- row$pos + 500000
  genes <- gene_exp[(gene_exp$Chr == as.character(row$chrom)) & (gene_exp$Coord >= start) & (gene_exp$Coord <= end),]
  genes_list <- c(genes_list, genes$Gene_Symbol)
}
genes_list <- unique(genes_list) # 1100 genes to be tested in eQTL
gene_exp_trunc <- gene_exp[gene_exp$Gene_Symbol %in% genes_list,]

#### Perform eQTL similarly to how we did it in HW3 in order to associate gene expression with SNPs #### 
# All GWAS is from individual of European ancestry so using that for the genotyping data too 
# not splitting into train/test patients because we just want to use eQTL to identify SNP-gene relationships
# for every gene, use the set of significant SNPs from the GWAS to run the eQTL analysis. only interested in the potential association of significant SNPs with genes

# load in just 1 chrom to identify intersect of patients in genotypes and gene exp
genotypes <- read_plink("~/OneDrive - UC San Diego/DSC291/1000G_Eur/1000G_eur_chr22")
patientID <-  genotypes$fam$V1[genotypes$fam$V1 %in% colnames(gene_exp)] # 358 patients in both genotypes and gene exp
n <- length(patientID)

geno_mat <- matrix(,nrow = n) # initialize so you can cbind
for (i in 1:22) {
  print(i)
  genotypes <- read_plink(paste0("~/OneDrive - UC San Diego/DSC291/1000G_Eur/1000G_eur_chr",i))
  
  snp_ind <- which(colnames(genotypes$bed) %in% gwas_sig_8$rsid) # only save genotype info for SNPs in gwas_sig_8
  patients <- unname(sapply(rownames(genotypes$bed), function(x) strsplit(x, ":")[[1]][1])) # remove colon from patient name in bed
  patient_ind <- which(patients %in% patientID) # indices of 358 patients
  snp_geno_matrix <- genotypes$bed[patient_ind,snp_ind] # extract genotype info for these SNPs and patients
  
  # add SNPs from this chrom onto overall geno_mat
  geno_mat <- cbind(geno_mat,snp_geno_matrix)
}
geno_mat <- geno_mat[,-1] # remove first empty column from geno_mat
rownames(geno_mat) <- patientID

dim(geno_mat) # 358 3588

geno_mat_std <- scale(geno_mat, scale = T) 
any(is.na(geno_mat_std)) # FALSE

save(file = "geno_mat_std_FINAL.RData",geno_mat_std)

# 1 - PCA on the genotypes
# dont need to downsample because there are only 3588 SNPs
# transpose and perform PCA 
pca_genotype <- prcomp(t(geno_mat_std), scale = T)

# 2 - PCA on the gene expression
# reorder gene expression matrix in order of genotyped individuals
gene_exp_mat <- as.data.frame(gene_exp_trunc)[,rownames(geno_mat_std)]
rownames(gene_exp_mat) <- gene_exp_trunc$Gene_Symbol
gene_exp_mat <- t(gene_exp_mat) # so patients are rows like in genotype mat
dim(gene_exp_mat) # 358 1100

# standardise: each gene has mean 0 and var 1
gene_exp_mat_std <- scale(gene_exp_mat, scale = T) 
any(is.na(gene_exp_mat_std)) # FALSE

# transpose and perform PCA
pca_gene_exp <- prcomp(t(gene_exp_mat_std), scale = T)

save(file = "pca_FINAL.RData",pca_gene_exp,pca_genotype)

# 3 - Regress out covariates
# combine 10 genotyping PCs and 10 gene expression PCs into one covariate matrix
cov_mat <- cbind(pca_genotype$rotation[,1:10], pca_gene_exp$rotation[,1:10])
dim(cov_mat) # 358  20

# for each gene
eqtl_sumstats_all <- data.frame(SNP=c(),Beta=c(),SE=c(),P=c(),Gene=c())
j = 1
for (gene in colnames(gene_exp_mat_std)) {
  # extract gene
  print(j)
  gene_vec_std <- gene_exp_mat_std[,grep(gene,colnames(gene_exp_mat_std))]
  lm_model <- lm(gene_vec_std ~ cov_mat)
  ge_vector_std_resid <- resid(lm_model)
  
  # perform eQTL analysis across the 3588 SNPs
  eqtl_sumstats <- matrix(0,ncol(geno_mat_std),4)
  eqtl_sumstats[,1] <- colnames(geno_mat_std)
  for (i in 1:ncol(geno_mat_std)){ 
    if(!all(is.na(geno_mat_std[,i]))){ # all NA should be removed but keeping just in case
      mod <- lm(ge_vector_std_resid ~ geno_mat_std[,i]) 
      eqtl_sumstats[i,2:4] <- summary(mod)$coefficients[2,c(1,2,4)]
    }else{
      eqtl_sumstats[i,2:4] <- NA
    }
  }
  colnames(eqtl_sumstats) <- c("SNP","Beta","SE","P")
  eqtl_sumstats <- eqtl_sumstats %>% as.data.frame()
  eqtl_sumstats$P <- as.numeric(eqtl_sumstats$P)
  eqtl_sumstats$Beta <- as.numeric(eqtl_sumstats$Beta)
  eqtl_sumstats$SE <- as.numeric(eqtl_sumstats$SE)
  eqtl_sumstats$Gene <- rep(gene, ncol(geno_mat_std))
  
  eqtl_sumstats_all <- rbind(eqtl_sumstats_all,eqtl_sumstats)
  j = j+1
}

save(file = "eqtl_sumstats_all.RData", eqtl_sumstats_all)

# use ENCODE to get ChIPseq data to identify regulatory regions and overlap those with our SNPs
# take H3K4me1 (active and poised enhancers) from brain organoid
# https://www.encodeproject.org/experiments/ENCSR844OQI/
chipseq <- fread("ENCFF026GBE.bed.gz")

# figure out which of our SNPs fall into a peak region in the chipseq data
enhancer_snps <- c()
enhancer_snps_df <- data.frame(SNP = c(), chr = c(), enhancer_start = c(), enhancer_end = c())
for (ind in 1:nrow(gwas_sig_8)) {
  row <- gwas_sig_8[ind,]
  chipseq_chrom <- chipseq[chipseq$V1 == paste0("chr",row$chrom)]
  a <- chipseq_chrom[(chipseq_chrom$V2 <= row$pos) & (chipseq_chrom$V3 >= row$pos),]
  if (dim(a)[1] != 0) { # if there are rows, ie if the SNP is in a peak
    enhancer_snps <- c(enhancer_snps, row$rsid)
    enhancer_snps_df <- rbind(enhancer_snps_df,data.frame(SNP = row$rsid, chr = a$V1, enhancer_start = a$V2, enhancer_end = a$V3))
  }
}
length(enhancer_snps) # 365
dim(enhancer_snps_df) # 365 4

# intersect significant eQTL associations with the enhancer SNP list
thresh <- 0.05/ncol(geno_mat_std) # divide by #SNPs tested for each gene
eqtl_sumstats_sig <- eqtl_sumstats_all[eqtl_sumstats_all$P < thresh,] # 1029

eqtl_enhancer <- eqtl_sumstats_sig[eqtl_sumstats_sig$SNP %in% intersect(eqtl_sumstats_sig$SNP,enhancer_snps),]
dim(eqtl_enhancer) # 74 5 
# 74 SNPs in enhancer regions that are significantly correlated with genes
length(unique(eqtl_enhancer$Gene)) # 33 unique genes

# map ENSEMBL id to interpretable gene name
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# remove version number from genes
ensembl_ids <- sub("\\.\\d+$", "", unique(eqtl_enhancer$Gene))
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = mart)

eqtl_enhancer$Gene <- sub("\\.\\d+$", "", eqtl_enhancer$Gene)
eqtl_enhancer <- merge(eqtl_enhancer, gene_mapping, by.x = "Gene", by.y = "ensembl_gene_id")

# gene_mapping maps each of the 33 unique ENSMBL ids to its gene name 
# eqtl_enhancer stores each of the 74 SNPs in enhancer regions that are significantly correlated with genes
# enhancer_snps_df has the chr, start and end of the enhancer that each of the 365 SNPs is in
save(file = "enhancer_gene_mapping.RData", gene_mapping, eqtl_enhancer, enhancer_snps_df)





