#generate deseq object and VST counts
library(tidyr)
library(dplyr)
library(DESeq2)
library(tximport)

directory <- "/path"

norm_rnaseq <- function(species, directory) {
  # Retrieve annotation file
  print("Retrieving annotation file..")
  gtf_file <- rtracklayer::import(paste0(species, "_annotation.gtf"))
  gtf_df <- as.data.frame(gtf_file)
  list_files <- list.files(directory)
  
  # Read metadata
  metadata <- read.csv(file = paste0(species, "_metadata.csv"))
  samples <- c(metadata$FileName)
  files <- file.path(directory, samples, "abundance.h5")
  names(files) <- paste(c(metadata$SampleName))
  
  # Set directory
  dirs <- file.path(directory, samples)
  datatable <- data.frame(sample = samples, condition = metadata$condition, path = dirs)
  
  # Import and process data
  print("Importing counts file using tximport..")
  txi <- tximport(files, type = "kallisto", tx2gene = gtf_df[, c("transcript_id", "gene_id")], txOut = FALSE)
  txi_counts <- as.data.frame(txi$counts)
  txi_counts <- txi_counts %>% mutate_if(is.numeric, round)
  txi_ab <- as.data.frame(txi$abundance)
  
  txi_counts$gene_id <- rownames(txi_counts)
  row.names(txi_counts) <- NULL
  txi_counts <- txi_counts[, c("gene_id", names(txi_counts)[-ncol(txi_counts)])]
  
  txi_ab$gene_id <- rownames(txi_ab)
  row.names(txi_ab) <- NULL
  txi_ab<- txi_ab[, c("gene_id", names(txi_ab)[-ncol(txi_ab)])]
  
  new_name <- paste0(species, "_counts")
  assign(new_name, txi_counts)
  
  write.table(txi_counts, file = paste0(species, "_counts.tsv"), sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(txi_ab, file = paste0(species,"_abundance.tsv"), sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Create DESeqDataSetFromTximport object
  dds <- DESeqDataSetFromTximport(txi, metadata, ~condition)
  dds$condition <- relevel(dds$condition, "prm") #use relevel to set your reference condition
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  # Normalize data
  print("Normalizing counts file using DESeq..")
  norm <- DESeq(dds)

  plotfile <- paste0(species, "_dispersion_plot.png")
  png(plotfile, width = 800, height = 400)
  plotDispEsts(norm)
  dev.off()
  
  return(list(norm = norm, counts = txi_counts, abundance = txi_ab, metadata = metadata, txi = txi, dds = dds))
}

vst_rnaseq <- function(norm, species, orthofile, prefix = "") {
  #transform counts using VST
  print("Transforming counts using Variance Stabilizing Transformation..")
  vst <- varianceStabilizingTransformation(norm, blind=F)
  vst_df <- as.data.frame(assay(vst))
  vst_df$gene_id <- paste0(prefix, rownames(vst_df))
  row.names(vst_df) <- NULL
  vst_df <- vst_df[, c("gene_id", names(vst_df)[-ncol(vst_df)])]
  vst_mean <- vst_df
  vst_mean$mean_value <- rowMeans(vst_df[, -c(1)], na.rm = TRUE)
  vst_mean$cluster_id <- orthofile$cluster_id[match(vst_mean$gene_id, orthofile$gene_id)]
  vst_mean <- na.omit(vst_mean)
  
  vst_ortho <- vst_mean %>%
    group_by(cluster_id) %>%
    slice(which.max(mean_value))
  vst_orthodb <- vst_ortho  %>%
    mutate(taxid = species)
  vst_orthodb <- vst_orthodb[, c("gene_id", "cluster_id", "taxid", "mean_value")]
  vst_orthodb$prot_id <- orthofile$prot_id[match(vst_orthodb$gene_id, orthofile$gene_id)]
  
  new_name1 <- paste0(species, "_vst_df")
  assign(new_name1, vst_df)
  new_name2 <- paste0(species, "_vst_mean")
  assign(new_name2, vst_mean)
  new_name3 <- paste0(species, "_vst_ortho")
  assign(new_name3, vst_ortho)
  new_name4 <- paste0(species, "_vst_orthodb")
  assign(new_name4, vst_orthodb)
  
  write.table(vst_mean, file = paste0(species, "_vst_mean.tsv"), sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(vst_df, file = paste0(species, "_vst.tsv"), sep = '\t', quote = FALSE, row.names = FALSE)
  return(list(vst_df = vst_df, vst_mean = vst_mean, vst_ortho = vst_ortho,vst_orthodb = vst_orthodb))
}

dmel_deseq <- norm_rnaseq(species, directory)
zcu_deseq <- norm_rnaseq(species, directory)
bacdorsa_instar1_deseq <- norm_rnaseq(species, directory)
bacdorsa_instar2_deseq <- norm_rnaseq(species, directory)
bacdorsa_deseq2 <- norm_rnaseq(species, directory)
lvannamei_deseq <- norm_rnaseq(species, directory)
scypa_deseq <- norm_rnaseq(species, directory)
bacdorsa_all_deseq2 <- norm_rnaseq(species, directory)

#generate vst counts
species_taxids <- c("6689", "7227", "27457", "28588", "85552", "72036")
kvpairtable <- read.csv(file = "kvpairtable.csv", header = TRUE, check.names = FALSE, sep = ',')
orthologer_030325 <- read.csv("orthology_file.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
orthologer_030325 <- orthologer_030325[, 1:3]
orthologer_030325$gene_id <- kvpairtable$gene_id[match(orthologer_030325$prot_id, kvpairtable$protein_id)]
orthologer_030325 <- na.omit(orthologer_030325)

dmel_vst <- vst_rnaseq(dmel_deseq$norm, "dmel", orthologer_030325)
zcu_vst <- vst_rnaseq(zcu_deseq$norm, "zcu", orthologer_030325, "zcu_")
bacdorsa_all_vst <- vst_rnaseq(bacdorsa_all_deseq$norm, "bacdorsa", orthologer_030325, "bd_")
lvannamei_vst <- vst_rnaseq(lvannamei_deseq$norm, "lvannamei", orthologer_030325,"lv_")
scypa_vst <- vst_rnaseq(scypa_deseq$norm, "scypa", orthologer_030325,"sc_")
all_orthodb <- rbind(dmel_vst$vst_orthodb, zcu_vst$vst_orthodb, bacdorsa_all_vst$vst_orthodb, lvannamei_vst$vst_orthodb, scypa_vst$vst_orthodb)
all_ortho_singlecopy <- all_orthodb %>% group_by(cluster_id) %>% filter(n_distinct(taxid) == 5)

#generate orthology matrix

bacdorsa_counts <- as.data.frame(bacdorsa_all_deseq$counts)
bacdorsa_counts <- bacdorsa_counts[, new_col_order]
bacdorsa_counts$gene_id <- paste0("bd_", bacdorsa_counts$gene_id)
bacdorsa_counts_orthogroup <- bacdorsa_counts %>% filter(bacdorsa_counts$gene_id %in% all_ortho_singlecopy$gene_id)
bacdorsa_counts_orthogroup$cluster_id <- all_ortho_singlecopy$cluster_id[match(bacdorsa_counts_orthogroup$gene_id, all_ortho_singlecopy$gene_id)]
rownames(bacdorsa_counts_orthogroup) <- bacdorsa_counts_orthogroup$gene_id

zcu_counts <- as.data.frame(zcu_deseq$counts)
zcu_counts$gene_id <- paste0("zcu_", zcu_counts$gene_id)
zcu_counts_orthogroup <- zcu_counts %>% filter(zcu_counts$gene_id %in% all_ortho_singlecopy$gene_id)
zcu_counts_orthogroup$cluster_id <- all_ortho_singlecopy$cluster_id[match(zcu_counts_orthogroup$gene_id, all_ortho_singlecopy$gene_id)]
rownames(zcu_counts_orthogroup) <- zcu_counts_orthogroup$gene_id

dmel_counts <- as.data.frame(dmel_deseq$counts)
dmel_counts_orthogroup <- dmel_counts %>% filter(dmel_counts$gene_id %in% all_ortho_singlecopy$gene_id)
dmel_counts_orthogroup$cluster_id <- all_ortho_singlecopy$cluster_id[match(dmel_counts_orthogroup$gene_id, all_ortho_singlecopy$gene_id)]
rownames(dmel_counts_orthogroup)<- dmel_counts_orthogroup$gene_id 

lvannamei_counts <- as.data.frame(lvannamei_deseq$counts)
lvannamei_counts$gene_id <- paste0("lv_", lvannamei_counts$gene_id)
lvannamei_counts_orthogroup <- lvannamei_counts %>% filter(lvannamei_counts$gene_id %in% all_ortho_singlecopy$gene_id)
lvannamei_counts_orthogroup$cluster_id <- all_ortho_singlecopy$cluster_id[match(lvannamei_counts_orthogroup$gene_id, all_ortho_singlecopy$gene_id)]
rownames(lvannamei_counts_orthogroup)<- lvannamei_counts_orthogroup$gene_id 

scypa_counts <- as.data.frame(scypa_deseq$counts)
scypa_counts$gene_id <- paste0("sc_", scypa_counts$gene_id)
scypa_counts_orthogroup <- scypa_counts %>% filter(scypa_counts$gene_id %in% all_ortho_singlecopy$gene_id)
scypa_counts_orthogroup$cluster_id <- all_ortho_singlecopy$cluster_id[match(scypa_counts_orthogroup$gene_id, all_ortho_singlecopy$gene_id)]
rownames(scypa_counts_orthogroup)<- scypa_counts_orthogroup$gene_id 

merge1 <- merge(dmel_counts_orthogroup, zcu_counts_orthogroup, by = "cluster_id")
merge2 <- merge(bacdorsa_counts_orthogroup,lvannamei_counts_orthogroup, by = "cluster_id")
merge3 <- merge(merge1,merge2, by = "cluster_id")
all_counts_5species <- merge(merge3,scypa_counts_orthogroup, by = "cluster_id")

write.table(bacdorsa_counts, file = "bacdorsa_counts.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(dmel_counts, file = "dmel_counts.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(zcu_counts, file = "zcu_counts.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(lvannamei_counts, file = "lvannamei_counts.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(scypa_counts, file = "scypa_counts.tsv", sep = '\t', quote = FALSE, row.names = FALSE)

all_metadata <- read.delim("all_metadata.txt", header = TRUE, sep = "\t")
rownames(all_metadata) <- all_metadata$name

all_counts_deseq <- DESeqDataSetFromMatrix(
  countData = all_counts_5species,
  colData = all_metadata,  
  design = ~ condition + species
)

all_species_norm <- DESeq(all_counts_deseq)

# Save the dispersion plot
plotfile <- "all_dispersion.png"
png(plotfile, width = 8, height = 4, units = "in", res = 400)
plotDispEsts(all_species_norm)
dev.off()

all_counts_deseq$condition <- relevel(all_counts_deseq$condition, "prm")
all_species_vst <- varianceStabilizingTransformation(all_species_norm, blind=F)
all_species_vst <- as.data.frame(assay(all_species_vst))

library(PCAtools)
p <- pca(all_species_vst, metadata = all_metadata, removeVar = 0.05)
pca_scores <- p$rotated
pca_scores$sample <- rownames(pca_scores)
pca_scores <- merge(pca_scores, all_metadata, by.x = "sample", by.y = "name")

percent_variance <- round(p$variance, 2)
pca_scores$groupings <- as.factor(pca_scores$groupings)

plotfile <- "PCA_PC1PC2.png"
png(plotfile, width = 6.4, height = 7, units = "in", res = 400)
ggplot(pca_scores, aes(x = PC1, y = PC2, color = stage, shape = species)) +
  geom_point(size = 3, alpha = 0.9) +  
  stat_ellipse(
    aes(group = groupings),   
    type = "norm",          
    level = 0.8,           
    alpha = 0.3,           
    color = "black",       
    linewidth = 0.5         
  ) +
  labs(
    x = paste0("PC1 (", percent_variance[1], "%)"),
    y = paste0("PC2 (", percent_variance[2], "%)")
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add black border
  )
dev.off()










