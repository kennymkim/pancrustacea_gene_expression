#generate deseq object and VST counts

directory <- "/path"
species <- "dmel"

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






