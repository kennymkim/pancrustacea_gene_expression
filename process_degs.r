#process DEGs

library(tidyr)
library(rtracklayer)
library(dplyr)
library(DESeq2)

dge_rnaseq <- function(norm, species) {
  print("Generating LFC results table..")
  res.wald <- results(norm)
  res.wald <- res.wald[order(res.wald$padj),]
  subset <- subset(res.wald, padj<0.05)
  summary(res.wald, alpha=0.05)
  summary_output <- capture.output(summary(res.wald, alpha = 0.05))
  sum(subset$padj < 0.05, na.rm=TRUE)
  resultsNames(norm)
  
  print("Running LFC shrinkage..")
  lfcshrink <- lfcShrink(dds=norm, type = "ashr", coef =2)
  lfcshrink <- lfcshrink[order(lfcshrink$padj),]
  
  print("Extracting significant differentially expressed genes...")
  DEGs <- lfcshrink[!is.na(lfcshrink$padj) & lfcshrink$padj < 0.05, ]
  upreg <- DEGs[DEGs$log2FoldChange > 0, ]
  downreg <- DEGs[DEGs$log2FoldChange < 0, ]
  deglist <- rownames(DEGs)
  
  write.table(as.data.frame(lfcshrink),file = paste0(species,"_lfcshrink.tsv"), quote = FALSE, sep = '\t')
  write.table(as.data.frame(res.wald),file = paste0(species,"_lfcwald.tsv"), quote = FALSE, sep = '\t')
  write.table(as.data.frame(DEGs),file = paste0(species,"_DEGs.tsv"), quote = FALSE,sep = '\t')
  write.table(as.data.frame(upreg),file = paste0(species,"_upDEGs.tsv"), quote = FALSE, sep = '\t')
  write.table(as.data.frame(downreg),file = paste0(species,"_downDEGs.tsv"), quote = FALSE,sep = '\t')
  writeLines(summary_output, con = paste0(species, "_wald_summary.txt"))
  
  return(list(lfcshrink = lfcshrink, deg_df = DEGs, deglist = deglist, upreg = upreg, downreg = downreg, wald = res.wald, wald_subset = subset ))
}

#get DEGs using DESeq
dmel_deg <- dge_rnaseq(dmel_deseq$norm, "dmel")
zcu_deg <- dge_rnaseq(zcu_deseq$norm, "zcu")
bacdorsa_deg <- dge_rnaseq(bacdorsa_deseq$norm, "bacdorsa")
lvannamei_deg <- dge_rnaseq(lvannamei_deseq$norm, "lvannamei")
scypa_deg <- dge_rnaseq(scypa_deseq$norm, "scypa")

bacdorsa_instar1_deg <- dge_rnaseq(bacdorsa_instar1_deseq$norm, "bacdorsa_instar1")
bacdorsa_instar2_deg <- dge_rnaseq(bacdorsa_instar2_deseq$norm,"bacdorsa_instar2")
bacdorsa_instar_deg_df <- rbind(bacdorsa_instar1_df, bacdorsa_instar2_df)
bacdorsa_instar_deg <- union(bacdorsa_instar1_deg$deglist, bacdorsa_instar2_deg$deglist)
bacdorsa_instar_deg_union <- bacdorsa_instar_deg_df %>% filter(rownames(bacdorsa_instar_deg_df )%in% bacdorsa_instar_deg)

get_coding_genes <- function(gtf_file) {
  gtf_data <- import(gtf_file)
  gtf_df <- as.data.frame(gtf_data)
  
  if (!"gene_biotype" %in% colnames(gtf_df)) {
    stop("Error: The GTF file does not contain a 'gene_biotype' column.")
  }
  coding_genes <- gtf_df %>%
    filter(type == "gene", gene_biotype == "protein_coding") %>%
    pull(gene_id) %>%
    unique()
  
  return(coding_genes)
}

dmel_coding_genes <- get_coding_genes("dmel_pupa_annotation.gtf")
bacdorsa_coding_genes <- get_coding_genes("bacdorsa_annotation.gtf")
zcu_coding_genes <- get_coding_genes("zcu_annotation.gtf")
lvannamei_coding_genes <- get_coding_genes("lvannamei_annotation.gtf")
scypa_coding_genes <- get_coding_genes("scypa_annotation.gtf")

filter_coding_DEGs <- function(deg_df, coding_genes, output_prefix = "", taxid = "") {
  filtered_deg <- deg_df %>%
    filter(rownames(deg_df) %in% coding_genes) %>%
    dplyr::select(log2FoldChange) %>%
    tibble::rownames_to_column("gene_id")

  if (output_prefix != "") {
    filtered_deg$gene_id <- paste0(output_prefix, filtered_deg$gene_id)
  }
  filtered_deg$taxid <- taxid
  
  return(filtered_deg)
}

dmel_coding_degs <- filter_coding_DEGs(as.data.frame(dmel_deg$deg_df), dmel_coding_genes, taxid = "dmel")
zcu_coding_degs <- filter_coding_DEGs(as.data.frame(zcu_deg$deg_df), zcu_coding_genes, "zcu_", taxid = "zcu")
bacdorsa_coding_degs <- filter_coding_DEGs(as.data.frame(bacdorsa_deg$deg_df), bacdorsa_coding_genes,"bd_", taxid = "bacdorsa")
bacdorsa_instar_coding_degs <- filter_coding_DEGs(bacdorsa_instar_deg_union, bacdorsa_coding_genes, "bd_", taxid = "bacdorsa_instar" )
lvannamei_coding_degs <- filter_coding_DEGs(as.data.frame(lvannamei_deg$deg_df), lvannamei_coding_genes,"lv_", taxid = "lvannamei")
scypa_coding_degs <- filter_coding_DEGs(as.data.frame(scypa_deg$deg_df), scypa_coding_genes, "sc_", taxid = "scypa")
