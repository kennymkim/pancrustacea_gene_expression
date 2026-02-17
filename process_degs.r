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

#all DEGs
all_degdf <- rbind(dmel_coding_degs,bacdorsa_coding_degs, bacdorsa_instar_coding_degs, zcu_coding_degs,lvannamei_coding_degs,scypa_coding_degs)
all_degdf$cluster_id <- all_orthodb$cluster_id[match(all_degdf$gene_id, all_orthodb$gene_id)]
all_degdf <- all_degdf %>% mutate(cluster_id = ifelse(is.na(cluster_id), gene_id, as.character(cluster_id)))

all_degdf_filtered_single_copy <- all_degdf %>% filter(all_degdf$cluster_id %in% all_ortho_singlecopy$cluster_id)
all_degdf <- all_degdf %>%
  mutate(lineage_taxid = case_when(
    taxid %in% c("bacdorsa", "bacdorsa_instar") ~ "bacdorsa",
    TRUE ~ taxid
  ))

#classify the DEGs 
crustacea <- c("lvannamei", "scypa")
insecta <- c("dmel", "zcu", "bacdorsa", "bacdorsa_instar")
metamorphic <- c("dmel", "zcu", "bacdorsa")
non_metamorphic <- c("lvannamei", "scypa", "bacdorsa_instar")

classify_and_summarize_clusters <- function(all_degdf, insecta, crustacea, custom_order) {
  
  cluster_classification <- all_degdf %>%
    group_by(cluster_id) %>%
    summarise(
      taxa_in_cluster = list(unique(taxid)),
      n_insecta = sum(taxid %in% insecta),
      n_crustacea = sum(taxid %in% crustacea),
      n_species = n_distinct(lineage_taxid),
      classification = case_when(
        n_species == 1 ~ "Species specific",
        n_crustacea >= 2 & n_insecta == 0 ~ "Malacostraca specific",
        n_insecta >= 2 & n_crustacea == 0 & !("bacdorsa_instar" %in% unlist(taxa_in_cluster)) ~ "Metamorphic",
        n_insecta >= 2 & n_crustacea == 0 & "bacdorsa_instar" %in% unlist(taxa_in_cluster) ~ "Insecta specific",
        n_crustacea >= 1 & n_insecta == 1 & "bacdorsa_instar" %in% unlist(taxa_in_cluster) ~ "Non metamorphic",
        n_crustacea >= 1 & n_insecta >= 1 ~ "Pancrustacea shared",
        TRUE ~ "Other"
      ),
      .groups = "drop"
    )
  
  all_classified <- all_degdf %>%
    left_join(cluster_classification %>% select(cluster_id, classification), by = "cluster_id")

  summary_df <- all_classified %>%
    group_by(taxid, classification) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    arrange(taxid, desc(n_genes))
  
  summary_long <- summary_df %>%
    rename(count = n_genes) %>%
    group_by(taxid) %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  summary_long$taxid <- factor(summary_long$taxid, levels = custom_order)
  summary_long$classification <- factor(
    summary_long$classification,
    levels = c(
      #"Orphan genes",
      "Species specific", 
      "Metamorphic", 
      "Insecta specific",
      "Malacostraca specific", 
      "Non metamorphic", 
      "Pancrustacea shared"
    )
  )
  
  # 7️⃣ Return both summary and classified tables
  return(list(summary = summary_long, classified = all_classified))
}

all_degdf_classified_all <- classify_and_summarize_clusters(
  all_degdf = all_degdf,
  insecta = insecta,
  crustacea = crustacea,
  custom_order = custom_order
)

all_degdf_metamorphic <- all_degdf_classified_df_all %>% filter(classification == "Metamorphic") %>% distinct(cluster_id)
all_degdf_nonmetamorphic <- all_degdf_classified_df_all %>% filter(classification == "Non metamorphic") %>% distinct(cluster_id)
all_degdf_classified_all_up <- all_degdf_classified_df_all %>% filter(log2FoldChange > 0)
all_degdf_classified_all_down <- all_degdf_classified_df_all %>% filter(log2FoldChange < 0)

all_degdf_classified_df_all  <- all_degdf_classified_all$classified
all_degdf_classified_summary_all <- all_degdf_classified_all$summary

all_degdf_classified_up_summary <- all_degdf_classified_all_up %>%
  group_by(taxid, classification) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = 100 * count / sum(count))

all_degdf_classified_down_summary <- all_degdf_classified_all_down %>%
  group_by(taxid, classification) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = 100 * count / sum(count))

insects <- c("dmel", "bacdorsa", "zcu")
malacostraca <- c("lvannamei", "scypa")

orthology_lookup <- all_orthodb %>%
  mutate(cluster_id = as.character(cluster_id)) %>%
  group_by(cluster_id) %>%
  summarise(
    has_insect = any(taxid %in% insects),
    has_malacostraca = any(taxid %in% malacostraca),
    lineage_status = case_when(
      has_insect & has_malacostraca ~ "Pancrustacea shared",
      has_insect & !has_malacostraca ~ "Insecta only",
      !has_insect & has_malacostraca ~ "Malacostraca only",
      TRUE ~ "No orthology / orphan"
    ),
    .groups = "drop"
  )

orthology_lookup <- orthology_lookup %>% mutate(cluster_id = as.character(cluster_id))
all_degdf_classified_df_all_annot <- all_degdf_classified_df_all %>% left_join(orthology_lookup, by = "cluster_id")
all_degdf_classified_df_all_annot <- all_degdf_classified_df_all_annot %>%
  mutate(
    lineage_status = ifelse(is.na(lineage_status),
                            "No orthology / orphan",
                            lineage_status)
  )

all_degdf_summary_orthology_level <- all_degdf_classified_df_all_annot %>%
  count(taxid, classification, lineage_status) %>%
  arrange(taxid, classification, lineage_status)

#plot the degdf summary
library(ggplot2)
library(patchwork)

plot_df <- all_degdf_summary_orthology_level %>%
  mutate(
    taxid = factor(taxid, levels = taxid_levels),
    classification = factor(classification, levels = class_levels),
    lineage_status = factor(lineage_status, levels = lineage_levels),
    x = paste(taxid, classification, sep = " | ")
  ) %>%
  arrange(taxid, classification) %>%
  mutate(x = factor(x, levels = rev(unique(x))))

right_species <- factor(c("lvannamei","dmel"), levels = c("lvannamei","dmel"))
middle_species <- factor(c("scypa","zcu"),levels = c("scypa","zcu"))
left_species <- factor(c("bacdorsa","bacdorsa_instar"),levels = c("bacdorsa","bacdorsa_instar"))

species_labels <- c(
  bacdorsa_instar = "B. dorsalis (instar)",
  bacdorsa        = "B. dorsalis",
  scypa            = "S. paramamosain",
  zcu              = "Z. cucurbitae",
  lvannamei        = "L. vannamei",
  dmel             = "D. melanogaster"
)

orthology_cols <- c(
  "All Pancrustacea"   = "#F8766D",
  "Insecta only"          = "#7CAE00",
  "Malacostraca only"     = "#00BFC4",
  "No detectable orthologs (orphan)" = "#C77CFF"
)
fill_scale <- scale_fill_manual(values = orthology_cols, drop = FALSE)

p_left <- plot_df %>%
  filter(taxid %in% left_species) %>%
  ggplot(aes(x = classification, y = n, fill = lineage_status)) +
  geom_col(color = "black", linewidth = 0.25) +
  facet_wrap(~ taxid, ncol = 1, scales = "free_y", labeller = labeller(taxid = species_labels)) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Orthology grouping"
  ) +
  fill_scale +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(color = "black", fill = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

p_middle <- plot_df %>%
  filter(taxid %in% middle_species) %>%
  ggplot(aes(x = classification, y = n, fill = lineage_status)) +
  geom_col(color = "black", linewidth = 0.25) +
  facet_wrap(~ taxid, ncol = 1, scales = "free_y", labeller = labeller(taxid = species_labels)) +
  coord_flip() +
  fill_scale +
  labs(
    x = NULL,
    y = "Number of DEGs",
    fill = "Orthology grouping"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(color = "black", fill = "grey90"),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_right <- plot_df %>%
  filter(taxid %in% right_species) %>%
  ggplot(aes(x = classification, y = n, fill = lineage_status)) +
  geom_col(color = "black", linewidth = 0.25) +
  facet_wrap(~ taxid, ncol = 1, scales = "free_y", labeller = labeller(taxid = species_labels)) +
  coord_flip() +
  fill_scale +
  labs(
    x = NULL,
    y = NULL,
    fill = "Orthology distribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(color = "black", fill = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

final_plot <- p_left|p_middle|p_right
final_plot

plotfile <- "Figure.png"
png(plotfile, width = 9, height = 5, units = "in", res = 600)
final_plot
dev.off()


