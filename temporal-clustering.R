process_data <- function(vst_df, cluster_genes, prefix = NULL) {
  
  #cluster_genes <- paste0(prefix, cluster_genes)
  vst_df <- vst_df[vst_df$gene_id %in% cluster_genes, ]
  vst_mat <- vst_df %>% select(-gene_id)
  cluster_zscores <- compute_zscores(vst_mat)
  
  zscores_df <- cbind(vst_df$gene_id, as.data.frame(cluster_zscores))
  colnames(zscores_df)[1] <- "gene_id"
  
  long_format <- zscores_df %>%
    pivot_longer(cols = -gene_id, names_to = "sample", values_to = "zscore")
  
  return(long_format)
}

pt_lvannamei_clust <- read.delim("lvannamei_clusters.tsv", header = TRUE, stringsAsFactors = FALSE)

tempo_lvannamei_C1_long <- process_data(pt_lvannamei_vst, pt_lvannamei_clust$C1) 
tempo_lvannamei_C2_long <- process_data(pt_lvannamei_vst, pt_lvannamei_clust$C2) 
tempo_lvannamei_C3_long <- process_data(pt_lvannamei_vst, pt_lvannamei_clust$C3) 
tempo_lvannamei_C4_long <- process_data(pt_lvannamei_vst, pt_lvannamei_clust$C4)
tempo_lvannamei_C5_long <- process_data(pt_lvannamei_vst, pt_lvannamei_clust$C5) 

tempo_lvannamei_stage_order <- c("early_premoult","mid_premoult","late_premoult","early_postmoult","late_postmoult","intermoult")

tempo_lvannamei_C1_long <- tempo_lvannamei_C1_long %>% mutate(cluster = "C1", stage = factor(sample, levels = tempo_lvannamei_stage_order))
tempo_lvannamei_C2_long <- tempo_lvannamei_C2_long %>% mutate(cluster = "C2", stage = factor(sample, levels = tempo_lvannamei_stage_order))
tempo_lvannamei_C3_long <- tempo_lvannamei_C3_long %>% mutate(cluster = "C3", stage = factor(sample, levels = tempo_lvannamei_stage_order))
tempo_lvannamei_C4_long <- tempo_lvannamei_C4_long %>% mutate(cluster = "C4", stage = factor(sample, levels = tempo_lvannamei_stage_order))
tempo_lvannamei_C5_long <- tempo_lvannamei_C5_long %>% mutate(cluster = "C5", stage = factor(sample, levels = tempo_lvannamei_stage_order))

cluster_order <- c("C1", "C2", "C3", "C4", "C5")
tempo_lvannamei_all <- bind_rows(
  tempo_lvannamei_C1_long,
  tempo_lvannamei_C2_long,
  tempo_lvannamei_C3_long,
  tempo_lvannamei_C4_long,
  tempo_lvannamei_C5_long
)

tempo_lvannamei_all <- tempo_lvannamei_all %>% mutate(cluster = factor(cluster, levels = cluster_order))
tempo_lvannamei_all$Phylostratum <- tempo_lvannamei_phyloset$Phylostratum[match(tempo_lvannamei_all$gene_id, tempo_lvannamei_phyloset$gene_id)]

phylostrata_pct <- tempo_lvannamei_all %>%
  filter(!is.na(Phylostratum)) %>%
  distinct(gene_id, Phylostratum, cluster) %>%
  count(cluster, Phylostratum) %>%
  group_by(cluster) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  ungroup()

phylostrata_pct <- phylostrata_pct %>% mutate(cluster = factor(cluster, levels = cluster_order))
phylostrata_pct <- phylostrata_pct %>%
  mutate(
    Phylostratum_label = case_when(
      Phylostratum == 1 ~ "Ancestral to Arthropoda",
      Phylostratum == 2 ~ "Arthropoda",
      Phylostratum == 3 ~ "Pancrustacea",
      Phylostratum == 4 ~ "Crustacea",
      Phylostratum == 5 ~ "Decapoda",
      TRUE ~ NA_character_
    )
  )
phylostrata_pct$Phylostratum <- phylostrata_pct$Phylostratum_label
phylostrata_pct$Phylostratum <- factor(phylostrata_pct$Phylostratum, levels = c("Ancestral to Arthropoda","Arthropoda", "Pancrustacea", "Crustacea", "Decapoda"))

lvannamei_cluster_stage_means <- tempo_lvannamei_all %>%
  group_by(cluster, stage) %>%
  summarise(mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop")

cluster_order <- c("C1","C2","C3","C4","C5")
tempo_lvannamei_all <- tempo_lvannamei_all %>% mutate(cluster = factor(cluster, levels = cluster_order))
phylostrata_pct <- phylostrata_pct %>% mutate(cluster = factor(cluster, levels = cluster_order))

p1 <- ggplot(tempo_lvannamei_all, aes(x = stage, y = zscore, group = gene_id)) +
  geom_line(color = "grey80", alpha = 0.5) +
  geom_line(
    data = lvannamei_cluster_stage_means,
    aes(x = stage, y = mean_zscore, group = 1),  # override group
    color = "red", size = 1.2
  ) +
  facet_wrap(~cluster, ncol = 1) +
  labs(
    title = "Gene Expression Cluster",
    x = "Stage",
    y = "Z-score"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

p2 <- ggplot(phylostrata_pct , aes(x = cluster, y = percentage, fill = as.factor(Phylostratum))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~cluster, ncol = 1, scales = "free_y") +  # one panel per cluster
  labs(
    title = "Phylostratum Composition",
    x = NULL,
    y = "Percentage",
    fill = "Phylostratum"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # adds border to each facet panel
    strip.background = element_blank(),
    strip.text = element_blank(),  # optionally remove facet labels since cluster is on x-axis
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12)
  )

combined_plot <- (p1 | p2) + plot_layout(widths = c(2.5, 1)) 
combined_plot 

plotfile <- "lvannamei_clusters.jpg"
png(plotfile, width = 7.5, height = 6.9, units = "in", res = 600)
combined_plot
dev.off()




