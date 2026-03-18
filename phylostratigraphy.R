#construct phyloset df

pt_lvannamei_gene_ages <- read.csv("6689_gene_ages.tsv", sep="\t", header=TRUE)
pt_lvannamei_phyloset <-  left_join(lvannamei_vst, pt_lvannamei_gene_ages[, c("gene_id","rank")], by = "gene_id")
pt_lvannamei_phyloset <- pt_lvannamei_phyloset %>% select(rank, everything())
pt_lvannamei_phyloset <- pt_lvannamei_phyloset %>% distinct()
pt_lvannamei_phyloset <- na.omit(pt_lvannamei_phyloset)
pt_lvannamei_phyloset$Phylostratum <- pt_lvannamei_phyloset$rank
pt_lvannamei_phyloset$rank <- NULL
pt_lvannamei_phyloset <- pt_lvannamei_phyloset %>% select(Phylostratum, everything())

pt_dmel_gene_ages <- read.csv("7227_gene_ages.tsv", sep="\t", header=TRUE)
pt_dmel_phyloset <-  left_join(dmel_vst , pt_dmel_gene_ages[, c("gene_id","rank")], by = "gene_id")
pt_dmel_phyloset <- pt_dmel_phyloset %>% select(rank, everything())
pt_dmel_phyloset <- pt_dmel_phyloset %>% distinct()
pt_dmel_phyloset <- na.omit(pt_dmel_phyloset)
pt_dmel_phyloset$Phylostratum <- pt_dmel_phyloset$rank
pt_dmel_phyloset$rank <- NULL
pt_dmel_phyloset <- pt_dmel_phyloset %>% select(Phylostratum, everything()) 

# Recode the Phylostratum column
pt_lvannamei_phyloset <- pt_lvannamei_phyloset %>%
  mutate(
    Phylostratum = case_when(
      Phylostratum >= 1 & Phylostratum <= 8 ~ 1,
      Phylostratum == 9 ~ 2,
      Phylostratum == 10 ~ 3,
      Phylostratum >= 11 & Phylostratum <= 14 ~ 4,
      Phylostratum >= 15 & Phylostratum <= 17 ~ 5,
      TRUE ~ Phylostratum # Preserve original if not in specified range
    )
  )

pt_dmel_phyloset <- pt_dmel_phyloset %>%
  mutate(
    Phylostratum = case_when(
      Phylostratum >= 1 & Phylostratum <= 8 ~ 1,
      Phylostratum == 9 ~ 2,
      Phylostratum == 10 ~ 3,
      Phylostratum >= 11 & Phylostratum <= 14 ~ 4,
      Phylostratum >= 15 & Phylostratum <= 26 ~ 5,
      TRUE ~ Phylostratum # Preserve original if not in specified range
    )
  )

pt_lvannamei_TAI <- PlotSignature(ExpressionSet = pt_lvannamei_phyloset,
              measure       = "TAI", 
              TestStatistic = "FlatLineTest",
              xlab          = "Ontogeny", 
              ylab          = "TAI" )

pt_dmel_TAI <- PlotSignature(ExpressionSet = pt_dmel_phyloset,
              measure       = "TAI", 
              TestStatistic = "FlatLineTest",
              xlab          = "Ontogeny", 
              ylab          = "TAI" )


ggplot(pt_dmel_TAI$data, aes(x = factor(Stage, levels = c("L3","wandering_early", "wandering_mid","wandering_late",  "p0", "p1", "p2")), y = TI, group = 1)) +
  geom_line(size = 2, color = "blue") +
  theme_minimal() +
  labs(
    title = "Transcriptome Age Index (TAI)",
    x = "Stage",
    y = "TAI"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12)
  )
