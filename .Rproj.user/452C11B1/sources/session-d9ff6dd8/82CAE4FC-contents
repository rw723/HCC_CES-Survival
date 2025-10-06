# Load the neccessary package
library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(ggtext)
library(scales)

# Load the CESAnalysis objects
clca_maf <- preload_maf(maf = "./hcc_clca_hg38.maf", refset = "ces.refset.hg38", 
                        keep_extra_columns = T)
clca_cesa <- load_cesa('clca_cesa.rds')

# Build data.table for visualization
CLCA_variant_CES <- clca_cesa$selection$CLCA_variant_CES
CLCA_variant_CES[ , variant_name := sub("\\s*\\(.*$", "", variant_name) ]
CLCA_variant_CES[, gene := fifelse(is.na(gene),clca_maf$Hugo_Symbol[match(variant_id, clca_maf$variant_id)], gene)]

# candidates revealed by the CLCA study
genes_interest <- c(
  "TP53","CTNNB1","ALB","AXIN1","ARID1A","RB1","TSC2","ARID2","JAK1","KEAP1",
  "BRD7","FGA","TSC1","ACVR2A","PTEN","RPS6KA3","HNF1A","PRDM11","CDKN2A",
  "CDKN1B","BMP5","RPL22","ECHS1","TERT","ZNF595","KCNJ12","KHNYN",
  "OR2A7","NEAT1","G035338","Z95704.4","RMRP","G085970","RN7SK","G032906",
  "RNU12","RP11-1151B14.3","ADH1B","PPP1R12B","SEC14L2","SERPINA1","ADH4",
  "RABGEF1","KCTD6","PPP1R10","HIST1H4C","POLR2A","SERBP1","HIST1H1E")

# Plot the CES for HCC-CLCA
# identify genes with at least 3 variants
gene_keep <- CLCA_variant_CES[ ,.(tot_variants = sum(included_with_variant, na.rm = TRUE)),
                               by = gene][ tot_variants >= 3, gene ] 
CLCA_keep  <- CLCA_variant_CES[ gene %in% gene_keep ]
# plot genes with at least 3 variants
p1 = plot_effects(
  effects   = CLCA_keep,
  topn      = NULL,
  group_by = "gene",  # group by gene
  label_individual_variants = T,
  order_by_effect = T,
  prevalence_method = "both",
  color_by = "variant_type"
) +
  scale_fill_manual(values = c("#4C72B0","#e7d159"), name = "Variant type") + 
  scale_x_log10()

p1_gene_levels <- p1$scales$get_scales("y")$get_labels()  # y axis order in the plot
p1_to_style    <- ifelse(p1_gene_levels %in% genes_interest, "italic", "bold.italic")
p1 <- p1 + theme(axis.text.y = element_text(face = p1_to_style))
p1
p2 = plot_effects(
  effects   = CLCA_keep[variant_type == "aac", ],
  topn      = 20,
  label_individual_variants = F,
  prevalence_method = "both",
  color_by  = "#4C72B0"
) + scale_x_log10() + theme(axis.text.y = element_text(face = "italic"))
p2
# TCGA variant CES
tcga_cesa <- load_cesa('tcga_cesa.rds')
tcga_maf <- preload_maf(maf = "./lihc_tcga_hg38.maf", refset = "ces.refset.hg38", 
                        keep_extra_columns = T)
TCGA_variant_CES <- tcga_cesa$selection$TCGA_variant_CES
TCGA_variant_CES[ , variant_name := sub("\\s*\\(.*$", "", variant_name) ]
TCGA_variant_CES[, gene := fifelse(is.na(gene), tcga_maf$Hugo_Symbol[match(variant_id, tcga_maf$variant_id)], gene)]
# Plot the CES for LIHC-TCGA
# plot genes with at least 3 variants
gene_keep_tcga <- TCGA_variant_CES[ ,
                               .(tot_variants = sum(included_with_variant, na.rm = TRUE)),
                               by = gene][ tot_variants >= 3, gene ] 
TCGA_keep  <- TCGA_variant_CES[ gene %in% gene_keep_tcga ]
p3 = plot_effects(
  effects = TCGA_keep[variant_type == "aac", ],
  topn = 20,
  label_individual_variants = F,
  prevalence_method = "both"
) + scale_x_log10() + theme(axis.text.y = element_text(face = "italic"))
p3

#Epistasis
CLCA_gene_epi = clca_cesa$epistasis$CLCA_gene_epi
CLCA_significant_epi <- CLCA_gene_epi[p_epistasis <= 0.05]
p4 = plot_epistasis(epistatic_effects = CLCA_significant_epi)

p4 = p4[[1]]+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "italic")
  )
p4


TCGA_gene_epi = tcga_cesa$epistasis$TCGA_gene_epi
significant_variants <- unique(
  unlist(TCGA_gene_epi[p_epistasis <= 0.05, .(variant_A, variant_B)])
)
heatmap_data <- TCGA_gene_epi[variant_A %in% significant_variants & 
                                variant_B %in% significant_variants]

p5 = ggplot(heatmap_data_upper, aes(x = variant_A, y = variant_B, fill = p_epistasis)) +
  geom_tile(color = "grey85") +
  # Use scale_fill_viridis_c() for a continuous viridis color scale
  scale_fill_viridis_c(
    name = "p-epistasis", 
    option = "D", # Option "D" is the default viridis palette
    limits = c(0, 0.05),
    oob = scales::squish, # Handles values outside the limits gracefully
    na.value = "grey50" # Color for any missing p-values on the diagonal
  ) +
  # Ensure the axes include all variants to maintain the square shape
  scale_x_discrete(limits = all_vars, drop = FALSE) +
  scale_y_discrete(limits = all_vars, drop = FALSE) +
  theme_minimal() +
  labs(
    x = "Gene A",
    y = "Gene B"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic")
  ) +
  coord_fixed()
p5

# Save the plots
ggsave(plot = p1, filename = 'Fig2_clca_all.png', width = 550 * 3, height = 800*2.5, units = 'px', dpi = 'retina')
ggsave(plot = p2, filename = 'Fig3_clca_aac.png', width = 500 * 3, height = 700*2, units = 'px', dpi = 'retina')
ggsave(plot = p3, filename = 'Fig3_tcga_aac.png', width = 500 * 3, height = 700*2, units = 'px', dpi = 'retina')
ggsave(plot = p4, filename = 'Fig4_clca_aac_epistasis.png', width = 800 * 3, height = 700*2, units = 'px', dpi = 'retina')
ggsave(plot = p5, filename = 'Fig4_tcga_aac_epistasis.png', width = 700 * 3, height = 700*3, units = 'px', dpi = 'retina')
