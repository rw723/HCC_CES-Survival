# Packages
library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(ggtext)
library(scales)
library(tidyverse)
library(ggrepel)
library(eulerr)
library(grid)

# Output folders and saving helpers
fig_dir <- "figures_plos"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plos_dpi <- 600
plos_max_width_in <- 7.5
plos_max_height_in <- 8.75

save_plos_ggplot <- function(plot, filename_base, width_in = 7.5, height_in = 6,
                             dpi = plos_dpi, cap_height = TRUE) {
  width <- min(width_in, plos_max_width_in)
  height <- if (cap_height) min(height_in, plos_max_height_in) else height_in
  
  ggsave(
    filename = file.path(fig_dir, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )
}

# Load CESAnalysis objects and MAF files
clca_maf <- preload_maf(
  maf = "./hcc_clca_hg38.maf",
  refset = "ces.refset.hg38",
  keep_extra_columns = TRUE
)
clca_cesa <- load_cesa("clca_cesa.rds")

tcga_maf <- preload_maf(
  maf = "./lihc_tcga_hg38.maf",
  refset = "ces.refset.hg38",
  keep_extra_columns = TRUE
)
tcga_cesa <- load_cesa("tcga_cesa.rds")

# Build variant-level data.tables for visualization
CLCA_variant_CES <- as.data.table(copy(clca_cesa$selection$CLCA_variant_CES))
CLCA_variant_CES[, variant_name := sub("\\s*\\(.*$", "", variant_name)]
CLCA_variant_CES[, gene := fifelse(
  is.na(gene) | gene == "",
  clca_maf$Hugo_Symbol[match(variant_id, clca_maf$variant_id)],
  gene
)]

TCGA_variant_CES <- as.data.table(copy(tcga_cesa$selection$TCGA_variant_CES))
TCGA_variant_CES[, variant_name := sub("\\s*\\(.*$", "", variant_name)]
TCGA_variant_CES[, gene := fifelse(
  is.na(gene) | gene == "",
  tcga_maf$Hugo_Symbol[match(variant_id, tcga_maf$variant_id)],
  gene
)]

# Candidates revealed by the CLCA study
genes_interest <- c(
  "TP53", "CTNNB1", "ALB", "AXIN1", "ARID1A", "RB1", "TSC2", "ARID2", "JAK1", "KEAP1",
  "BRD7", "FGA", "TSC1", "ACVR2A", "PTEN", "RPS6KA3", "HNF1A", "PRDM11", "CDKN2A",
  "CDKN1B", "BMP5", "RPL22", "ECHS1", "TERT", "ZNF595", "KCNJ12", "KHNYN",
  "OR2A7", "NEAT1", "G035338", "Z95704.4", "RMRP", "G085970", "RN7SK", "G032906",
  "RNU12", "RP11-1151B14.3", "ADH1B", "PPP1R12B", "SEC14L2", "SERPINA1", "ADH4",
  "RABGEF1", "KCTD6", "PPP1R10", "HIST1H4C", "POLR2A", "SERBP1", "HIST1H1E"
)

# Plot the CES for HCC-CLCA: gene-level grouped plot
# Identify genes with at least 3 variants.
gene_keep <- CLCA_variant_CES[
  , .(tot_variants = sum(included_with_variant, na.rm = TRUE)),
  by = gene
][tot_variants >= 3, gene]

CLCA_keep <- CLCA_variant_CES[gene %in% gene_keep]

# Set variant-type labels for a manuscript-style legend.
# In cancereffectsizeR output, amino-acid-changing variants are usually coded as "aac"
# and non-coding/single-nucleotide variants as "snv".
variant_type_colors <- c(
  "aac" = "#4C72B0",
  "snv" = "#e7d159"
)

variant_type_labels <- c(
  "aac" = "Coding",
  "snv" = "Non-coding"
)

p1 <- plot_effects(
  effects = CLCA_keep,
  topn = NULL,
  group_by = "gene",
  label_individual_variants = FALSE,
  order_by_effect = TRUE,
  prevalence_method = "both",
  color_by = "variant_type",
  legend_size_name = "Variant prevalence\n(percent of samples)"
) +
  scale_fill_manual(
    values = variant_type_colors,
    breaks = names(variant_type_colors),
    labels = variant_type_labels,
    name = "Variant type"
  ) +
  scale_x_log10() +
  labs(
    x = "Cancer effect size\n(scaled selection coefficient)",
    y = "Gene"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )

p1_gene_levels <- p1$scales$get_scales("y")$get_labels()
p1_to_style <- ifelse(p1_gene_levels %in% genes_interest, "italic", "bold.italic")
p1 <- p1 + theme(axis.text.y = element_text(face = p1_to_style))

p1
save_plos_ggplot(
  plot = p1,
  filename_base = "Fig_gene_CLCA_CES_p1",
  width_in = 7.5,
  height_in = max(6, 0.25 * length(p1_gene_levels)),
  cap_height = TRUE
)

#  Prepare CLCA/TCGA AAC data for cross-cohort panels
# Identify TCGA genes with at least 3 variants
gene_keep_tcga <- TCGA_variant_CES[
  , .(tot_variants = sum(included_with_variant, na.rm = TRUE)),
  by = gene
][tot_variants >= 3, gene]

TCGA_keep <- TCGA_variant_CES[gene %in% gene_keep_tcga]

CLCA_aac <- copy(as.data.table(CLCA_keep[variant_type == "aac"]))
TCGA_aac <- copy(as.data.table(TCGA_keep[variant_type == "aac"]))

CLCA_aac[, cohort := "CLCA"]
TCGA_aac[, cohort := "TCGA"]

CLCA_aac[, variant_key := as.character(variant_id)]
TCGA_aac[, variant_key := as.character(variant_id)]

CLCA_aac[, cohort := as.character(cohort)]
TCGA_aac[, cohort := as.character(cohort)]
CLCA_aac[, gene := as.character(gene)]
TCGA_aac[, gene := as.character(gene)]
CLCA_aac[, variant_name := as.character(variant_name)]
TCGA_aac[, variant_name := as.character(variant_name)]

# Merge AAC tables for overlap and rank analysis.
AAC_merged <- rbindlist(
  list(CLCA_aac, TCGA_aac),
  use.names = TRUE,
  fill = TRUE
)
setnames(AAC_merged, make.unique(names(AAC_merged)))

AAC_merged <- AAC_merged[
  !is.na(selection_intensity) &
    !is.na(variant_key) &
    variant_key != ""
]
AAC_merged[, cohort := as.character(cohort)]
AAC_merged[, variant_key := as.character(variant_key)]

# Create clean labels: TP53 R249S, not TP53 TP53 R249S.
AAC_merged[, clean_variant_name := as.character(variant_name)]
AAC_merged[, clean_variant_name := sub("\\s*\\(.*$", "", clean_variant_name)]

AAC_merged[
  !is.na(gene) & gene != "" &
    !is.na(clean_variant_name) & clean_variant_name != "",
  clean_variant_name := mapply(
    function(g, v) {
      v2 <- sub(pattern = paste0("^\\s*", g, "\\s+"), replacement = "", x = v)
      v2 <- sub(pattern = paste0("^\\s*", g, "\\s+"), replacement = "", x = v2)
      paste(g, v2)
    },
    gene,
    clean_variant_name,
    USE.NAMES = FALSE
  )
]

AAC_merged[
  is.na(clean_variant_name) | clean_variant_name == "",
  clean_variant_name := variant_key
]

# Create exactly one label per variant_key.
variant_label_map <- AAC_merged[
  ,
  {
    clca_label <- clean_variant_name[
      cohort == "CLCA" & !is.na(clean_variant_name) & clean_variant_name != ""
    ]
    tcga_label <- clean_variant_name[
      cohort == "TCGA" & !is.na(clean_variant_name) & clean_variant_name != ""
    ]
    chosen_label <- if (length(clca_label) > 0) {
      clca_label[1]
    } else if (length(tcga_label) > 0) {
      tcga_label[1]
    } else {
      variant_key[1]
    }
    .(variant_label_fixed = chosen_label)
  },
  by = variant_key
]

AAC_merged <- merge(AAC_merged, variant_label_map, by = "variant_key", all.x = TRUE)
AAC_merged[is.na(variant_label_fixed) | variant_label_fixed == "", variant_label_fixed := variant_key]

# Safety check: should be empty.
label_problem_check <- AAC_merged[
  , .(n_labels = uniqueN(variant_label_fixed), labels = paste(unique(variant_label_fixed), collapse = " | ")),
  by = variant_key
][n_labels > 1]
print(label_problem_check)

# Create shared / cohort-specific status.
variant_status <- AAC_merged[
  , .(in_CLCA = any(cohort == "CLCA"), in_TCGA = any(cohort == "TCGA")),
  by = variant_key
][
  , status := fifelse(
    in_CLCA & in_TCGA,
    "Shared",
    fifelse(in_CLCA, "CLCA-specific", "TCGA-specific")
  )
]

if ("status" %in% names(AAC_merged)) AAC_merged[, status := NULL]
AAC_merged <- merge(AAC_merged, variant_status[, .(variant_key, status)], by = "variant_key", all.x = TRUE)
AAC_merged[, status := factor(status, levels = c("Shared", "CLCA-specific", "TCGA-specific"))]

overlap_counts <- unique(AAC_merged[, .(variant_key, status)])[, .N, by = status]
print(overlap_counts)


# Overlap Venn diagrams
save_two_set_venn <- function(
    set1,
    set2,
    set1_name,
    set2_name,
    filename_base,
    main_title,
    fill_colors = c("#4C72B0", "darkseagreen4"),
    width_in = 3,
    height_in = 3,
    dpi = plos_dpi
) {
  set1 <- unique(na.omit(as.character(set1)))
  set2 <- unique(na.omit(as.character(set2)))
  set1 <- set1[set1 != ""]
  set2 <- set2[set2 != ""]
  
  n_set1_only <- length(setdiff(set1, set2))
  n_set2_only <- length(setdiff(set2, set1))
  n_shared <- length(intersect(set1, set2))
  
  cat("\n", main_title, "\n", sep = "")
  cat(set1_name, " total: ", length(set1), "\n", sep = "")
  cat(set2_name, " total: ", length(set2), "\n", sep = "")
  cat(set1_name, "-specific: ", n_set1_only, "\n", sep = "")
  cat(set2_name, "-specific: ", n_set2_only, "\n", sep = "")
  cat("Shared: ", n_shared, "\n", sep = "")
  
  if ((n_set1_only + n_set2_only + n_shared) == 0) {
    stop(paste0("Venn input is empty for ", filename_base))
  }
  
  venn_counts <- stats::setNames(
    c(n_set1_only, n_set2_only, n_shared),
    c(set1_name, set2_name, paste0(set1_name, "&", set2_name))
  )
  
  venn_fit <- eulerr::euler(venn_counts)
  
  draw_venn <- function() {
    grid::grid.newpage()
    venn_plot <- plot(
      venn_fit,
      fills = list(
        fill = fill_colors,
        alpha = 0.55
      ),
      edges = list(
        col = fill_colors,
        lwd = 1.2
      ),
      labels = list(
        font = 2,
        cex = 1.2
      ),
      quantities = list(
        type = "counts",
        cex = 1.2
      ),
      main = main_title
    )
    
    # Draw the returned object only if it is drawable.
    if (inherits(venn_plot, "grob") || inherits(venn_plot, "gTree") || inherits(venn_plot, "gList")) {
      grid::grid.draw(venn_plot)
    }
  }
  
  png(
    filename = file.path(fig_dir, paste0(filename_base, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi
  )
  draw_venn()
  dev.off()
  
  invisible(list(
    fit = venn_fit,
    set1_total = length(set1),
    set2_total = length(set2),
    set1_only = n_set1_only,
    set2_only = n_set2_only,
    shared = n_shared
  ))
}

# AAC variant-level overlap.
venn_variant_summary <- save_two_set_venn(
  set1 = CLCA_aac$variant_key,
  set2 = TCGA_aac$variant_key,
  set1_name = "CLCA",
  set2_name = "TCGA",
  filename_base = "Panel_A1_AAC_variant_level_venn",
  main_title = "Variant overlap"
)

# AAC gene-level overlap.
venn_gene_summary <- save_two_set_venn(
  set1 = CLCA_aac$gene,
  set2 = TCGA_aac$gene,
  set1_name = "CLCA",
  set2_name = "TCGA",
  filename_base = "Panel_A2_AAC_gene_level_venn",
  main_title = "Gene overlap"
)

clca_variant_keys <- unique(as.character(CLCA_aac$variant_key))
tcga_variant_keys <- unique(as.character(TCGA_aac$variant_key))

clca_specific_keys <- setdiff(clca_variant_keys, tcga_variant_keys)
tcga_specific_keys <- setdiff(tcga_variant_keys, clca_variant_keys)

clca_specific_labels <- unique(AAC_merged[
  variant_key %in% clca_specific_keys,
  variant_label_fixed
])

tcga_specific_labels <- unique(AAC_merged[
  variant_key %in% tcga_specific_keys,
  variant_label_fixed
])


# Panel B: CLCA AAC plot_effects()
CLCA_plot <- copy(CLCA_aac)

pB <- plot_effects(
  effects = CLCA_plot,
  topn = nrow(CLCA_plot),
  label_individual_variants = TRUE,
  prevalence_method = "both",
  color_by = "#4C72B0",
  order_by_effect = TRUE,
  show_ci = TRUE,
  legend_size_name = "Variant prevalence\n(percent of samples)"
) +
  scale_x_log10() +
  labs(
    x = "Cancer effect size\n(scaled selection coefficient)",
    y = "Variants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.y = element_text(face = "italic", size = 10),
    legend.position = "right"
  )

pB_y_levels <- pB$scales$get_scales("y")$get_labels()
pB_to_style <- ifelse(pB_y_levels %in% clca_specific_labels, "bold.italic", "italic")
pB <- pB + theme(axis.text.y = element_text(face = pB_to_style, size = 10))

pB
save_plos_ggplot(
  plot = pB,
  filename_base = "Panel_B_CLCA_AAC_plot_effects",
  width_in = 7.5,
  height_in = 10,
  cap_height = FALSE
)

# Panel C: TCGA AAC plot_effects()
TCGA_plot <- copy(TCGA_aac)

pC <- plot_effects(
  effects = TCGA_plot,
  topn =  nrow(TCGA_plot),
  label_individual_variants = TRUE,
  prevalence_method = "both",
  color_by = "darkseagreen4",
  order_by_effect = TRUE,
  show_ci = TRUE,
  legend_size_name = "Variant prevalence\n(percent of samples)"
) +
  scale_x_log10() +
  labs(
    x = "Cancer effect size\n(scaled selection coefficient)",
    y = "Variants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.y = element_text(face = "italic", size = 10),
    legend.position = "right"
  )

pC_y_levels <- pC$scales$get_scales("y")$get_labels()
pC_to_style <- ifelse(pC_y_levels %in% tcga_specific_labels, "bold.italic", "italic")
pC <- pC + theme(axis.text.y = element_text(face = pC_to_style, size = 10))

pC
save_plos_ggplot(
  plot = pC,
  filename_base = "Panel_C_TCGA_AAC_plot_effects",
  width_in = 7.5,
  height_in = 7.2,
  cap_height = FALSE
)

# Rank-rank scatter for shared AAC variants
rank_dt <- as.data.table(copy(AAC_merged))
setnames(rank_dt, make.unique(names(rank_dt)))

rank_dt[, variant_key := as.character(variant_key)]
rank_dt[, cohort := as.character(cohort)]
rank_dt[, variant_label_fixed := as.character(variant_label_fixed)]

rank_dt <- rank_dt[
  !is.na(variant_key) & variant_key != "" &
    cohort %in% c("CLCA", "TCGA") &
    !is.na(selection_intensity) &
    !is.na(variant_label_fixed) & variant_label_fixed != ""
]

# Keep one row per variant per cohort.
rank_dt <- rank_dt[
  order(variant_key, cohort, -selection_intensity)
][
  , .SD[1],
  by = .(variant_key, cohort)
]

shared_keys <- rank_dt[, .(n_cohort = uniqueN(cohort)), by = variant_key][n_cohort == 2, variant_key]
shared_dt <- copy(rank_dt[variant_key %in% shared_keys])
cat("Number of shared AAC variants:", uniqueN(shared_dt$variant_key), "\n")

# Rank shared variants within each cohort. Rank 1 = strongest selection.
shared_dt[, shared_rank := frank(-selection_intensity, ties.method = "first"), by = cohort]

rank_wide <- dcast(
  shared_dt,
  variant_key + variant_label_fixed ~ cohort,
  value.var = "shared_rank"
)
rank_wide <- rank_wide[!is.na(CLCA) & !is.na(TCGA)]
rank_wide[, abs_rank_shift := abs(TCGA - CLCA)]

rank_wide[, rank_direction := fifelse(
  TCGA > CLCA,
  "Higher rank in CLCA",
  fifelse(TCGA < CLCA, "Higher rank in TCGA", "Same rank")
)]
rank_wide[, rank_direction := factor(
  rank_direction,
  levels = c("Higher rank in CLCA", "Higher rank in TCGA", "Same rank")
)]

# Label every shared variant.
rank_wide[, label_plot := variant_label_fixed]
missing_label_check <- rank_wide[
  is.na(label_plot) | label_plot == "",
  .(variant_key, variant_label_fixed, CLCA, TCGA, rank_direction)
]
print(missing_label_check)

rank_colors <- c(
  "Higher rank in CLCA" = "#4C72B0",
  "Higher rank in TCGA" = "darkseagreen4",
  "Same rank" = "black"
)

max_rank <- max(c(rank_wide$CLCA, rank_wide$TCGA), na.rm = TRUE)
rank_axis_max <- ceiling(max_rank / 5) * 5
rank_breaks <- seq(0, rank_axis_max, by = 5)

pD <- ggplot(rank_wide, aes(x = CLCA, y = TCGA)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.6
  ) +
  geom_point(
    aes(fill = rank_direction),
    shape = 21,
    size = 3.2,
    color = "black",
    stroke = 0.3,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    aes(label = label_plot, color = rank_direction),
    size = 3.1,
    fontface = "italic",
    na.rm = FALSE,
    max.overlaps = Inf,
    max.time = 10,
    max.iter = 200000,
    force = 1.5,
    force_pull = 0.10,
    box.padding = 0.45,
    point.padding = 0.30,
    min.segment.length = 0,
    segment.color = "gray65",
    segment.size = 0.25,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = rank_colors, name = "Relative rank") +
  scale_color_manual(values = rank_colors, guide = "none") +
  scale_x_reverse(
    limits = c(rank_axis_max, 0),
    breaks = rank_breaks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_reverse(
    limits = c(rank_axis_max, 0),
    breaks = rank_breaks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    x = "CLCA within-cohort rank",
    y = "TCGA within-cohort rank"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 15),
    axis.text = element_text(color = "black"),
    plot.margin = margin(8, 8, 8, 8)
  )

save_plos_ggplot(
  plot = pD,
  filename_base = "Panel_D_Shared_AAC_rank_rank_CLCA_TCGA",
  width_in = 10,
  height_in = 6.8,
  cap_height = FALSE
)


# Epistasis
CLCA_gene_epi <- as.data.table(copy(clca_cesa$epistasis$CLCA_gene_epi))
CLCA_gene_epi[, p_epistasis_adj := p.adjust(p_epistasis, method = "BH")]

CLCA_significant_epi <- CLCA_gene_epi[p_epistasis <= 0.05]

CLCA_significant_epi[, `:=`(
  pair_id = .I,
  A_change_dir = ifelse(ces_A_on_B >= ces_A0, "rising", "lowering"),
  B_change_dir = ifelse(ces_B_on_A >= ces_B0, "rising", "lowering"),
  p_label = paste0(
    "Nominal P = ", scales::pvalue(p_epistasis, accuracy = 0.001),
    "\nBH-adjusted P = ", scales::pvalue(p_epistasis_adj, accuracy = 0.001)
  )
)]

plot_data_A <- CLCA_significant_epi[, .(
  pair_id,
  p_label,
  gene = variant_A,
  change_dir = A_change_dir,
  selection_independent = ces_A0,
  selection_epistatic = ces_A_on_B,
  p_epistasis,
  p_epistasis_adj
)]

plot_data_B <- CLCA_significant_epi[, .(
  pair_id,
  p_label,
  gene = variant_B,
  change_dir = B_change_dir,
  selection_independent = ces_B0,
  selection_epistatic = ces_B_on_A,
  p_epistasis,
  p_epistasis_adj
)]

plot_data <- rbind(plot_data_A, plot_data_B)

plot_data[, gene_order := seq_len(.N), by = pair_id]
plot_data[, x_offset := ifelse(gene_order == 1, -0.08, 0.08), by = pair_id]
plot_data[, x_ind := 0.5 + x_offset]
plot_data[, x_epi := 2.5 + x_offset]
plot_data[, x_lab := 2.0 + x_offset]

p4 <- ggplot(plot_data) +
  geom_segment(
    aes(
      x = x_ind,
      xend = x_epi,
      y = selection_independent,
      yend = selection_epistatic,
      color = change_dir
    ),
    linewidth = 1.2,
    arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed")
  ) +
  geom_point(
    aes(
      x = x_ind,
      y = selection_independent,
      color = change_dir
    ),
    size = 4
  ) +
  geom_text_repel(
    aes(
      x = x_lab,
      y = selection_epistatic,
      label = gene
    ),
    data = . %>% filter(gene != "ENSG00000291313"),
    color = "black",
    nudge_x = 0.7,
    direction = "y",
    hjust = 0,
    size = 4,
    fontface = "italic",
    segment.color = NA
  ) +
  geom_text_repel(
    aes(
      x = x_lab,
      y = selection_epistatic,
      label = gene
    ),
    data = . %>% filter(gene == "ENSG00000291313"),
    color = "black",
    nudge_x = 0.7,
    nudge_y = 0.5,
    direction = "y",
    hjust = 0,
    size = 4,
    segment.color = NA
  ) +
  facet_wrap(
    ~p_label,
    ncol = 3,
    scales = "free_x"
  ) +
  scale_y_log10(
    labels = scales::label_scientific(digits = 2),
    n.breaks = 8
  ) +
  scale_color_manual(
    values = c(
      "rising" = "#e41a1c",
      "lowering" = "#377eb8"
    )
  ) +
  scale_x_continuous(
    name = "Condition",
    breaks = c(0.5, 2.5),
    labels = c("Independent", "Epistatic"),
    limits = c(0, 4)
  ) +
  labs(
    y = "Estimated selection coefficient (log scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, hjust = 0),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 11),
    panel.spacing = unit(1.5, "lines")
  )

p4


# TCGA epistasis
TCGA_gene_epi <- as.data.table(copy(tcga_cesa$epistasis$TCGA_gene_epi))
TCGA_gene_epi[, p_epistasis_adj := p.adjust(p_epistasis, method = "BH")]

TCGA_significant_epi <- TCGA_gene_epi[p_epistasis <= 0.05]

cat(
  "Number of nominally significant TCGA epistasis pairs:",
  nrow(TCGA_significant_epi),
  "\n"
)

if (nrow(TCGA_significant_epi) == 0) {
  p5 <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "TCGA: no gene-level epistasis pairs with nominal P <= 0.05",
      size = 5
    ) +
    theme_void()
  
  plot_data_TCGA <- data.table(pair_id = integer())
  
} else {
  TCGA_significant_epi[, `:=`(
    pair_id = .I,
    A_change_dir = ifelse(ces_A_on_B >= ces_A0, "rising", "lowering"),
    B_change_dir = ifelse(ces_B_on_A >= ces_B0, "rising", "lowering"),
    p_label = paste0(
      "Nominal P = ", scales::pvalue(p_epistasis, accuracy = 0.001),
      "\nBH-adjusted P = ", scales::pvalue(p_epistasis_adj, accuracy = 0.001)
    )
  )]
  
  plot_data_A_TCGA <- TCGA_significant_epi[, .(
    pair_id,
    p_label,
    gene = variant_A,
    change_dir = A_change_dir,
    selection_independent = ces_A0,
    selection_epistatic = ces_A_on_B,
    p_epistasis,
    p_epistasis_adj
  )]
  
  plot_data_B_TCGA <- TCGA_significant_epi[, .(
    pair_id,
    p_label,
    gene = variant_B,
    change_dir = B_change_dir,
    selection_independent = ces_B0,
    selection_epistatic = ces_B_on_A,
    p_epistasis,
    p_epistasis_adj
  )]
  
  plot_data_TCGA <- rbind(plot_data_A_TCGA, plot_data_B_TCGA)
  
  plot_data_TCGA <- plot_data_TCGA[
    !is.na(selection_independent) &
      !is.na(selection_epistatic) &
      selection_independent > 0 &
      selection_epistatic > 0
  ]
  
  plot_data_TCGA[, gene_order := seq_len(.N), by = pair_id]
  plot_data_TCGA[, x_offset := ifelse(gene_order == 1, -0.08, 0.08), by = pair_id]
  plot_data_TCGA[, x_ind := 0.5 + x_offset]
  plot_data_TCGA[, x_epi := 2.5 + x_offset]
  plot_data_TCGA[, x_lab := 2.0 + x_offset]
  
  p5 <- ggplot(plot_data_TCGA) +
    geom_segment(
      aes(
        x = x_ind,
        xend = x_epi,
        y = selection_independent,
        yend = selection_epistatic,
        color = change_dir
      ),
      linewidth = 1.2,
      arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed")
    ) +
    geom_point(
      aes(
        x = x_ind,
        y = selection_independent,
        color = change_dir
      ),
      size = 4
    ) +
    geom_text_repel(
      aes(
        x = x_lab,
        y = selection_epistatic,
        label = gene
      ),
      color = "black",
      nudge_x = 0.7,
      direction = "y",
      hjust = 0,
      size = 4,
      fontface = "italic",
      segment.color = NA,
      max.overlaps = Inf
    ) +
    facet_wrap(
      ~p_label,
      ncol = 3,
      scales = "free_x"
    ) +
    scale_y_log10(
      labels = scales::label_scientific(digits = 2),
      n.breaks = 8
    ) +
    scale_color_manual(
      values = c(
        "rising" = "#e41a1c",
        "lowering" = "#377eb8"
      )
    ) +
    scale_x_continuous(
      name = "Condition",
      breaks = c(0.5, 2.5),
      labels = c("Independent", "Epistatic"),
      limits = c(0, 4)
    ) +
    labs(
      y = "Estimated selection coefficient (log scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 11, hjust = 0),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 11),
      panel.spacing = unit(1.5, "lines")
    )
}

p5


save_plos_ggplot(
  plot = p4,
  filename_base = "Panel_E_CLCA_gene_epistasis_p4",
  width_in = 12,
  cap_height = FALSE
)

save_plos_ggplot(
  plot = p5,
  filename_base = "Panel_F_TCGA_gene_epistasis_p5",
  width_in = 12,
  cap_height = FALSE
)

