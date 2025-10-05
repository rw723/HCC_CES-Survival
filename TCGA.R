library(cancereffectsizeR)
library(data.table)
library(dplyr)

# Prepare data
maf <- preload_maf(maf = "cohortMAF.2025-04-23.maf.gz", refset = "ces.refset.hg38")

# Create cancereffectsizeR analysis and load data
cesaT <- CESAnalysis(refset = "ces.refset.hg38")
cesaT <- load_maf(cesa = cesaT, maf = maf)

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LIHC", treatment_naive = TRUE)
cesaT <- trinuc_mutation_rates(
  cesa = cesaT, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
cesaT <- gene_mutation_rates(cesaT, covariates = ces.refset.hg38$covariates$LIHC)

head(cesaT$gene_rates)
dndscv_results  = cesaT$dNdScv_results[[1]]
sig_genes = dndscv_results[qallsubs_cv <= .05][order(qallsubs_cv)]
sig_genes

cesaT <- ces_variant(cesaT, run_name = "example")

# Visualize top-effect variants.
select_gene = c("P4HA1","JAK1", "CTNNB1","ALB","TP53","FGA","HNF1A","PRDM11","CDKN1B","BMP5","ECHS1","AXIN1",
                "ARID1A","TLE1")
plot_effects(effects = cesaT$selection$example, #%>%
               #filter(gene  %in% select_gene),
             #group_by = "gene",
             label_individual_variants = F)

# Attribute effects to mutational signatures
mut_effects <- mutational_signature_effects(cesaT, cesaT$selection$example)

# Plot a comparison of how signatures contribute to mutation vs. selection
plot_signature_effects(mut_effects, viridis_option = "F")



