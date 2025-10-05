# Load the ces and reference genome package
library(cancereffectsizeR)
library(ces.refset.hg38)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(annotatr)


# Create cancereffectsizeR analysis object
#HCC-CLCA
clca_maf <- preload_maf(maf = "./hcc_clca_hg38.maf", refset = "ces.refset.hg38", 
                        keep_extra_columns = T)
clca_cesa <- CESAnalysis(refset = "ces.refset.hg38")
clca_cesa <- load_maf(cesa = clca_cesa, maf =  clca_maf, coverage = "genome") #coverage definition


#LIHC-TCGA
tcga_maf <- preload_maf(maf = "./lihc_tcga_hg38.maf",refset  = "ces.refset.hg38", 
                        keep_extra_columns = T)
tcga_cesa <- CESAnalysis(refset = "ces.refset.hg38")
tcga_cesa <- load_maf(cesa = tcga_cesa, maf  = tcga_maf, coverage = "exome")

# Mutational signature extraction with suggested exclusion
# Sample are treatment-free
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LIHC", treatment_naive = TRUE) 

# Inferred relative rates of mutation, produced by matrix-multiplying biological_weights and signature definitions.
#HCC-CLCA
clca_cesa <- trinuc_mutation_rates(
  cesa = clca_cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)
#LIHC-TCGA
tcga_cesa <- trinuc_mutation_rates(
  cesa = tcga_cesa, signature_set = ces.refset.hg38$signatures$COSMIC_v3.4,
  signature_exclusions = signature_exclusions
)

# Estimate neutral gene mutation rates using dNdScv, with tissue-specific mutation rate covariates.
#HCC-CLCA
clca_cesa <- gene_mutation_rates(clca_cesa, 
                                 dndscv_args = list(max_muts_per_gene_per_sample = Inf),
                                 covariates = ces.refset.hg38$covariates$LIHC)

#LIHC-TCGA
tcga_cesa <- gene_mutation_rates(tcga_cesa,
                                 dndscv_args = list(max_muts_per_gene_per_sample = Inf),
                                 covariates = ces.refset.hg38$covariates$LIHC)

### CES inference
## Basic model:
# HCC-CLCA
clca_cesa <- ces_variant(cesa = clca_cesa, run_name = "CLCA_variant_CES")
# LIHC-TCGA
tcga_cesa <- ces_variant(cesa = tcga_cesa, run_name = "TCGA_variant_CES")

## Epistatic model:
# Build data.table for visualization
CLCA_variant_CES <- clca_cesa$selection$CLCA_variant_CES
clca_gene_keep <- CLCA_variant_CES[variant_type == "aac",.(tot_variants = sum(included_with_variant, na.rm = TRUE)),
                               by = gene][ tot_variants >= 3, gene ]
TCGA_variant_CES <- tcga_cesa$selection$TCGA_variant_CES
tcga_gene_keep <- TCGA_variant_CES[variant_type == "aac",.(tot_variants = sum(included_with_variant, na.rm = TRUE)),
                               by = gene][ tot_variants >= 3, gene ]
# gene-level epistasis for genes with at least 3 variants
clca_cesa <- ces_gene_epistasis(
  clca_cesa,
  genes      = clca_gene_keep,
  run_name   = "CLCA_gene_epi")

tcga_cesa <- ces_gene_epistasis(
  tcga_cesa,
  genes      = clca_gene_keep, # use the same genes as in CLCA
  run_name   = "TCGA_gene_epi")
TCGA_gene_epi = tcga_cesa$epistasis$TCGA_gene_epi
# Save the CESAnalysis objects
save_cesa(clca_cesa, 'clca_cesa.rds')
save_cesa(tcga_cesa, 'tcga_cesa.rds')

