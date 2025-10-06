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
                        keep_extra_columns = c("Consequence", "Variant_Classification"))
clca_cesa <- CESAnalysis(refset = "ces.refset.hg38")
clca_cesa <- load_maf(cesa = clca_cesa, maf =  clca_maf, coverage = "genome") #coverage definition


#LIHC-TCGA
tcga_maf <- preload_maf(maf = "./lihc_tcga_hg38.maf",refset  = "ces.refset.hg38", 
                        keep_extra_columns = c("Consequence", "Variant_Classification"))
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

# CES inference
#HCC-CLCA
clca_cesa <- ces_variant(cesa = clca_cesa, run_name = "CLCA_variant_CES")
#LIHC-TCGA
tcga_cesa <- ces_variant(cesa = tcga_cesa, run_name = "TCGA_variant_CES")

# Save the CESAnalysis objects
save_cesa(clca_cesa, 'clca_cesa.rds')
save_cesa(tcga_cesa, 'tcga_cesa.rds')
# Load the CESAnalysis objects
clca_cesa <- load_cesa('clca_cesa.rds')
tcga_cesa <- load_cesa('tcga_cesa.rds')

# Build data.table for visualization
CLCA_variant_CES <- clca_cesa$selection$CLCA_variant_CES
TCGA_variant_CES <- tcga_cesa$selection$TCGA_variant_CES

CLCA_variant_CES$chr = clca_cesa$variants$chr[match(CLCA_variant_CES$variant_id, clca_cesa$variants$variant_id)]
CLCA_variant_CES$start = clca_cesa$variants$start[match(CLCA_variant_CES$variant_id, clca_cesa$variants$variant_id)]
CLCA_variant_CES$end = clca_cesa$variants$end[match(CLCA_variant_CES$variant_id, clca_cesa$variants$variant_id)]
## For HCC-CLCA
# Build a 1-bp GRanges of variant position
gr <- GRanges(
  seqnames   = CLCA_variant_CES$chr,
  ranges     = IRanges(start = CLCA_variant_CES$start,
                       end   = CLCA_variant_CES$end  ),
  strand     = "*",
  variant_id = CLCA_variant_CES$variant_id  # needed for the join
)
# annotation for regions
annots <- c(
  "hg38_basicgenes",
  "hg38_enhancers_fantom" ,
  "hg38_cpgs",
  "hg38_lncrna_gencode"
)
anno_gr <- build_annotations(genome = "hg38", annotations = annots)

# Use NCBI style for seqlevels
seqlevelsStyle(gr) <- "NCBI"
seqlevelsStyle(anno_gr) <- "NCBI"
ov <- annotate_regions(gr, anno_gr, ignore.strand = TRUE)
region_vec <- sub(":.*$", "", ov$annot$id) 


# assign genes
CLCA_variant_CES[, gene := fifelse(is.na(gene),ov$annot$symbol[match(variant_id, ov$variant_id)], gene)]
CLCA_variant_CES$Consequence = clca_maf$Consequence[match(CLCA_variant_CES$variant_id, clca_maf$variant_id)]
CLCA_variant_CES$Classification = clca_maf$Variant_Classification[match(CLCA_variant_CES$variant_id, clca_maf$variant_id)]
CLCA_variant_CES[  , region_type := region_vec[ match(variant_id, ov$variant_id) ]]

TCGA_variant_CES$Consequence = tcga_maf$Consequence[match(TCGA_variant_CES$variant_id,  tcga_maf$variant_id)]
TCGA_variant_CES$Classification = tcga_maf$Variant_Classification[match(TCGA_variant_CES$variant_id,  tcga_maf$variant_id)]

# Plot the CES for HCC-CLCA
plot_effects(
  effects   = CLCA_variant_CES,
  topn      = 10,
  group_by = "gene",  # group by gene
  label_individual_variants = T,
  prevalence_method = "both",
  #order_by_effect = F,
  #color_by  = "#4C72B0"  ,
  #color_by = "Consequence",
  color_by = "region_type",
  title = "HCC-CLCA: Top 10 Genes by CES (all variants)"
)

plot_effects(
  effects   = CLCA_variant_CES[variant_type == "aac", ],
  topn      = 10,
  group_by = "gene",  # group by gene
  label_individual_variants = T,
  prevalence_method = "both",
  #order_by_effect = F,
  color_by  = "#4C72B0",
  title = "HCC-CLCA: Top 10 Genes by CES (Coding)"
)

plot_effects(
  effects = TCGA_variant_CES,
  topn = 10,
  group_by = "gene",  # group by gene
  label_individual_variants = T,
  prevalence_method = "both",
  #order_by_effect = F,
  #color_by = "Consequence",
  title = "LIHC-TCGA: Top 10 Genes by CES (all variants)"
)



