# Load the ces and reference genome package
library(cancereffectsizeR)
library(ces.refset.hg38)
library(data.table)

# Load the data
# Get file for liftOver
chain_file = 'hg19ToHg38.over.chain'
if(! file.exists(chain_file)) {
  download.file(url = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
                destfile = 'hg19ToHg38.over.chain.gz' )
  writeLines(readLines('hg19ToHg38.over.chain.gz'), chain_file)
}

#HCC-CLCA
# somatic variant MAF data loading (WGS)
clca_maf_file <- fread("./hcc_clca_2024/data_mutations.txt", sep = "\t", quote = "")
clca_maf <- preload_maf(maf = clca_maf_file, chain_file = chain_file, 
                        refset = "ces.refset.hg38", keep_extra_columns = "Variant_Classification")
# leave only somatic mutations
# and remove all repetitive regions calls but those at COSMIC-annotated sites 
# (where variants are assigned into tiers 1-3).
clca_maf = clca_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]
clca_maf = clca_maf[is.na(problem)]


#LIHC-TCGA
# somatic variant MAF data loading (WXS)
tcga_maf_file <- fread("./lihc_tcga/data_mutations.txt", sep = "\t", quote = "")
tcga_maf = preload_maf(maf = tcga_maf_file, chain_file = chain_file, 
                       coverage_intervals_to_check = ces.refset.hg38$default_exome,
                       refset = "ces.refset.hg38", keep_extra_columns = "Variant_Classification")
#remove problematic
tcga_maf = tcga_maf[is.na(problem)]
#we'll remove records >100 bp out-of-coverage.
tcga_maf = tcga_maf[dist_to_coverage_intervals <= 100]
# Again, leave only somatic mutations
# and remove all repetitive regions calls but those at COSMIC-annotated sites. 
tcga_maf = tcga_maf[germline_variant_site == F][repetitive_region == F | cosmic_site_tier %in% 1:3]

# clean maf output
fwrite(clca_maf, 'hcc_clca_hg38.maf', sep = "\t")
fwrite(tcga_maf, 'lihc_tcga_hg38.maf', sep = "\t")
