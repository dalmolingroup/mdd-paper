
# GWAS risk variants ------------------------------------------------------

library(dplyr)
library(purrr)

# Get studies related to 'major depressive disorder' and 'unipolar depression'.
efo_id <- c(
  unipolar_depression        = "EFO_0003761",
  major_depressive_disorder  = "MONDO_0002009"
)
mdd_results <- get_studies(efo_id = efo_id, set_operation = "intersection") 

mdd_filtered <- mdd_results@publications %>% 
  filter(publication_date > "2018-01-01")

# From the results above, we manually removed studies that involved different themes,
# such as other psychiatric conditions. The results are shown in 'mdd_filtered.csv'
mdd_filtered <- read.csv("results/tables/mdd_filtered.csv")
studies <- unique(mdd_filtered$study_id)

# Gather results from association and variation annotation from the GWAS Catalog.
map_dfr(studies, function(study_id) {
  
  variants <- get_variants(study_id = study_id, set_operation = "intersection")

  variant_gc_data <- variants@genomic_contexts %>%
    dplyr::select(variant_id, chromosome_name, chromosome_position)

  variant_var_data <- variants@variants %>%
    dplyr::select(variant_id, functional_class)

  variant_data <- inner_join(variant_gc_data, variant_var_data, by = c("variant_id"))
  
  associations <- get_associations(study_id = study_id, set_operation = "intersection")
  
  assoc_1 <- associations@risk_alleles %>% 
    dplyr::select(association_id, variant_id, risk_allele, risk_frequency, genome_wide)
  assoc_2 <- associations@associations %>% 
    dplyr::select(association_id, pvalue, pvalue_description, range, beta_number, beta_direction)

  associations_data <- reduce(list(assoc_1, assoc_2), inner_join, by = "association_id")
  
  associations_data

  inner_join(associations_data, variant_data, by = "variant_id") %>% 
    distinct()
  
}) -> res_gwas

# Filter alleles by their level of significance < 1e-6 
chrom <- c(1:22, "X", "Y", "MT")
res_gwas %>% 
  filter(!is.na(risk_allele), pvalue < 1e-6, chromosome_name %in% chrom) -> risk_alleles

# Get ensembl variation information
ensembl <- useEnsembl("snps", "hsapiens_snp")

# Get the correct gene names associated to each variant.
# This is important because the GWAS Catalog API doesn't provide consistent
# information about the gene-variant annotation. 
alleles <- unique(risk_alleles$variant_id)

dict <- getBM(attributes = c("ensembl_gene_name", "refsnp_id"),
              filters = "snp_filter",
              values = alleles,
              mart = ensembl)

# Join information from Ensembl variation
risk_alleles <- risk_alleles %>% 
  left_join(dict, by = c("variant_id" = "refsnp_id"))

# Information about the risk alleles
n_distinct(risk_alleles$variant_id)

# Gather information of TAGs and genes annotated for risk alleles
load("results/diff_exp/diff_df.rda")
intersection <- inner_join(risk_alleles, diff_df, by = c("ensembl_gene_name" = "gene"))

# Save 
save(intersection, file = "results/diff_exp/gwas_intersections.rda")
write.csv(intersection, file = "results/tables/gwas_intersection.csv", row.names = F, quote = F)
