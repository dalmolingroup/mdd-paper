library(edgeR)
library(ppcseq)
library(tidyverse)
library(magrittr)

# Differentially expressed genes from edgeR -------------------------------

load("results/diff_exp/edger_gene_rin_ph_diff.rda")

# Sample annotation  ------------------------------------------------------

load("results/important_variables/ann.rda")

# Counts for each gene normalized in TMM ----------------------------------

load("results/txi/txi_gene.rda")

# PPCSEQ has methods to normalize counts, use raw estimates
counts <- txi$counts

##########################################################################

# Identify outliers genes in each comparison group

ann$run <- rownames(ann)

comp <- paste(rep(unique(ann$region), each = 2), unique(ann$gender), sep = "_")

map(comp, function(c) {
  
  # Diff genes for this comparison
  diff_genes <- df_edger_ph_rin_group_gene %>% 
    filter(group == c) %>% 
    pull(gene)
  
  # Select information for each region/sex
  ann2 <- ann %>% 
    mutate(group = paste(region, gender, sep = "_")) %>% 
    filter(group == c) %>% 
    dplyr::select(run, group, ph, rin, phenotype) %>% 
    mutate(phenotype = as.factor(phenotype))
  
  # Select counts of diff expressed genes from sex/region comparison
  counts2 <- counts %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% diff_genes)
  
  # Transform count data
  counts2 <- counts2 %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene") %>% 
    gather(key = "run", value = "count", -gene)
  
  # Join data and annotation info
  counts2 <- inner_join(counts2, ann2, by = "run")
  
  # Results from edgeR analysis
  df_res <- df_edger_ph_rin_group_gene %>% 
    filter(group == c) %>% 
    dplyr::select(gene, pval = PValue)
  
  ppc_table <- inner_join(counts2, df_res, by = "gene") %>%
    dplyr::select(run, gene, count, phenotype, ph, rin, pval) %>% 
    filter(!is.na(pval))
  
  counts.ppc <- ppc_table %>%
    mutate(is_significant = pval < 0.05,
           count = as.integer(count+1)) %>% # Adding 1 to prevent division by 0
    identify_outliers(
      formula = ~ phenotype,
      .sample = run, 
      .transcript = gene,
      .abundance = count,
      .significance = pval,
      .do_check = is_significant,
      percent_false_positive_genes = 5, 
      approximate_posterior_inference = FALSE,
      cores = 8
   )
  
  return(counts.ppc)
  
}) -> ppc_outliers

names(ppc_outliers) <- comp

# Select only the genes with outliers for further removal
imap_dfr(ppc_outliers, function(x, y) {

  idx <- which(x[["tot_deleterious_outliers"]] != 0)

  if(length(idx) > 0) {

    res <- x %>%
      slice(idx) %>%
      select(gene, sample_wise_data) %>%
      unnest(cols = sample_wise_data) %>%
      select(gene, run, deleterious_outliers) %>%
      filter(deleterious_outliers) %>%
      mutate(group = y) %>%
      select(gene, run, group)

  }

}) -> outliers_samples_dge

save(outliers_samples, file = "results/diff_exp/outliers_samples_dge.rda")
