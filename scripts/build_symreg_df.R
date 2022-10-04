library(dplyr)

load("results/diff_exp/diff_df.rda")

load("results/txi/txi_gene.rda")

load("results/important_variables/ann.rda")

selected_genes <- c("ATAT1", "DDX39B", "ABR")

genes_interest <- diff_df %>%
  filter(hgnc_symbol == selected_genes) %>%
  select(gene, hgnc_symbol) %>%
  distinct()

tpms <-
  txi$abundance[rownames(txi$abundance) %in% genes_interest$gene,] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  setNames(., c("run", selected_genes))

full_data <- tpms %>%
  left_join(ann, by = "run") %>%
  mutate(# ??
    ph = ph[, 1],
    rin = rin[, 1]) %>%
  mutate(
    phenotype_reg = ifelse(phenotype == "MDD", 1, -1),
    gender_num = ifelse(gender == "female", 1, 0),
    # value = 1
  ) %>%
  select(-c(group))

  # tidyr::spread(region, value, fill = 0) %>%
  # select(-c(group, `<NA>`))

readr::write_tsv(full_data, "results/sym_reg/selected_genes_for_reg.tsv")


