library(dplyr)

load("results/diff_exp/diff_df.rda")

load("results/txi/txi_gene.rda")

load("results/important_variables/ann.rda")

genes_interest <- diff_df %>%
  filter(type == "DGE") %>%
  select(gene, hgnc_symbol) %>%
  distinct()

tpms <-
  txi$abundance[rownames(txi$abundance) %in% genes_interest$gene, ] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  setNames(., c("run", genes_interest$hgnc_symbol))

full_data <- tpms %>%
  left_join(ann, by = "run") %>%
  mutate(# ??
    ph = ph[, 1],
    rin = rin[, 1]) %>%
  mutate(phenotype_reg = ifelse(phenotype == "MDD", 1, 0),) %>%
  select(-c(group))

readr::write_tsv(full_data, "results/sym_reg/genes_for_reg.tsv")
