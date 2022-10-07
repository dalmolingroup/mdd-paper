library(dplyr)

load("results/diff_exp/diff_df.rda")

load("results/txi/txi_gene.rda")

load("results/important_variables/ann.rda")

get_expr_from_gene_list <- function(genelist, symbol_list) {
  txi$abundance[rownames(txi$abundance) %in% genelist,] %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    setNames(., c("run", symbol_list))
}

metadata <- ann %>%
  mutate(# ??
    ph = ph[, 1],
    rin = rin[, 1]) %>%
  mutate(phenotype_reg = ifelse(phenotype == "MDD", 1, 0),) %>%
  select(-c(group))

genes_interest <- diff_df %>%
  filter(type == "DGE") %>%
  select(gene, hgnc_symbol) %>%
  distinct()

dge_tpms <-
  get_expr_from_gene_list(genes_interest$gene, genes_interest$hgnc_symbol)

full_data <- dge_tpms %>%
  left_join(metadata, by = "run")

selected_genes <- c("ATAT1", "DDX39B", "ABR")

three_genes_df <- diff_df %>%
  filter(hgnc_symbol == selected_genes) %>%
  select(gene, hgnc_symbol) %>%
  distinct()

selected_tpms <-
  get_expr_from_gene_list(three_genes_df$gene, selected_genes)

three_gene_data <- selected_tpms %>%
  left_join(metadata, by = "run")

readr::write_tsv(three_gene_data, "results/sym_reg/selected_genes_for_reg.tsv")
readr::write_tsv(full_data, "results/sym_reg/genes_for_reg.tsv")
