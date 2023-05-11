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
  mutate(phenotype_reg = ifelse(phenotype == "MDD", 1, 0),
         run = rownames(.)) %>%
  dplyr::select(-c(group))

genes_interest <- diff_df %>%
  filter(type == "DGE") %>%
  dplyr::select(gene, hgnc_symbol) %>%
  distinct()

dge_tpms <-
  get_expr_from_gene_list(genes_interest$gene, genes_interest$hgnc_symbol)

full_data <- dge_tpms %>%
  left_join(metadata, by = "run")

selected_genes <- c("ATAT1", "DDX39B")

two_genes_df <- diff_df %>%
  filter(hgnc_symbol == selected_genes) %>%
  dplyr::select(gene, hgnc_symbol) %>%
  distinct()

selected_tpms <-
  get_expr_from_gene_list(two_genes_df$gene, selected_genes)

two_gene_data <- selected_tpms %>%
  left_join(metadata, by = "run")

if(!dir.exists("results/sym_reg/")) {
  dir.create("results/sym_reg/")
}

readr::write_tsv(two_gene_data, "results/sym_reg/selected_genes_for_reg.tsv")
readr::write_tsv(full_data, "results/sym_reg/genes_for_reg.tsv")
