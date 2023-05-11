library(tidyverse)
library(edgeR)

load("results/diff_exp/diff_tx_corrected.rda")
load("results/diff_exp/outliers_samples_dge.rda")
load("results/diff_exp/outliers_samples_dte.rda")

# Organize DTE data -------------------------------------------------------

# Get logFC data
# This information is not stored because the stageR step requires only the p-values
load("results/diff_exp/edger_tx_rin_ph_diff.rda")
imap_dfr(lrt_comp, function(x, y) {
    df <- topTags(x, n = Inf)$table
    df <- df %>%
        tibble::rownames_to_column("tx") %>% 
        mutate(tx = gsub("\\.+\\d+", "", tx),
               group = y) %>% 
        dplyr::select(tx, logFC, group)
    return(df)
}) -> logFC_tx

df_res_padj_tx <- inner_join(df_res_padj_tx, logFC_tx, by = c("txID" = "tx", "group"))

# Filter out the outliers identified by ppcseq
df_res_padj_tx_out_filtered <- anti_join(df_res_padj_tx, outliers_samples_dte, by = c("txID" = "tx", "group"))

# Organize DGE data -------------------------------------------------------

load("results/diff_exp/edger_gene_rin_ph_diff.rda")
df_res_padj_gene_out_filtered <- anti_join(df_edger_ph_rin_group_gene, outliers_samples_dge, by = c("gene", "group"))

save(df_res_padj_gene_out_filtered, df_res_padj_tx_out_filtered, file = "results/diff_exp/df_res_dge_dte.rda")





