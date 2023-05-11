library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)

# 1. Load GTF ---------------------------------------
gtf <- "./data/genome/Homo_sapiens.GRCh38.97.gtf.gz"

gtf_data <- import(gtf)

# 2. Import DGE and DTE results and get feature biotypes -------------
load("./results/diff_exp/diff_df.rda")

dge_genes <- diff_df %>%
  filter(type == "DGE")

# 2.1 Get DGE biotypes -----------------------------------------------
dge_w_biotype <- gtf_data[, c("gene_id", "gene_biotype")] %>%
  as.data.frame() %>%
  filter(gene_id %in% dge_genes$gene) %>%
  dplyr::select(gene_id, gene_biotype) %>%
  right_join(dge_genes, by = c("gene_id" = "gene")) %>%
  distinct() %>% 
  dplyr::select(gene_id, gene_biotype, group)

readr::write_csv(dge_w_biotype, "results/diff_exp/dge_w_biotype.csv")

# 2.2. Get DTE biotypes ----------------------------------------------
load("./results/diff_exp/diff_tx_corrected.rda")

dte_genes <- df_res_padj_tx %>%
  dplyr::select(txID, transcript, group) %>%
  filter(transcript < 0.05)

dte_w_biotype <-
  gtf_data[, c("transcript_id", "transcript_biotype")] %>%
  as.data.frame() %>%
  filter(transcript_id %in% dte_genes$txID) %>%
  dplyr::select(transcript_id, transcript_biotype) %>%
  right_join(dte_genes, by = c("transcript_id" = "txID")) %>%
  distinct() %>% 
  dplyr::select(transcript_id, transcript_biotype, group)

readr::write_csv(dte_w_biotype, "results/diff_exp/dte_w_biotype.csv")
 
# 2.3. Get DTU biotypes ---------------------------------------------

dtu_w_biotype <- readr::read_csv("results/ISA/dtu_w_biotype.csv")

# Plot  ---------------------------------------
plot_biotype_bar <- function(data, id_col, n_col, color) {
  
  id_col <- enquo(id_col)
  n_col <- enquo(n_col)
  
  data %>% 
    ggplot(aes(x = reorder(!!id_col, dplyr::desc(!!n_col)), y = !!n_col)) +
    geom_col(fill = color) + 
    scale_y_continuous(labels = scales::percent_format(scale = 1), name = "", limits = c(0, 100)) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(-1, 0, -1, 0))
}

dge_plot <- dge_w_biotype %>% 
  group_by(gene_biotype) %>% 
  summarise(biotype_n = n() / length(unique(dge_w_biotype$gene_id)) * 100) %>% 
  ungroup() %>% 
  mutate(type = "DGE") %>% 
  dplyr::rename(biotype = gene_biotype)

dte_plot <- dte_w_biotype %>%
  group_by(transcript_biotype) %>%
  summarise(biotype_n = n() / length(unique(dte_w_biotype$transcript_id))* 100) %>%
  ungroup() %>%
  mutate(type = "DTE") %>% 
  dplyr::rename(biotype = transcript_biotype)

dtu_plot <- dtu_w_biotype %>%
  group_by(iso_biotype) %>%
  summarise(biotype_n = n() / length(unique(dtu_w_biotype$isoform_id))* 100) %>%
  ungroup() %>%
  mutate(type = "DTU") %>%
  dplyr::rename(biotype = iso_biotype)

df_plot <- Reduce(bind_rows, list(dge_plot, dte_plot, dtu_plot))
df_plot$biotype <- gsub("_", " ", df_plot$biotype)

# Plot feature biotypes

color_scale <- c("DGE" = "#0ac80aff", "DTE" = "#4f4affff", "DTU" = "#ff822fff")

ggplot(df_plot, aes(x = reorder(biotype, dplyr::desc(biotype_n)), y = biotype_n, fill = type)) +
  geom_col(show.legend = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  facet_wrap(nrow = 3, ncol = 1, facets = vars(type), scales = "free_y", strip.position = "right") +
  scale_fill_manual(values = color_scale) + 
  labs(x = "Feature biotypes", y = "% of feature biotype by the total features") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) -> biotype_plot


# Save
ggsave(biotype_plot, file = "results/plots_paper/biotype_plot.pdf", width = 7, height = 4)

