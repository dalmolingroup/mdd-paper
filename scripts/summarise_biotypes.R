library(dplyr)
library(purrr)
library(tidyr)
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
  dplyr::select(gene_id, gene_biotype, group) %>% 
  dplyr::rename(biotype = gene_biotype) %>% 
  mutate(type = "DGE")

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
  dplyr::select(transcript_id, transcript_biotype, group) %>% 
  dplyr::rename(biotype = transcript_biotype) %>% 
  mutate(type = "DTE")

readr::write_csv(dte_w_biotype, "results/diff_exp/dte_w_biotype.csv")
 
# 2.3. Get DTU biotypes ---------------------------------------------

dtu_w_biotype <- readr::read_csv("results/ISA/dtu_w_biotype.csv")
dtu_w_biotype <- dtu_w_biotype %>% 
  dplyr::rename(biotype = iso_biotype) %>% 
  mutate(type = "DTU")


# Plot for both sexes -----------------------------------------------------

dge_plot <- dge_w_biotype %>% 
  group_by(biotype) %>% 
  dplyr::summarise(biotype_n = dplyr::n() / length(unique(dge_w_biotype$gene_id)) * 100) %>% 
  ungroup() %>% 
  dplyr::mutate(type = "DGE")

dte_plot <- dte_w_biotype %>%
  group_by(biotype) %>%
  summarise(biotype_n = dplyr::n() / length(unique(dte_w_biotype$transcript_id))* 100) %>%
  ungroup() %>% 
  mutate(type = "DTE")

dtu_plot <- dtu_w_biotype %>%
  group_by(biotype) %>%
  summarise(biotype_n = dplyr::n() / length(unique(dtu_w_biotype$isoform_id))* 100) %>%
  ungroup() %>% 
  mutate(type = "DTU")

df_plot <- Reduce(bind_rows, list(dge_plot, dte_plot, dtu_plot))
df_plot$biotype <- gsub("_", " ", df_plot$biotype)

# Plot feature biotypes
color_scale <- c("DGE" = "#0ac80aff", "DTE" = "#4f4affff", "DTU" = "#ff822fff")

# Female and male together
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


# Plot for feature types for each sex -------------------------------------

dge_plot <- dge_w_biotype %>% 
  separate(group, into = c("region", "sex")) %>% 
  group_by(biotype, sex) %>% 
  summarise(biotype_n = dplyr::n()) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  mutate(prop = biotype_n / sum(biotype_n) * 100,
         type = "DGE")

dte_plot <- dte_w_biotype %>% 
  separate(group, into = c("region", "sex")) %>% 
  group_by(biotype, sex) %>% 
  summarise(biotype_n = dplyr::n()) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  mutate(prop = biotype_n / sum(biotype_n) * 100,
         type = "DTE")

dtu_plot <- dtu_w_biotype %>% 
  separate(group, into = c("region", "sex")) %>% 
  group_by(biotype, sex) %>% 
  summarise(biotype_n = dplyr::n()) %>% 
  ungroup() %>% 
  group_by(sex) %>% 
  mutate(prop = biotype_n / sum(biotype_n) * 100,
         type = "DTU")

df_plot <- Reduce(bind_rows, list(dge_plot, dte_plot, dtu_plot))
df_plot$biotype <- gsub("_", " ", df_plot$biotype)

# Plot feature biotypes
color_scale <- c("DGE" = "#0ac80aff", "DTE" = "#4f4affff", "DTU" = "#ff822fff")

# Female and male together
ggplot(df_plot, aes(x = reorder(biotype, dplyr::desc(prop)), y = prop, fill = type)) +
  geom_col(show.legend = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  facet_grid(rows = vars(type), cols = vars(sex), scales = "free_y",
             labeller = labeller(sex = as_labeller(c("female" = "Female", "male" = "Male")))) +
  scale_fill_manual(values = color_scale) + 
  labs(x = "Feature biotypes", y = "% of feature biotype by the total features") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white")) -> biotype_plot_by_sex

# Save
ggsave(biotype_plot_by_sex, filename = "results/plots_paper/fig3.pdf", width = 7, height = 4)

# Test feature prevalence differences between female and male -------------

biotypes_by_sex <- Reduce(bind_rows, list(
  dge_w_biotype %>% dplyr::select(-gene_id), 
  dte_w_biotype %>% dplyr::select(-transcript_id), 
  dtu_w_biotype %>% dplyr::select(-isoform_id)))

biotypes_by_sex %>% 
  separate(group, into = c("region", "sex"), sep = "_") %>% 
  arrange(type, biotype) %>% 
  group_by(type) %>% 
  group_map(~ {
    cat(.y$type, sep = "\n")
    cont_table <- table(.x$biotype, .x$sex)
    return(list(fisher = fisher.test(cont_table, simulate.p.value = T), count_table = cont_table))
  }) -> biot_tests_fisher

biotypes_by_sex %>% 
  separate(group, into = c("region", "sex"), sep = "_") %>% 
  arrange(type, biotype) %>% 
  group_by(type) %>% 
  group_map(~ {
    cat(.y$type, sep = "\n")
    cont_table <- table(.x$biotype, .x$sex)
    return(chisq.test(cont_table))
  }) -> biot_tests_chisq

names(biot_tests_fisher) <- c("DGE", "DTE", "DTU")
names(biot_tests_chisq) <- c("DGE", "DTE", "DTU")

# Divide plot by protein coding and non-coding ----------------------------
df_plot %>% 
  mutate(
    coding = case_when(
      biotype == "protein coding" ~ "coding",
      .default = "non-coding"
  )) %>% 
  group_by(coding, sex, type) %>% 
  mutate(prop_coding = sum(prop))  %>%
  dplyr::select(-biotype_n, -biotype, -prop) %>% 
  unique() -> df_plot_coding

ggplot(df_plot_coding, aes(x = reorder(coding, dplyr::desc(prop_coding)), y = prop_coding, fill = type)) +
  geom_col(show.legend = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  facet_grid(rows = vars(type), cols = vars(sex), scales = "free_y",
             labeller = labeller(sex = as_labeller(c("female" = "Female", "male" = "Male")))) +
  scale_fill_manual(values = color_scale) + 
  labs(x = "Feature biotypes", y = "% of feature biotype by the total features") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white")) -> biotype_plot_coding

ggsave(biotype_plot_coding, filename = "results/plots_paper/biotype_sex_coding.pdf", width = 7, height = 4)

# Check biotypes by region in females -------------------------------------

biotypes_by_sex %>% 
  separate(group, into = c("region", "sex"), sep = "_") %>% 
  arrange(type, biotype) %>% 
  filter(sex == "female") %>%  
  group_by(region, type) %>% 
  mutate(n1 = dplyr::n()) %>% 
  ungroup() %>% 
  group_by(biotype, type, region) %>% 
  mutate(n2 = dplyr::n(),
         prop_by_region = (n2 / n1) * 100) %>% 
  arrange(desc(type), desc(region)) %>% 
  ungroup() %>% 
  dplyr::select(biotype, region, type, prop_by_region) %>% 
  distinct() -> biotypes_female

ggplot(biotypes_female, aes(x = reorder(biotype, dplyr::desc(prop_by_region)), y = prop_by_region, fill = type)) +
  geom_col(show.legend = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  facet_grid(rows = vars(region), cols = vars(type), scales = "free_y") +
  scale_fill_manual(values = color_scale) + 
  labs(x = "Feature biotypes", y = "% of feature biotype by the total features") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white")) -> plot_biotypes_female

ggsave(plot_biotypes_female, filename = "results/plots_paper/biotype_female.pdf",  width = 10, height = 7)


# Check biotypes by region male -------------------------------------------

biotypes_by_sex %>% 
  separate(group, into = c("region", "sex"), sep = "_") %>% 
  arrange(type, biotype) %>% 
  filter(sex == "male") %>%  
  group_by(region, type) %>% 
  mutate(n1 = dplyr::n()) %>% 
  ungroup() %>% 
  group_by(biotype, type,region) %>% 
  mutate(n2 = dplyr::n(),
         prop_by_region = (n2 / n1) * 100) %>% 
  arrange(desc(type), desc(region)) %>% 
  ungroup() %>% 
  dplyr::select(biotype, region, type, prop_by_region) %>% 
  distinct() -> biotypes_male

ggplot(biotypes_male, aes(x = reorder(biotype, dplyr::desc(prop_by_region)), y = prop_by_region, fill = type)) +
  geom_col(show.legend = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  facet_grid(rows = vars(region), cols = vars(type), scales = "free_y") +
  scale_fill_manual(values = color_scale) + 
  labs(x = "Feature biotypes", y = "% of feature biotype by the total features") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white")) -> plot_biotypes_male

ggsave(plot_biotypes_male, filename = "results/plots_paper/biotype_male.pdf",  width = 10, height = 7)


# -------------------------------------------------------------------------
