library(dplyr)
library(purrr)
library(ggplot2)
library(rtracklayer)
library(GenomicFeatures)


# 1. Carregando GTF ---------------------------------------
gtf <- "./data/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz"

gtf_data <- import(gtf)

# 2. Extraindo métricas do GTF ----------------------------

n_genes <- length(unique(gtf_data$gene_id))

n_biotipos_gene <- gtf_data[, c("gene_id", "gene_biotype")] %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_biotype) %>%
  distinct() %>%
  group_by(gene_biotype) %>%
  count() %>%
  ungroup()

n_transcripts <- length(unique(gtf_data$transcript_id))

n_biotipos_transcripts <-
  gtf_data[, c("transcript_id", "transcript_biotype")] %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, transcript_biotype) %>%
  distinct() %>%
  group_by(transcript_biotype) %>%
  count() %>%
  ungroup()

rm(gtf_data)
gc()

# 3. Lendo dados com biotipos ----------------------------
dge_w_biotypes <-
  readr::read_csv("results/diff_exp/dge_w_biotype.csv") %>%
  dplyr::select(gene_id, gene_biotype) %>%
  distinct()
dte_w_biotypes <-
  readr::read_csv("results/diff_exp/dte_w_biotype.csv") %>%
  dplyr::select(transcript_id, transcript_biotype) %>%
  distinct()
dtu_w_biotypes <-
  readr::read_csv("results/diff_exp/dtu_w_biotype.csv") %>%
  dplyr::select(isoform_id, iso_biotype) %>%
  distinct()

# 2. Calculando os pvalores com phyper ----------------------------
# N -> Número de genes total = Universo
# M -> Número de genes de uma via
# n -> Número de genes da sua lista de interesse
# k -> Número de genes na interseção entre os genes da via e os genes da sua lista

calculate_phyper <-
  function(N_total,
           biotipos_df,
           tag_biotipos_df,
           biotipo_col,
           id_col,
           biotipo) {
    biotipo_col <- enquo(biotipo_col)
    id_col <- enquo(id_col)
    
    M <- biotipos_df %>%
      filter(!!biotipo_col == biotipo) %>%
      pull(n)
    
    n_lista <- tag_biotipos_df %>%
      pull(!!id_col) %>%
      unique() %>%
      length()
    
    k <- tag_biotipos_df %>%
      group_by(!!biotipo_col) %>%
      summarise(biotype_n = n()) %>%
      ungroup() %>%
      filter(!!biotipo_col == biotipo) %>%
      pull(biotype_n)
    
    
    data.frame(biotipo_name = biotipo,
               pvalor = 1 - phyper(k - 1, M, N_total - M, n_lista))
  }

dge_enriched <- dge_w_biotypes %>%
  purrr::map_dfr(
    .x = unique(dge_w_biotypes$gene_biotype),
    .f = ~ calculate_phyper(
      N_total = n_genes,
      biotipos_df = n_biotipos_gene,
      tag_biotipos_df = dge_w_biotypes,
      biotipo_col = gene_biotype,
      id_col = gene_id,
      biotipo = .x
    )
  )

dte_enriched <- dte_w_biotypes %>%
  purrr::map_dfr(
    .x = unique(dte_w_biotypes$transcript_biotype),
    .f = ~ calculate_phyper(
      N_total = n_transcripts,
      biotipos_df = n_biotipos_transcripts,
      tag_biotipos_df = dte_w_biotypes,
      biotipo_col = transcript_biotype,
      id_col = transcript_id,
      biotipo = .x
    )
  )

dtu_enriched <- dtu_w_biotypes %>%
  purrr::map_dfr(
    .x = unique(dtu_w_biotypes$iso_biotype),
    .f = ~ calculate_phyper(
      N_total = n_transcripts,
      biotipos_df = n_biotipos_transcripts %>%
        # Pq o resultado do isoform vem com um nome de coluna diferente
        dplyr::rename(iso_biotype = transcript_biotype),
      tag_biotipos_df = dtu_w_biotypes,
      biotipo_col = iso_biotype,
      id_col = isoform_id,
      biotipo = .x
    )
  )


all_enriched <- dge_enriched %>%
  mutate(group = "DGE") %>%
  bind_rows(dte_enriched %>% mutate(group = "DTE")) %>%
  bind_rows(dtu_enriched %>% mutate(group = "DTU")) %>%
  mutate(padjust = p.adjust(pvalor, method = "BH"))

save(dge_enriched, dte_enriched, dtu_enriched, all_enriched, file = "results/diff_exp/enriched_biotypes.rda")

all_enriched %>%
  mutate(log10p = -log10(padjust)) %>%
  ggplot(aes(x = log10p, y = biotipo_name, fill = group)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis_d() +
  labs(y = NULL, x = "-log10(padj)")
