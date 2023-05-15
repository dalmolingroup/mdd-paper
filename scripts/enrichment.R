# ENRICHMENT ANALYSIS -----------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(TreeAndLeaf)
library(RedeR)
library(igraph)
library(GOSemSim)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
library(ggpubr)

# CREATE FUNCTIONS

# Convert IDS
convert_id <- function(data, col) {
  data %>% 
    pull({{col}}) %>% 
    unique() %>% 
    bitr(fromType = "ENSEMBL",
         toType = "ENTREZID",
         OrgDb = "org.Hs.eg.db") %>% 
    mutate(group = data$group[1]) %>%
    mutate(sex = data$sex[1])
}

# Enrichment with enrichGO function
enrich <- function(entrezgenes, col, id) {
  
  res <- enrichGO(gene          = entrezgenes[[col]],
                  OrgDb         = org.Hs.eg.db,
                  keyType       = id,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE) # Uncomment to change ENSEMBL id to SYMBOL
}

# Save plots
create_plots <- function(enrichment, annotation, directory){
    names(enrichment) <- lapply(annotation, function(x){
     unique(as.vector(x$group))
    }) %>% unlist() %>% unname()
    
    if(!dir.exists(directory)) {
      dir.create(directory, recursive = T)
    }
    
    for(i in 1:12){
      if(is.null(enrichment[i][[1]])){
        message(paste(names(enrichment)[i], "is NULL",  sep= ' '))
        next()
      }
      else if(sum(enrichment[[i]]@result$p.adjust < 0.05) == 0){
        message(paste(names(enrichment)[i], "doesn't have enriched terms",  sep= ' '))
        next()
      } else{
        suppressMessages(dotplot(enrichment[[i]], showCategory=20) + scale_color_viridis_c())
        ggsave(filename = paste0(directory, "enrichment_", names(enrichment)[i], ".pdf"),
               width=12, height=8, plot = last_plot())
      }
    }

    return(enrichment)
  }

# Transform results to Dataframe
transform_to_df <- function(enrichment, padj_cutoff = 0.05){
  list_of_df <- list()
  for(i in 1:12){
    if(is.null(enrichment[i][[1]])){
      print(NULL)
    } else{
      enrichment_df <- as_tibble(enrichment[[i]]@result) %>% filter(p.adjust < padj_cutoff) %>% 
      mutate(group = names(enrichment)[i]) %>% list()
      list_of_df <- append(list_of_df, enrichment_df)
    }
  }
  return(list_of_df)
}

# Save dataframes to csv files
save_df <- function(enrichment, file){
  enrichment %>%
    bind_rows() %>%
    mutate(qvalue = scales::scientific(qvalue)) %>% 
    write_csv(file=file)
}


# ENRICHMENT OF DIFFERENTIAL GENES ---------------------------------------
load("results/diff_exp/diff_df.rda")
# 
# split_genes <- diff_df %>%group_split(group)

split_genes <- split(diff_df, diff_df$group)

converted <- lapply(split_genes, convert_id, col = "gene")

enriched <- lapply(converted, enrich, col = "ENTREZID", id = "ENTREZID")

enriched_ann <- create_plots(enrichment = enriched, directory = 'results/tx_enrich/', annotation = converted)

enriched_df_diff <- transform_to_df(enriched_ann)

names(enriched_df_diff) <- names(split_genes)

save_df(enriched_df_diff, file = 'results/tx_enrich/go_terms.csv')

save(enriched_df_diff, file = "results/tx_enrich/go_terms_tx_by_group.rda")

# ## INTERSECTION ---------------------------------------------------------------
# load("results/intersects/intersects.rda")
# 
# split_genes <- genes_by_sex %>% 
#   group_split(intersect)
# 
# converted <- lapply(split_genes, convert_id, col = "genes")
# 
# enriched <- lapply(converted, enrich, col = "ENTREZID", id = "ENTREZID")
# 
# enriched_ann <- create_plots(enrichment = enriched, directory = 'results/tx_enrich/intersection_by_sex/', annotation = converted)
# 
# enriched_df_by_sex <- transform_to_df(enriched_ann)
# 
# names(enriched_df_by_sex) <- c('female', 'female:male', 'male')
# 
# save_df(enriched_df_by_sex, file = 'results/tx_enrich/intersection_by_sex/go_terms_by_sex.csv')
# 
# save(enriched_df_by_sex, file = "results/tx_enrich/intersection_by_sex/go_terms_by_sex.rda")


# PLOT TREE AND LEAF ------------------------------------------------------

# Organize data

cutoff <- 0
sizeIntervals <- 5 # usar 3 para o idx 8 e 5 para o 9

organize_data <- function(enriched_df_diff, contrast) {
  
  data <- enriched_df_diff[[contrast]]
  terms <- data$ID
  
  semData <- godata(ont = "BP")
  data$path_lenght <- as.integer(sapply(strsplit(data$GeneRatio, "/"), "[", 2))
  data$ratio <- data$Count/data$path_lenght
  mSim <- mgoSim(terms, terms, semData = semData, measure = "Wang", combine = NULL)
  
  if(cutoff > 0){
    RS <- rowSums(mSim) - 1
    GOs <- names(sort(RS)[-1:-cutoff])
    mSim <- mSim[rownames(mSim) %in% GOs, colnames(mSim) %in% GOs]
    terms <- terms[terms %in% GOs]
    data <- data[data$ID %in% GOs,]
  }
  
  return(list(msim = mSim, data = data))
  
}

cols <- brewer.pal(9, "OrRd")

# Nac male ----------------------------------------------------------------

nac_male <- organize_data(enriched_df_diff, "Nac_male")

# Create tree and leaf object
hc <- hclust(dist(nac_male$msim), "average")
tal_nac <- treeAndLeaf(hc)
class(tal_nac) <- "igraph"

# Modify graph
gg <- as_tbl_graph(tal_nac, directed = F) %>% 
  activate(nodes) %>% 
  left_join(nac_male$data[, c("ID", "Count", "ratio", "Description")], by = c("name" = "ID")) %>% 
  mutate(Count = ifelse(is.na(Count), 0, Count),
         ratio = ifelse(is.na(ratio), 0, ratio))

# plot
g_nac_male <- ggraph(gg, "unrooted", length = 1) + 
  geom_edge_link(edge_color = "gray", width = 0.7, n = 10) + 
  geom_node_point(aes(size = Count), stroke = 1) + 
  geom_node_point(aes(size = Count, col = ratio)) + 
  geom_node_text(aes(label = Description),nudge_y = 0.01) +
  scale_color_gradient(low = cols[1], high = cols[length(cols)]) + 
  scale_size_continuous(breaks = seq(2, 18, 6), limits = c(2, 18)) + 
  labs(title = "Nac: Male") + 
  theme(element_text(family = "Arial")) +
  theme_graph()


# OFC female --------------------------------------------------------------

ofc_female <- organize_data(enriched_df_diff, "OFC_female")

# Create tree and leaf object
hc <- hclust(dist(ofc_female$msim), "average")
tal_ofc <- treeAndLeaf(hc)
class(tal_ofc) <- "igraph"

# Modify graph
gg <- as_tbl_graph(tal_ofc, directed = F) %>% 
  activate(nodes) %>% 
  left_join(ofc_female$data[, c("ID", "Count", "ratio", "Description")], by = c("name" = "ID")) %>% 
  mutate(Count = ifelse(is.na(Count), 0, Count),
         ratio = ifelse(is.na(ratio), 0, ratio))

# plot
g_ofc_female <- ggraph(gg, "unrooted", length = 0.1) + 
  geom_edge_link(edge_color = "gray",  width = 0.7, n = 10) + 
  geom_node_point(aes(size = Count), stroke = 1) + 
  geom_node_point(aes(size = Count, col = ratio)) + 
  geom_node_text(aes(label = Description), nudge_y = 0.01) +
  scale_color_gradient(low = cols[1], high = cols[length(cols)]) + 
  scale_size_continuous(breaks = seq(2, 18, 6), limits = c(2, 18)) + 
  theme(element_text(family = "Arial")) +
  labs(title = "OFC: Female") + 
  theme_graph()


# aINS Female -------------------------------------------------------------

ains_female <- organize_data(enriched_df_diff, "aINS_female")

# Create tree and leaf object
hc <- hclust(dist(ains_female$msim), "average")
tal_ains <- treeAndLeaf(hc)
class(tal_ains) <- "igraph"

# Modify graph
gg <- as_tbl_graph(tal_ains, directed = F) %>% 
  activate(nodes) %>% 
  left_join(ains_female$data[, c("ID", "Count", "ratio", "Description")], by = c("name" = "ID")) %>% 
  mutate(Count = ifelse(is.na(Count), 0, Count),
         ratio = ifelse(is.na(ratio), 0, ratio))

# plot
g_ains_female <- ggraph(gg, "unrooted", length = 0.1) + 
  geom_edge_link(edge_color = "gray",  width = 0.7, n = 10) + 
  geom_node_point(aes(size = Count), stroke = 1) + 
  geom_node_point(aes(size = Count, col = ratio)) + 
  geom_node_text(aes(label = Description), nudge_y = 0.01) +
  scale_color_gradient(low = cols[1], high = cols[length(cols)]) + 
  scale_size_continuous(breaks = seq(2, 18, 6), limits = c(2, 18)) + 
  theme(element_text(family = "Arial")) +
  labs(title = "aINS: Female") + 
  theme_graph()


# Arrange graphs
ggarrange(g_nac_male, g_ofc_female, g_ains_female, common.legend = T)

ggsave("results/plots_paper/tal_enrichment.svg", width = 10, height = 7, dpi = 300)
