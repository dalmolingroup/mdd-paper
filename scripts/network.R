
# PPI networks for each group ---------------------------------------------

library(biomaRt)
library(tidyverse)
library(igraph)
library(tidygraph)
library(scatterpie)
library(ggraph)
library(magrittr)
library(RedeR)

# ------------------ GET GENES AND INTERACTIONS ---------------------------
# Load diff genes table
load("results/diff_exp/diff_df.rda")
gwas_intersections <- read_csv("results/tables/gwas_intersection.csv")
gwas_genes <- unique(gwas_intersections$hgnc_symbol)

# FUNCTIONS ---------------------------------------------------------------

# Get interaction network from string database
get_map  <- function(ids) {
  s_map <- read.table(
    text = RCurl::postForm(
      "https://string-db.org/api/tsv/network"
      ,identifiers = ids
      ,echo_query  = "1"
      ,required_score = "0"
      ,show_query_node_labels = "1"
      ,species     = "9606"
    )[[1]], sep = "\t", header = T)
}

# Function to combine scores from different channels
combinescores <- function(dat, evidences = "all", confLevel = 0.4){
  
  if(evidences[1] == "all"){
    edat <- dat[, -c(1, 2, ncol(dat))]
  } else {
    if(!all(evidences %in% colnames(dat))) {
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat <- dat[, evidences]
  }
  edat <- edat/1000
  
  edat <- 1 - edat
  sc <- apply(X = edat, MARGIN = 1, FUN = function(x) 1 - prod(x))
  dat <- cbind(dat[, c(1, 2)], combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <- dat[idx,]
  return(dat)
  
}

# Get interaction network
identifiers <- paste0(unique(diff_df$hgnc_symbol), collapse = "%0d")
s_map <- get_map(identifiers)

# Select channels (from experimental data, coexpression, and database) and combine their scores
s_map %>% 
  mutate(across(contains("score"), function(i) round(i * 1000, digits = 0))) %>% 
  dplyr::select(contains("preferredName"), contains("score")) %>% 
  combinescores(evidences = c("escore", "ascore", "dscore"), confLevel = 0.4) %>% 
  unique() -> int

save(int, file = "results/networks/int.rda")

# RedeR

# Here we select proper coordinates to our network
g <- graph_from_edgelist(as.matrix(int[,1:2]), directed = F)
rdp <- RedeR::RedPort()
calld(rdp)
addGraph(rdp, g)

# Using REDER, we selected the best visualization for our network.
# Nodes and edges coordinates are saved as nodes2 and edges2 in 'results/networks' directory.

# Plot network layout ----------------------------------------------------

nodes <- read_tsv("results/networks/model_nodes.txt")
edges <- read_delim("results/networks/model_edges.txt")

# Import nodes coordinates determined by vivagraph
layout <- read.csv("results/networks/layout.csv")

# Mark nodes as GWAS
nodes %<>% 
  mutate(gwas = ifelse(alias %in% gwas_genes, "gwas", "not_gwas"))

# Plot
g <- graph_from_data_frame(edges, nodes, directed = F)

V(g)$x <- layout[,1]
V(g)$y <- layout[,2]
V(g)$size <- 1
V(g)$a <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
V(g)$b <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
V(g)$c <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)

ggraph(g, x = x, y = y) + 
  geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.2) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "gwas"),
    colour = NA,
    n = 5,
    pie_scale = 0.9,
  show.legend = F) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "not_gwas"),
    colour = NA,
    pie_scale = 0.3,
    show.legend = F
  ) +
  #geom_node_text(aes(label = alias), size = 1.7, nudge_x = 2, nudge_y = 4) + 
  #geom_node_label(aes(label = alias)) + 
  scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
  coord_fixed() +
  theme_graph() -> p

# Save 

if(!dir.exists("results/plots_paper/")) {
  dir.create("results/plots_paper")
}

pdf(file = "~/Área de Trabalho/network.pdf", height = 15, width = 15)
print(p)
dev.off()

# Percentage of total genes in the network: 51,24%
n_distinct(nodes$alias) / n_distinct(diff_df$hgnc_symbol)


# Plot PPI netwoks for each brain region and sex --------------------------

# Create igraph object to selected genes
build_graph <- function(int, genes, neighboors = T) {
  
  if (neighboors) {
    int_temp <- int[int$preferredName_A %in% genes | int$preferredName_B %in% genes, ]
  } else {
    int_temp <- int[int$preferredName_A %in% genes & int$preferredName_B %in% genes, ]
  }
  
  # Remove duplicates (create undirected graph)
  idx <- !duplicated(t(apply(int_temp, 1, sort)))
  int_temp <- int_temp[idx,]
  
  # Create graph
  g <- graph_from_edgelist(as.matrix(int_temp[,1:2]), directed = F)
  
  return(g)
}

# Filter genes in differential df by each group
filter_diff <- function(diff_df, genes, gr) {
  diff_df %>% 
    filter(group %in% gr & hgnc_symbol %in% genes) %>% 
    unique() -> diff_df
  
  d_sym <- anyDuplicated(diff_df$hgnc_symbol)
  if (d_sym != 0) {
    message(paste("Dataframe with duplicated genes!", "Check if graph has genes with multiple types", gr))
  }
  return(diff_df)
}

# Create pie chart for nodes
# Say which genes are which types
create_ann_pie <- function(diff_df, genes) {
  
  # Create list with types
  tps <- lapply(split(diff_df$type, diff_df$hgnc_symbol), unique)
  tps <- tps[diff_df$hgnc_symbol]
  
  if(all(names(tps) == diff_df$hgnc_symbol)) {
    message("label types and genes are ok.")
  }
  
  # Create list to pie chart
  ls <- list()
  f_ls <- lapply(seq_along(genes), function(i) {
    if (genes[i] %in% names(tps)) {
      ls[[i]] <- as.integer(c("NONE", "DGE", "DTE", "DTU") %in% tps[names(tps) %in% genes[i]][[1]])
    } else {
      ls[[i]] <- c(1, 0, 0, 0)
    }
  })
  
  return(f_ls)
}

# Set param for igraph object before plotting
set_graph_params <- function(g, dict, f_ls) {
  
  genes <- V(g)$name
  
  # Set node labels
  labels <- sapply(seq_along(genes), function(i) {
    if(genes[i] %in% diff_df$hgnc_symbol) {
      #if(genes[i] %in% dict$hgnc_symbol) {
      return(genes[i])
    } else {
      return(NA)
    }
  })
  
  # Set params
  V(g)$pie.color <- list(c("#666666ff", "#008000ff", "#04009eff", "#ff6600ff"))
  V(g)$shape <- "pie"
  V(g)$pie <- f_ls
  V(g)$size <- 2
  V(g)$label <- labels
  V(g)$label.cex <- 0.6
  V(g)$label.dist <- 0.5
  V(g)$label.color <- "black"
  V(g)$frame.color <- NA
  
  return(g)
  
}

l_groups <- map(split(diff_df$hgnc_symbol, diff_df$group), unique)

graphs_by_group <- imap(l_groups, function(x, i) {
  
  # CHANGE HERE IF YOU WANT FIRST NEIGHBORS
  # (also remember to change the path)
  include_first_neighbors <- F
  
  g <- build_graph(int, x, neighboors = include_first_neighbors)
  genes <- V(g)$name
  
  diff_temp <- filter_diff(diff_df, genes, i)
  f_ls <- create_ann_pie(diff_temp, genes)
  
  # To set all gene labels, set diff_df instead of diff_temp
  g <- set_graph_params(g, dc, f_ls)
  
  if(include_first_neighbors) {
    colors_list <- V(g)$pie
    degrees <- degree(g,v=V(g))
    filter <- !(paste(colors_list) == "c(1, 0, 0, 0)" & degrees == 1)
    g <- induced_subgraph(g, filter)
  }
  
  if(length(V(g)$pie.color) != 0) {
    pdf(paste0("results/networks/", i, ".pdf"), width = 10, height = 10)
    plot(g)
    dev.off()
  }
  
  return(list(graph = g, diff_temp = diff_temp, genes_from_graph = V(g)$name))
  
})

save(graphs_by_group, file = "results/networks/graphs_by_group_wo_nbor.rda")

  
  
  

