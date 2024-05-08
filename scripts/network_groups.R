
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(ggraph)


# Create individual networks ----------------------------------------------

load("results/networks/graphs_by_group_wo_nbor.rda")
load("results/diff_exp/diff_df.rda")
load("results/networks/int.rda")

make_graph_by_group <- function(genes, group_name, diff_df, int, nx = 0.3, ny = 0.3, ps = 0.5) {
  
  if(!is.null(genes)) {
    
    diff_df_temp <- diff_df %>% 
      dplyr::select(hgnc_symbol, type, group) %>% 
      dplyr::filter(hgnc_symbol %in% genes, group %in% group_name) %>% 
      distinct()
    
    int_temp <- int %>% 
      filter(preferredName_A %in% diff_df_temp$hgnc_symbol, preferredName_B %in% diff_df_temp$hgnc_symbol)
    
    g <- graph_from_edgelist(as.matrix(int_temp[,c(1,2)]))
    
    # lay <- create_layout(g, layout = "nicely")
    # V(g)$x <- lay[,"x"]
    # V(g)$y <- lay[,"y"]
    V(g)$a <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
    V(g)$b <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
    V(g)$c <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)
    
    layout <- easylayout::vivagraph(graphs_by_group[[group_name]]$graph, pin_nodes = FALSE)
    layout <- easylayout::vivagraph(graphs_by_group[[group_name]]$graph, layout = layout, pin_nodes = TRUE, pinned_cols = 10, lcc_margin_left = 500)
    V(g)$x <- layout[, 1]
    V(g)$y <- layout[, 2]
    
    # buffer <- c(-2, 2)
    # xlims <- ceiling(range(V(g)$x)) + buffer
    # ylims <- ceiling(range(V(g)$y)) + buffer
    
    #xlims <- ylims <- c(-50, 50)
    
    ggraph(g, x = x, y = y) +
      geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.2) +
      geom_scatterpie(
        data = as_data_frame(g, "vertices"),
        cols = c("a", "b", "c"),
        color = NA,
        pie_scale = ps,
        show.legend = F) +
      geom_node_text(aes(label = name), size = 1,  nudge_x = nx, nudge_y = ny) + 
      scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
      coord_fixed() +
      theme_void() +
      theme(plot.margin = margin(3,1,1.5,3, "cm")) -> plot_g
    
    return(plot_g)
    
  } else {
    
    return(NULL)
    
  }
  
}

imap(graphs_by_group, ~ {
  
make_graph_by_group(.x$genes_from_graph, .y, diff_df, int)
  
}) %>% 
  discard(is.null) -> ls_graphs


design_female <- "
AB
CD
EE
EE
"

wrap_plots(ls_graphs[c(1,2,4,5,7)], design = design_female)
ggsave("results/networks/ppi_female.pdf", height = 10, width = 10)


# -------------------------------------------------------------------------

make_graph_by_group(graphs_by_group$aINS_female$genes_from_graph, 
                    group_name = "aINS_female",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.1, ny = 0.1)
ggsave(filename = "results/networks/aINS_female.pdf", width = 5, height = 5)


make_graph_by_group(graphs_by_group$Cg25_female$genes_from_graph, 
                    group_name = "Cg25_female",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.1, ny = 0.1)
ggsave(filename = "results/networks/Cg25_female.pdf", width = 5, height = 5)

make_graph_by_group(graphs_by_group$dlPFC_female$genes_from_graph, 
                    group_name = "dlPFC_female",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.1, ny = 0.1)
ggsave(filename = "results/networks/dlPFC_female.pdf", width = 5, height = 5)

make_graph_by_group(graphs_by_group$Nac_female$genes_from_graph, 
                    group_name = "Nac_female",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.2, ny = 0.2)
ggsave(filename = "results/networks/Nac_female.pdf", width = 5, height = 5)

make_graph_by_group(graphs_by_group$OFC_female$genes_from_graph, 
                    group_name = "OFC_female",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.2, ny = 0.2, ps = 0.3)
ggsave(filename = "results/networks/OFC_female.pdf", width = 14, height = 14)

###
make_graph_by_group(graphs_by_group$Nac_male$genes_from_graph, 
                    group_name = "Nac_male",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.2, ny = 0.2)
ggsave(filename = "results/networks/Nac_male.pdf", width = 5, height = 5)

make_graph_by_group(graphs_by_group$Cg25_male$genes_from_graph, 
                    group_name = "Cg25_male",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.2, ny = 0.2, ps = 0.3)
ggsave(filename = "results/networks/Cg25_male.pdf", width = 10, height = 10)


make_graph_by_group(graphs_by_group$Sub_male$genes_from_graph, 
                    group_name = "Sub_male",
                    diff_df = diff_df,
                    int = int,
                    nx = 0.2, ny = 0.2)
ggsave(filename = "results/networks/Sub_male.pdf", width = 7, height = 7)
