
# Create network layout with vivagraph ------------------------------------

library(easylayout)

# Read igraph data
load("results/networks/int.rda")

# Create graph
g <- graph_from_edgelist(as.matrix(int[,1:2]), directed = F)

# Organize main layout
layout <- easylayout::vivagraph(g)

# Pin unconnected nodes
layout <- easylayout::vivagraph(g, layout = layout, pin_nodes = TRUE, pinned_cols = 9, pinned_rows = 10)

# Save layout
colnames(layout) <- c("x", "y")

layout %>%
  write.csv("results/networks/layout.csv", quote = F, row.names = F)
