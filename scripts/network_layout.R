
# Create network layout with vivagraph ------------------------------------

library(easylayout)

# Organize main layout
layout <- easylayout::vivagraph(g)

# Pin unconnected nodes
layout <- easylayout::vivagraph(g, layout = layout, pin_nodes = TRUE, pinned_cols = 9, pinned_rows = 10)

# Save layout
colnames(layout) <- c("x", "y")

layout %>%
  write.csv("results/networks/layout.csv", quote = F, row.names = F)
