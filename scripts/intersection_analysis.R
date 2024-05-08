
# Quantitative analysis ---------------------------------------------------

library(DESeq2)
library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(UpSetR)
library(gplots)

# Load results from DGE, DTE and DTU --------------------------------------
load("results/diff_exp/diff_df.rda")

diff_df %<>% 
  separate(col = "group", into = c("region", "sex"), sep = "_", remove = F)

# Intersections -----------------------------------------------------------
# Create intersections between all possible comparisons

# By analysis type
pdf("results/intersects/genes_by_type.pdf", height = 5, width = 10)
l_type <- lapply(split(diff_df$gene, diff_df$type), unique)
upset(fromList(l_type), text.scale = c(1.5, 2, 1, 1.5, 2, 3), point.size = 5, nintersects = NA)
dev.off()

png("results/intersects/genes_by_type.png",  width = 20, height = 14, units = "cm", res = 300)
l_type <- lapply(split(diff_df$gene, diff_df$type), unique)
upset(fromList(l_type), text.scale = c(1.5, 2, 1, 1.5, 2, 3), point.size = 5, nintersects = NA)
dev.off()

# By groups in each sex 

# Female
pdf("results/intersects/genes_by_group_female.pdf", height = 7, width = 15)
l_group_female <- lapply(split(diff_df$gene[diff_df$sex == "female"], diff_df$group[diff_df$sex == "female"]), unique)
upset(fromList(l_group_female), nsets = length(l_group_female),  text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

png("results/intersects/genes_by_group_female.png",  width = 20, height = 14, units = "cm", res = 300)
upset(fromList(l_group_female), nsets = length(l_group_female),  text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

# Male
pdf("results/intersects/genes_by_group_male.pdf", height = 7, width = 15)
l_group_male <- lapply(split(diff_df$gene[diff_df$sex == "male"], diff_df$group[diff_df$sex == "male"]), unique)
upset(fromList(l_group_male), nsets = length(l_group_male), text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

png("results/intersects/genes_by_group_male.png",  width = 20, height = 14, units = "cm", res = 300)
upset(fromList(l_group_male), nsets = length(l_group_male),  text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

# By regions
pdf("results/intersects/genes_by_regions.pdf", height = 7, width = 15)
l_regions <- lapply(split(diff_df$gene, diff_df$region), unique)
upset(fromList(l_regions), nsets = length(l_regions), text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

png("results/intersects/genes_by_region.png",  width = 20, height = 14, units = "cm", res = 300)
upset(fromList(l_regions), nsets = length(l_regions),  text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

# By sex
pdf("results/intersects/genes_by_sex.pdf", height = 5, width = 10)
l_sex <- lapply(split(diff_df$gene, diff_df$sex), unique)
upset(fromList(l_sex), nsets = length(l_sex), text.scale = c(1.5, 2, 1, 1.5, 2, 3))
dev.off()

png("results/intersects/genes_by_sex.png",  width = 20, height = 14, units = "cm", res = 300)
upset(fromList(l_sex), nsets = length(l_sex),  text.scale = c(1.5, 1.5, 1, 1, 2, 2), 
      point.size = 3, nintersects = NA)
dev.off()

# Get the list of genes for each intersection
get_intersect_genes <- function(ls) {
  tmp <- gplots::venn(ls, show.plot = F)
  tmp <- attr(tmp, "intersections")
  df <- data.frame(intersect = rep(names(tmp), sapply(tmp, length)), genes = unlist(tmp))
  rownames(df) <- NULL
  return(df)
}
genes_by_sex <- get_intersect_genes(l_sex)
genes_by_type <- get_intersect_genes(l_type)
genes_by_regions <- get_intersect_genes(l_regions)
genes_by_group_female <- get_intersect_genes(l_group_female)
genes_by_group_male <- get_intersect_genes(l_group_male)

if(!dir.exists("results/intersects")) {
  dir.create("results/intersects", recursive = T)
}
save(list = ls()[grep("genes_by", ls())], file = "results/intersects/intersects.rda")

# Save csv files for each comparison
lapply(c(
  "genes_by_sex",
  "genes_by_type",
  "genes_by_regions",
  "genes_by_group_female",
  "genes_by_group_male"
  ), 
function(x) {
  write.csv(get(x), file = paste0("results/intersects/", x, ".csv"), row.names = F, quote = F)
})

# Number of genes ---------------------------------------------------------
# Plot number of genes in each female group
diff_df %>% 
  filter(sex == "female") %>% 
  group_by(region, type) %>% 
  dplyr::count(name = "n") %>% 
ggplot(aes(x = type, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ region) + 
  scale_y_continuous(name = "Number of genes", limits = c(0, 250)) +
  scale_fill_manual(name = "", values = c("DGE" = "#008000ff", "DTE" = "#04009eff", "DTU" = "#ff6600ff")) +
  scale_x_discrete(name = "") +
  ggtitle("Expression analysis for female regions") +
  guides(fill = F) +
  theme_bw()
ggsave("results/intersects/number_of_genes_by_sex_female.pdf", width = 7, height = 3)
ggsave("results/intersects/number_of_genes_by_sex_female.png", width = 7, height = 3)

# Plot number of genes in each male group
diff_df %>% 
  filter(sex == "male") %>% 
  group_by(region, type) %>% 
  dplyr::count(name = "n") %>% 
  ggplot(aes(x = type, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ region) + 
  scale_y_continuous(name = "Number of genes", limits = c(0, 250)) +
  scale_fill_manual(name = "", values = c("DGE" = "#008000ff", "DTE" = "#04009eff", "DTU" = "#ff6600ff")) +
  scale_x_discrete(name = "") +
  ggtitle("Expression analysis for male regions") +
  guides(fill = F) +
  theme_bw()
ggsave("results/intersects/number_of_genes_by_sex_male.pdf", width = 7, height = 3)
ggsave("results/intersects/number_of_genes_by_sex_male.png", width = 7, height = 3)

# Gather the number of significant genes in male and female for each group
diff_df %>% 
  group_by(sex, region, type) %>% 
  dplyr::count(name = "n") %>% 
  ggplot(aes(x = type, y = n, fill = type)) +
    geom_bar(stat = "identity") +
    facet_grid(sex ~ region, scales = "free_y") + 
    scale_y_continuous(name = "Number of genes", limits = c(0, 250)) +
    scale_fill_manual(name = "", values = c("DGE" = "#008000ff", "DTE" = "#04009eff", "DTU" = "#ff6600ff")) +
    scale_x_discrete(name = "") +
    guides(fill = F) +
    theme_bw()
ggsave("results/intersects/number_of_genes_by_sex.pdf", width = 9, height = 5)
ggsave("results/intersects/number_of_genes_by_sex.png", width = 9, height = 5)

# Raw number of altered genes in each analysis
diff_df %>% 
  dplyr::select(gene, type) %>% 
  group_by(type) %>% 
  unique() %>% 
  dplyr::count(name = "n") %>% 
  ggplot(aes(x = type, y = n, fill = type)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "Number of genes", limits = c(0, 700)) +
    scale_fill_manual(name = "", values = c("DGE" = "#008000ff", "DTE" = "#04009eff", "DTU" = "#ff6600ff")) +
    scale_x_discrete(name = "") + 
    guides(fill = F) +
    ggtitle("Raw number of genes for each analysis") +
    theme_bw()
ggsave("results/intersects/number_of_genes_by_type.pdf", width = 6, height = 3)
ggsave("results/intersects/number_of_genes_by_type.png", width = 6, height = 3)





    
  
  


