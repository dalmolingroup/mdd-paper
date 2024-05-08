library(tidyverse)
library(PCAtools)
library(DESeq2)
library(pals)
library(pheatmap)
library(RColorBrewer)
# library(animation)
# library(gganimate)
# library(magick)
# library(ExpressionNormalizationWorkflow)
# library(PoiClaClu)

# Load predicted metadata
load("results/important_variables/ann_complete.rda")

ann_complete <- ann_complete %>% 
  dplyr::rename(sex = gender)

# Organize predicted metadata
ann_regression <- ann_complete
ann_pca <- ann_complete
lapply(seq_along(ann_pca), function(i) {
  if(colnames(ann_pca)[i] %in% c("alcool", "drugs", "medication", "smoking")) {
    ann_pca[,i] <<- factor(ann_pca[,i], levels = c("1", "2"), labels = c("no", "yes"))
  }
})
ann_pca$phenotype <- factor(ann_pca$phenotype, levels = c("1", "2"), labels = c("CTRL", "MDD"))
ann_pca$region <- sapply(strsplit(ann_pca$group, split = "_"), "[[", 1)
ann_pca$sex <- factor(ann_pca$sex, levels = c("1", "2"), labels = c("female", "male"))


rm("ann_complete")

# Load counts from kallisto
load("results/txi/txi_gene.rda")

# Exploratory analysis ----------------------------------------------------

# DESeq2 object
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = ann_pca, 
                                design = ~ group)

# Filter genes that do not have a count value lesser than 10 in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# VST normalization - See DESeq workflow for details
vds <- vst(dds, blind = T)

# Plot PCA individually based on the 500 genes with highest variance 
get_plots <- function(n_genes) {
  vsd_data <- DESeq2::plotPCA(vds, intgroup = "group", returnData = T, ntop = n_genes)
  vsd_data$n <- n_genes
  return(vsd_data)
}

# # Plot animation 
# df <- do.call(rbind, lapply(c(500, 2000, 5000, 10000, 20000, nrow(vds)),  get_plots))
# df$gender <- sapply(strsplit(as.character(df$group), split =  "_"), "[[", 3)
# df$group <- paste(ann_pca$region, ann_pca$phenotype, sep = "_")
# df$n <- as.factor(df$n)
# 
# anim1 <- ggplot(df, aes(x = PC1, y = PC2, col = group, shape = gender)) +
#   geom_point(size = 2.5) +
#   transition_states(n, transition_length = 3, state_length = 1)  + 
#   labs(title = "PCA with the {closest_state} genes with highest variance")
# 
# anim2 <- animate(anim1, width = 800, height = 500)
# 
# if(!dir.exists("results/rank_plots/")) {
#   dir.create("results/rank_plots/")
# }
# anim_save(filename = "results/rank_plots/pca_high_var_genes.gif", animation = anim2)

# PCA analysis of variants ------------------------------------------------

# Compute and plot PCA
rownames(ann_regression) <- ann_regression$run
p_kallisto <- pca(assay(vds), metadata = ann_regression)

df <- data.frame(PC1 = p_kallisto$rotated[,1],
                 PC2 = p_kallisto$rotated[,2],
                 sex = ann_pca$sex,
                 group = paste(ann_pca$region, ann_pca$phenotype, sep = "_"),
                 stringsAsFactors = F)
cols <- brewer.paired(12)
cols[11] <- "#e0e05cf8"
names(cols) <- unique(df$group)[order(sapply(strsplit(unique(df$group), split = "_"), "[[", 1))]
ggplot(df, aes(x = PC1, y = PC2, col = group, shape = sex)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols, name = "Group") +
  scale_shape_manual(name = "Sex", values = c(15, 17), labels = c("female" = "Female", "male" = "Male")) +
  labs(x = paste0("PC1", " (", round(p_kallisto$variance[1], 2), "%)"),
       y = paste0("PC2", " (", round(p_kallisto$variance[2], 2), "%)")) +
  guides(colour = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) + 
  theme_bw() + 
  theme(legend.key.size = unit(10, units = "mm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

if(!dir.exists("results/rank_plots/")) {
  dir.create("results/rank_plots/")
}

ggsave(plot = last_plot(), "results/rank_plots/pca_vst.pdf", width = 10, height = 7)
ggsave(plot = last_plot(), "results/rank_plots/pca_vst.png", width = 10, height = 7)

# Screeplot
sc_plot <- screeplot(p_kallisto, vline = findElbowPoint(p_kallisto$variance), components = 1:20)
ggsave(plot = sc_plot, file = "results/rank_plots/screeplot.pdf", width = 10, height = 7)
ggsave(plot = sc_plot, file = "results/rank_plots/screeplot.png", width = 10, height = 7)

# Eigenplot - spearman correlation between variables and principal components
pdf(file = "results/rank_plots/eigencorplot.pdf", width = 10, height = 7)
eigencorplot(p_kallisto, corFUN = "spearman",
                         corUSE = 'pairwise.complete.obs',
                         metavars = c("age", "alcool", "drugs",
                                      "ph", "medication", "sex", 
                                      "region", "pmi", "rin", "phenotype"))
dev.off()

png(file = "results/rank_plots/eigencorplot.png", width = 20, height = 14, units = "cm", res = 300)
eigencorplot(p_kallisto, corFUN = "spearman",
             corUSE = 'pairwise.complete.obs',
             metavars = c("age", "alcool", "drugs",
                          "ph", "medication", "sex", 
                          "region", "pmi", "rin", "phenotype"))
dev.off()

# # Poisson distance - See similarities between regions
# poisd <- PoissonDistance(t(counts(dds)))
# samplePoisDistMatrix <- as.matrix(poisd$dd)
# rownames(samplePoisDistMatrix) <- ann_pca$group
# colnames(samplePoisDistMatrix) <- NULL
# colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# 
# pdf("results/rank_plots/poisson_dist.pdf", height = 20, width = 10)
# pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = colors, 
#          fontsize_row = 5,
#          cellwidth = 1, cellheight = 5, legend = F)
# dev.off()

# # Principal Variance Component Analysis -----------------------------------
# # From ExpressionNormalizationWorkflow package
# rownames(ann_pca) <- ann_pca$run
# inpData <- expSetobj(assay(vds), ann_pca)
# cvrts_eff_var <- c("age", "alcool", "drugs", "ph", "medication", 
#                    "gender", "region", "pmi", "rin", "phenotype")
# pct_thrsh <- 0.95
# pdf("results/rank_plots/pvca.pdf", width = 10, height = 7)
# pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)
# dev.off()

# Rank variables by their correlation to each components ------------------
vars <- c("age", "alcool", "drugs", "ph", "medication", "smoking",
          "sex", "region", "pmi", "rin", "phenotype")
m_cor <- matrix(0, nrow = length(vars), ncol = length(p_kallisto$rotated))

# Create correlation matrix
for (i in 1:nrow(m_cor)) {
  for (j in 1:ncol(m_cor)) {
    m_cor[i, j] <- cor.test(x = ann_regression[, vars[i] ], 
                            y = p_kallisto$rotated[,j],
                            method = "spearman")$p.value
  }
}
m_cor <- m_cor[,1:10]
dimnames(m_cor) <- list(vars, paste0("PC", 1:10))

# Rank variables by PCs
ls <- list()
for (i in 1:nrow(m_cor)) {
  temp <- 0
  for (j in 1:ncol(m_cor)) {
    if (m_cor[i,j] <= 0.05) {
      temp[j] <-  p_kallisto$variance[j]
      ls[[i]] <- temp
    } else {
      temp[j] <-  NA
      ls[[i]] <- temp
    }
  }
}
names(ls) <- vars
df_vars <- as.data.frame(do.call(rbind, ls))
df_vars$vars <- vars
names(df_vars)[1:10] <- paste0("PC", 1:10)
df_vars <- gather(df_vars, key = "component", value = "value", PC1:PC10)
df_vars <- df_vars[complete.cases(df_vars),]
df_vars$vars <- factor(df_vars$vars, levels = names(sort(sapply(ls, sum, na.rm = T), decreasing = T)),
                       labels = names(sort(sapply(ls, sum, na.rm = T), decreasing = T)))
df_vars$component <- as.numeric(gsub("PC", "", df_vars$component))

# Plot cumulative importance for each significantly correlated variable
ggplot(df_vars, aes(x = vars, y = value, fill = factor(component))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = unname(polychrome())[3:12], labels = as.factor(1:10)) +
  labs(x = "", y = "Cumulative variance of correlated components %", fill = "PC") +
  ggtitle("Variance explained by the cumulative sum of components correlated to each variable",
          subtitle = paste0("Variance explained by the first 10 components = ", 
                            round(sum(p_kallisto$variance[1:10]), 2), "%")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14))
ggsave("results/rank_plots/rank_vars.pdf", width = 10, height = 7)
ggsave("results/rank_plots/rank_vars.png", width = 10, height = 7, dpi = 300)

