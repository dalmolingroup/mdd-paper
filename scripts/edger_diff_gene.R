# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS

# GENE
library(edgeR)
library(dplyr)
library(purrr)

# Load metadata -----------------------------------------------------------
load("results/important_variables/ann.rda")

# Set comparisons
comp <- paste(rep(unique(ann$region), each = 2), unique(ann$gender), sep = "_")

# Load gene counts --------------------------------------------------------
load("results/txi/txi_gene.rda")

# Differential expression -------------------------------------------------

# Creating DGElist object with design matrix including selected covariates ('rin' and 'ph').
# See https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR
counts <- txi$counts
norm_mat <- txi$length
norm_mat <- norm_mat/exp(rowMeans(log(norm_mat)))
norm_counts <- counts / norm_mat
eff_library_size <- calcNormFactors(norm_counts, method = "TMM") * colSums(norm_counts)
norm_mat <- sweep(norm_mat, 2, eff_library_size, "*")
norm_mat <- log(norm_mat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(counts, group = ann$group)
y <- scaleOffset(y, norm_mat)

# Design matrix
design <- model.matrix(~ 0 + ph + rin + group, data = ann)
colnames(design)[3:ncol(design)] <- gsub("group", "", colnames(design)[3:ncol(design)])

# design <- model.matrix(~ 0 + ph + rin + region:phenotype:gender, data = ann)
# colnames(design)[3:ncol(design)] <- sapply(strsplit(colnames(design)[3:ncol(design)], ":"), function(x) {
#   paste(gsub("region", "", x[1]), gsub("phenotype", "", x[2]), gsub("gender", "", x[3]), sep = "_")
# })

ct <- makeContrasts(
  aINS_female = aINS_MDD_female-aINS_CTRL_female,
  aINS_male = aINS_MDD_male-aINS_CTRL_male,
  Cg25_female = Cg25_MDD_female-Cg25_CTRL_female,
  Cg25_male = Cg25_MDD_male-Cg25_CTRL_male,
  dlPFC_female = dlPFC_MDD_female-dlPFC_CTRL_female,
  dlPFC_male = dlPFC_MDD_male-dlPFC_CTRL_male,
  Nac_female = Nac_MDD_female-Nac_CTRL_female,
  Nac_male = Nac_MDD_male-Nac_CTRL_male,
  OFC_female = OFC_MDD_female-OFC_CTRL_female,
  OFC_male = OFC_MDD_male-OFC_CTRL_male,
  Sub_female = Sub_MDD_female-Sub_CTRL_female,
  Sub_male = Sub_MDD_male-Sub_CTRL_male, 
  levels = design
)

# Filter genes of low expression
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep,]

# TMM normalization
# y <- calcNormFactors(y)

# Estimate dispersions with 'estimateGLMRobustDisp'
y <- estimateGLMRobustDisp(y, design = design, verbose = T)

# Fit model
fit <- glmFit(y, design = design)

# Extract results for each comparison

# Summary dataframe
map_df(comp, function(c) {
  
  cat("Comparison: ", c, )
  
  lrt <- glmLRT(fit, contrast = ct[, c])
  df <- topTags(lrt, n = Inf, adjust.method = "BH", p.value = 0.05)$table
  if(!is.null(df)) {
    df %>%
      mutate(group = c,
             gene = rownames(df)) %>% 
      relocate(gene)
  }
}) -> df_edger_ph_rin_group_gene

# LRT DGElist objects for each comparison
map(comp, function(c) {
  lrt <- glmLRT(fit, contrast = ct[, c])
}) -> lrt_comp
names(lrt_comp) <- comp

if(!dir.exists("results/diff_exp/")) {
  dir.create("results/diff_exp/")
}

save(df_edger_ph_rin_group_gene, lrt_comp, file = "results/diff_exp/edger_gene_rin_ph_diff.rda")

