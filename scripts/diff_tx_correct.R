# Perform gene-wise pvalue correction with stageR -------------------------

library(dplyr)
library(stageR)
library(GenomicFeatures)

# Load TX data from differential expression
load("results/diff_exp/edger_tx_rin_ph_diff.rda")
df_res <- df_edger_ph_rin_group_tx
colnames(df_res)[1] <- "tx"

# remove transcript version
df_res$tx <- gsub("\\.+\\d+", "", rownames(df_res))

# Load transcript-gene info -----------------------------------------------
gtf <- "data/genome/Homo_sapiens.GRCh38.97.gtf.gz"
txdb.filename <- "data/genome/Homo_sapiens.GRCh38.97.gtf.sqlite"
#gtf <- "Homo_sapiens.GRCh38.97.gtf.gz"
#txdb.filename <- "Homo_sapiens.GRCh38.97.gtf.sqlite"

# Load db
txdb <- loadDb(txdb.filename)
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
tx2gene <- data.frame(tx = txdf$TXNAME, gene = txdf$GENEID, stringsAsFactors = F)

# Filter by the transcripts annotated in the gtf
df_res1 <- df_res[df_res$tx %in% tx2gene$tx,]

# Set regions and gender variables 
regions <- unique(sapply(strsplit(df_edger_ph_rin_group_tx$group, split = "_"), "[[", 1))
gender <- unique(sapply(strsplit(df_edger_ph_rin_group_tx$group, split = "_"), "[[", 2))

# Correct p-values for each contrast
ls_res <- list()
for (i in 1:length(regions)) {
  ls_temp <- list()
  for (j in 1:length(gender)) {
    
    # Create a subset for each region/gender
    tmp <- df_res1[df_res1$group %in% paste(regions[i], gender[j], sep = "_"),]
    tx2gene_tmp <- unique(merge(tmp, tx2gene, by = "tx")[, c("tx", "gene")])
    tx2gene_tmp <- tx2gene_tmp[order(tx2gene_tmp$gene),]
    
    tmp <- tmp[match(tx2gene_tmp$tx, tmp$tx),]
    
    stopifnot(identical(tmp$tx, tx2gene_tmp$tx))
      
    # Check if all transcripts are in the annotation df
    if (all(tx2gene_tmp$tx %in% tmp$tx)) {
      print("All ok!")
    } else {
      print("Check annotation")
    }
    
    # Perform stageR gene-wise p-value correction
    # pScreen and pConfirmation with padj
    pScreen <- tmp$FDR
    names(pScreen) <- tx2gene_tmp$gene
    pConfirmationTx <- matrix(tmp$PValue, ncol = 1)
    rownames(pConfirmationTx) <- tmp$tx
    
    # Create stageR object
    stageRObj <- stageRTx(pScreen = pScreen, 
                          pConfirmation = pConfirmationTx, 
                          tx2gene = data.frame(transcripts = tx2gene_tmp$tx, genes = tx2gene_tmp$gene, stringsAsFactors = F), 
                          pScreenAdjusted = T)
    
    # Correct p-values with gene-wise approach  
    stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dte", alpha = 0.05)
    
    # Get the corrected values
    padj <- getAdjustedPValues(stageRObj, order = TRUE, onlySignificantGenes = T)
    # padj <- padj[!padj$transcript == 0,]
    
    if (nrow(padj) == 0) {
      ls_temp[[j]] <- NULL
    } else {
      padj$group <- paste(regions[i], gender[j], sep = "_")
      ls_temp[[j]] <- padj
    }
  }
  
  df_temp <- do.call(rbind, ls_temp)
  ls_res[[i]] <- df_temp
}

df_res_padj_tx <- do.call(rbind, ls_res)

if(!dir.exists("results/diff_exp/")) {
  dir.create("results/diff_exp/")
}

# Save results
save(df_res_padj_tx, file = "results/diff_exp/diff_tx_corrected.rda")



