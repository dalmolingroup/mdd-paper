
# Read transcript counts --------------------------------------------------

library(tximport)

# Load annotation 
load("results/important_variables/ann.rda")

# Load counts from kallisto -----------------------------------------------

# Read files
files <- list.files(path = "data/kallisto/", pattern="tsv", recursive = TRUE, full.names = TRUE)
files <- files[sapply(rownames(ann), function(x) grep(x, files))]

# Samples to remove 

names(files) <- rownames(ann)
txi <- tximport(files = files, type = "kallisto", txOut = T)

# Save
if(!dir.exists("results/txi")) {
  dir.create("results/txi")
}
save(txi, file = "results/txi/txi_tx.rda")
