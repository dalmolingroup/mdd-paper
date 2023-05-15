
# Variables correlation study ---------------------------------------------

library(tidyverse)
library(mice)
library(lattice)
library(directlabels)
library(PCAtools)
library(magrittr)

# Read and organize table 
ann <- read.table("data/meta/SraRunTable.txt", header = T, stringsAsFactors = F, sep = ",")
colnames(ann) <- tolower(colnames(ann))
ann <- ann %>% 
  dplyr::select(run, age, alcool, drugs, drug_type, gender, medication, medication_type, 
                organism, smoking, ph, phenotype, pmi, rin, tissue) %>% 
  filter(organism == "Homo sapiens") %>% 
  mutate(region = case_when(
    tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
    tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
    tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
    tissue == "Anterior Insula (aINS)" ~ "aINS",
    tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
    tissue == "Subiculum (Sub)" ~ "Sub"
  ), 
  group = paste(region, phenotype, gender, sep = "_")) %>% 
  dplyr::rename(sample_id = run)

# Start analysis of important variables -----------------------------------

# Transform categorical data into numeric for imputation steps
invisible(lapply(seq_along(ann), function(i) {
  if(colnames(ann)[i] %in% c("alcool", "drugs", "drug_type", "medication", "medication_type", "smoking")) {
    ann[,i] <<- as.character(tolower(ann[,i]))
    ann[,i] <<- as.numeric(as.factor(ann[,i]))
  }
}))

# Check how many NA values are there in each category
sapply(ann, function(i) sum(is.na(i)) )
ann_miss <- VIM::aggr(ann, col = mdc(1:2), numbers = TRUE, sortVars = TRUE, 
                      labels = names(ann), cex.axis = 0.7, gap = 3, 
                      ylab = c("Proportion of missingness", "Missingness Pattern"))

# Predict imput values (5 predictions)
ann_imputed <- mice(ann, m = 5, maxit = 40, method = "pmm", seed = 500)

# Check predicted values
xyplot(ann_imputed, alcool ~ drugs | .imp, pch = 20, cex = 1.4, alpha = 0.4)
direct.label(densityplot(x = ann_imputed, data = ~ alcool + drugs + medication + medication_type))

# Choose dataset 2 (one can choose any)
ann_complete <- mice::complete(ann_imputed, 2)

# Remove column 
ann_complete <- ann_complete %>% 
  dplyr::select(-tissue, -drug_type, -medication_type)

# Check the number of missing values 
sum(sapply(ann_complete, function(i) sum(is.na(i)) )) == 0

# Scale continuos variables
ann_complete$age <- scale(ann_complete$age)[,1]
ann_complete$pmi <- scale(ann_complete$pmi)[,1]
ann_complete$rin <- scale(ann_complete$rin)[,1]
ann_complete$ph <- scale(ann_complete$ph)[,1]

# Transform categorical data
ann_complete$gender <- as.numeric(as.factor(ann_complete$gender))
ann_complete$region <- as.numeric(as.factor(ann_complete$region))
ann_complete$phenotype <- as.numeric(as.factor(ann_complete$phenotype))

# Remove sample names
rownames(ann_complete) <- ann_complete$sample_id
ann_complete$sample_id <- NULL
ann_complete %<>% 
  tibble::rownames_to_column("run") 

# Remove outliers
ann_complete %<>%
  filter(!(run %in% c("SRR5961961", "SRR5961809")))

# Save
if(!dir.exists("results/important_variables")) {
  dir.create("results/important_variables/", recursive = T)
}

save(ann_complete, file = "results/important_variables/ann_complete.rda")



