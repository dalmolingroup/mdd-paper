
# Relate patients to their respective samples -----------------------------

library(dplyr)
library(tidyr)
library(janitor)

# Get organized metadata about each patient -------------------------------
# This metadata comes from the supplementary material that can be found here:
# http://neuroscience.mssm.edu/nestler/contecenter/Cohort%201%20Metadata.pdf
patients <- read.csv("~/Downloads/Metadata.csv")

# Organize data for further merge
colnames(patients) <- make_clean_names(colnames(patients))
patients <- patients %>% 
  pivot_longer(cg25:sub, names_to = "region", values_to = "rin") %>% 
  select(patients = samples, gender, region, group, age, pmi, ph, rin) %>% 
  mutate(across(where(is.numeric), ~ as.character(.)),
         gender = tolower(gender))

# Get samples metadata ----------------------------------------------------

ann <- read.csv("data/meta/SraRunTable.txt", stringsAsFactors = F, header = T)
colnames(ann) <- tolower(colnames(ann))
ann <- ann %>% 
  dplyr::select(run, ph, rin, pmi, age, phenotype, gender, tissue, organism) %>% 
  filter(organism == "Homo sapiens") %>% 
  mutate(region = case_when(
    tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
    tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
    tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
    tissue == "Anterior Insula (aINS)" ~ "aINS",
    tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
    tissue == "Subiculum (Sub)" ~ "Sub"
  ),
  region = tolower(region)) %>% 
  dplyr::rename(sample_id = run) %>% 
  dplyr::select(sample_id, ph, rin, pmi, age, group = phenotype, gender, region) %>% 
  mutate(across(where(is.numeric), ~ as.character(.)))


patients_metadata_table <- inner_join(ann, patients, by = c("gender", "group", "region", "ph", "rin", "pmi", "age"))

write.csv(patients_metadata_table, file = "results/tables/patients_metadata.csv", quote = F, row.names = F)
