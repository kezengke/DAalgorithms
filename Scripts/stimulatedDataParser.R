rm(list = ls())
library(tidyr)
library(dplyr)

# Assuming your data is in a data frame called `df`
df <- read.table("mouse_transformed_simulated_counts.txt", sep = "\t", header = T)

# Combine parent_lineage and lineage
df <- df %>%
  mutate(species_taxon = paste0(parent_lineage, "_", lineage))

# Reshape data to wide format
otu_table <- df %>%
  select(sample_id, species_taxon, taxa_count_raw) %>%
  pivot_wider(names_from = sample_id, values_from = taxa_count_raw, values_fill = list(taxa_count_raw = 0))

# Set row names as species_taxon
otu_table <- as.data.frame(otu_table)
rownames(otu_table) <- otu_table$species_taxon
otu_table$species_taxon <- NULL

write.table(otu_table, "CountsTables/WGSRaw/mouse_simulated.txt", sep = "\t", row.names = T)
