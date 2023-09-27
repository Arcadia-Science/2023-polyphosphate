library(tidyverse)

#################################################
# Filter Uniprot accessions and create metadata table
#################################################

# bacteria ppk1 uniprot accessions 
ppk1_uniprot_bacteria_accessions <- read.delim("metadata/uniprot-ppk1-bacteria-accessions.tsv", sep="\t", header=TRUE, na.strings = "")

ppk1_uniprot_bacteria_accessions %>% 
  ggplot(aes(x=Length)) +
  geom_histogram() # majority with length of ~700 AAs, filter above 500 to get rid of trash

ppk1_uniprot_bacteria_accessions_filtered <- ppk1_uniprot_bacteria_accessions %>% 
  filter(Length > 500) %>% 
  filter(!is.na(AlphaFoldDB)) %>% 
  select(Entry, Protein.names, Length, Organism, Taxonomic.lineage) %>% 
  mutate(Phylum = str_extract(Taxonomic.lineage, "(?<=\\(class\\),\\s)(.*?)(?=\\s\\(phylum\\))"))

# archaea ppk1 accessions 
ppk1_uniprot_archaea_accessions <- read.delim("metadata/uniprot-ppk1-archaea-accessions.tsv", sep="\t", header=TRUE, na.strings="")

ppk1_uniprot_archaea_accessions %>% 
  ggplot(aes(x=Length)) + 
  geom_histogram() # same length distribution as well, filter above 500 AAs

ppk1_uniprot_archaea_accessions_filtered <- ppk1_uniprot_archaea_accessions %>% 
  filter(Length > 500) %>% 
  filter(!is.na(AlphaFoldDB)) %>% 
  select(Entry, Protein.names, Length, Organism, Taxonomic.lineage) %>% 
  mutate(Phylum = str_extract(Taxonomic.lineage, "(?<=\\(class\\),\\s)(.*?)(?=\\s\\(phylum\\))"))

ppk1_uniprot_archaea_accessions_filtered$Phylum <- gsub("Stenosarchaea group \\(no rank\\), ", "", ppk1_uniprot_archaea_accessions_filtered$Phylum)

ppk1_uniprot_archaea_accessions_filtered$Phylum <- gsub("Methanomada group \\(no rank\\), ", "", ppk1_uniprot_archaea_accessions_filtered$Phylum)


# combine bacteria and archaea accessions into one metadata file
all_filtered_ppk1_accessions <- rbind(ppk1_uniprot_archaea_accessions_filtered, ppk1_uniprot_bacteria_accessions_filtered)

write.table(all_filtered_ppk1_accessions, "metadata/all-filtered-ppk1-accessions.tsv", sep="\t", row.names = FALSE, quote = FALSE)

# save accessions in a list for downloading protein FASTAs and AlphafoldDB structures
filtered_ppk1_accessions <- all_filtered_ppk1_accessions %>% 
  select(Entry)

write.table(filtered_ppk1_accessions, "metadata/filtered-ppk1-accessions-list.txt", quote = FALSE, row.names = FALSE)
