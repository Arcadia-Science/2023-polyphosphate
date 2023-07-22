library(tidyverse)
library(ArcadiaColorBrewer)
library(ggpubr)

# metadata
ppk1_metadata <- all_filtered_ppk1_accessions %>% 
  mutate(accession = Entry) %>% 
  select(accession, Taxonomic.lineage, Organism, Phylum)

#################################################
# Accumulibacter phosphatis ppk1 comparisons
#################################################

# mmseqs easy-search results 
ppk1_mmseqs <- read.table("results/CAP_ppk1_sequences_search.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))

# foldseek easy-search results
ppk1_foldseek <- read.table("results/CAP_ppk1_structures_search.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

# histogram of frequencies of seq id
ppk1_mmseqs %>% 
  ggplot(aes(x=seqid)) +
  geom_histogram()

# histogram of frequencies of alntmscores
ppk1_foldseek %>% 
  ggplot(aes(x=alntmscore)) +
  geom_histogram()

# reorganize and combine by accession
ppk1_mmseqs_filtered <- ppk1_mmseqs %>% 
  select(mmseqs_target, seqid) %>% 
  mutate(accession = mmseqs_target) %>% 
  select(accession, seqid)

ppk1_foldseek_filtered <-ppk1_foldseek %>% 
  select(foldseek_target, alntmscore) %>% 
  mutate(accession = gsub("-F1-model_v1.pdb", "", foldseek_target)) %>% 
  select(accession, alntmscore)

ppk1_results <- left_join(ppk1_mmseqs_filtered, ppk1_foldseek_filtered)

# merge with metadata to get taxonomy information
ppk1_results_metadata <- left_join(ppk1_results, ppk1_metadata) %>% 
  mutate(Phylum = if_else(is.na(Phylum), "Other", Phylum))

# plot comparison of protein seqid to alntmscore against Candidatus Accumulibacter ppk1 query 
top_filtered_phyla <- ppk1_results_metadata %>% 
  count(Phylum) %>% 
  top_n(10, n) %>% 
  pull(Phylum)

all_CAP_ppk1_comps_plot <- ppk1_results_metadata %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>% 
  ggplot(aes(x=seqid, y=alntmscore)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr(legend = c(0.8, 0.5)) +
  labs(x="Protein Sequence Identity", y="Protein Strucutre Alignment (Tm score)") +
  ggtitle("Comparisons of Protein Sequence Identity and Structure Alignment to Candidatus Accumulibacter Ppk1")

all_CAP_ppk1_comps_plot

ggsave("figs/all-CAP-ppk1-seq-structure-comps.png", all_CAP_ppk1_comps_plot, width=30, height=20, units=c("cm"))

ppk1_results_metadata %>% 
  filter(alntmscore > 0.95) %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>% 
  ggplot(aes(x=seqid, y=alntmscore)) + 
  geom_point(aes(color=Phylum)) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr(legend = c("bottom"))

top_hits_info <- ppk1_results_metadata %>% 
  filter(alntmscore > 0.95) %>% 
  arrange(desc(alntmscore)) %>% 
  select(accession, seqid, alntmscore, Taxonomic.lineage, Phylum)

stringent_top_hits_info <- ppk1_results_metadata %>% 
  filter(alntmscore > 0.98) %>% 
  arrange(desc(alntmscore))

# retrieve all pseudomonadota results 
ppk1_pseud_accession_list <- ppk1_results_metadata %>% 
  filter(Phylum == 'Pseudomonadota') %>% 
  filter(seqid > 0.6) %>% 
  select(accession)

ppk1_pseud_metadata <- ppk1_results_metadata %>% 
  filter(Phylum == 'Pseudomonadota') %>% 
  filter(seqid > 0.6) %>% 
  select(accession, alntmscore, Organism, Taxonomic.lineage)

# cluster information
structure_clusters <- read.table("results/ppk1_clustering_results/struclusters_features.tsv", header = TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

leiden_clusters <- read.table("results/ppk1_clustering_results/leiden_features.tsv", header = TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

ppk1_pseud_metadata_clusters <- left_join(ppk1_pseud_metadata, structure_clusters) %>% 
  left_join(leiden_clusters)

write.table(ppk1_pseud_accession_list, "metadata/pseudomonadota-ppk1-list.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(ppk1_pseud_metadata_clusters, "metadata/pseudomonadota-ppk1-metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE)

#################################################
# Tetrasphaera japonica ppk1 comparisons
#################################################

# mmseqs easy-search results 
tet_ppk1_mmseqs <- read.table("results/TET_ppk1_sequences_search.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))

# foldseek easy-search results
tet_ppk1_foldseek <- read.table("results/TET_ppk1_structures_search.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

# histograms of seq id and alntm
tet_ppk1_mmseqs %>% 
  ggplot(aes(x=seqid)) +
  geom_histogram()

tet_ppk1_foldseek %>% 
  ggplot(aes(x=alntmscore)) +
  geom_histogram()

# reorganize and combine by accession
tet_ppk1_mmseqs_filtered <- tet_ppk1_mmseqs %>% 
  select(mmseqs_target, seqid) %>% 
  mutate(accession = mmseqs_target) %>% 
  select(accession, seqid)

tet_ppk1_foldseek_filtered <-tet_ppk1_foldseek %>% 
  select(foldseek_target, alntmscore) %>% 
  mutate(accession = gsub("-F1-model_v1.pdb", "", foldseek_target)) %>% 
  select(accession, alntmscore)

tet_ppk1_results <- left_join(tet_ppk1_mmseqs_filtered, tet_ppk1_foldseek_filtered)

# merge with metadata to get taxonomy information
tet_ppk1_results_metadata <- left_join(tet_ppk1_results, ppk1_metadata) %>% 
  mutate(Phylum = if_else(is.na(Phylum), "Other", Phylum))

# plot comparison of protein seqid to alntmscore against Candidatus Accumulibacter ppk1 query 
top_tet_filtered_phyla <- tet_ppk1_results_metadata %>% 
  count(Phylum) %>% 
  top_n(10, n) %>% 
  pull(Phylum)

tet_ppk1_results_metadata %>% 
  mutate(Phylum = if_else(Phylum %in% top_tet_filtered_phyla, Phylum, "Other")) %>% 
  ggplot(aes(x=seqid, y=alntmscore)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr(legend = c(0.9, 0.4)) +
  labs(x="Protein Sequence Identity", y="Protein Strucutre Alignment (tm_score)") +
  ggtitle("Comparisons of Protein Sequence Identity and Structure Alignment to Tetrasphaera japonica Ppk1")
