library(tidyverse)
library(ggpubr)

#################################################
# Explore and plot structure clusters from foldseek
#################################################

# metadata
ppk1_metadata <- read_tsv("metadata/all-filtered-ppk1-accessions.tsv") %>%
  mutate(accession = Entry) %>%
  select(accession, Taxonomic.lineage, Organism, Phylum) %>%
  mutate(Phylum = replace_na(Phylum, "Other"))

# Reconstruct ppk1_results (used below in left_join) and ppk1_results_metadata (used to
# derive top_filtered_phyla) from the committed mmseqs/foldseek outputs. Previously these
# objects were carried over from scripts/ppk1-seq-vs-structure-comps.R in a shared
# interactive session; rebuilding them here from the same on-disk files (mirroring
# ppk1-seq-vs-structure-comps.R) makes this script self-contained and reproduces the
# identical objects.
ppk1_mmseqs <- read.table("results/A0A369XMZ4_CAP_ppk1_sequences_search.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))
ppk1_foldseek <- read.table("results/A0A369XMZ4_CAP_ppk1_structures_search.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

ppk1_mmseqs_filtered <- ppk1_mmseqs %>%
  select(mmseqs_target, seqid) %>%
  mutate(accession = mmseqs_target) %>%
  select(accession, seqid)

ppk1_foldseek_filtered <- ppk1_foldseek %>%
  select(foldseek_target, alntmscore) %>%
  mutate(accession = gsub("-F1-model_v1.pdb", "", foldseek_target)) %>%
  select(accession, alntmscore)

ppk1_results <- left_join(ppk1_mmseqs_filtered, ppk1_foldseek_filtered)

ppk1_results_metadata <- left_join(ppk1_results, ppk1_metadata) %>%
  mutate(Phylum = if_else(is.na(Phylum), "Other", Phylum))

# clustering results
pca_tsne <- read.table("results/ppk1_clustering_results/all_by_all_tmscore_pivoted_pca_tsne.tsv", header = TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

pca_umap <- read.table("results/ppk1_clustering_results/all_by_all_tmscore_pivoted_pca_umap.tsv", header=TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

# merge clustering results with metadata
pca_tsne_info <- left_join(pca_tsne, ppk1_metadata)

pca_umap_info <- left_join(pca_umap, ppk1_metadata)

# top phylum for colors 
top_filtered_phyla <- ppk1_results_metadata %>% 
  count(Phylum) %>% 
  top_n(10, n) %>% 
  pull(Phylum)

# plotting
pao_wwtp_points <- data.frame(x=c(93.18268), y=c(-16.13617))

pca_tsne_plot <- pca_tsne_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>%
  ggplot(aes(x=tSNE1, y=tSNE2)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr()

pca_tsne_plot +
  geom_point(data = pao_wwtp_points, aes(x=x, y=y), shape = 4, size = 5, color = "black")

pca_umap_plot <- pca_umap_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>%
  ggplot(aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr()

# cluster information
structure_clusters <- read.table("results/ppk1_clustering_results/struclusters_features.tsv", header = TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

structure_clusters_info <- left_join(structure_clusters, pca_tsne_info) %>% 
  left_join(ppk1_results)

structure_clusters_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>%
  ggplot(aes(x=tSNE1, y=tSNE2)) +
  geom_point(aes(color=StruCluster), alpha=0.5) +
  theme_pubr()

# accumulibacter clusters info
acc_clusters_info <- structure_clusters_info %>% 
  filter(StruCluster == 'SC59' | StruCluster == 'SC13' | StruCluster == 'SC21')

write.table(acc_clusters_info, "results/accumulibacter_clusters.tsv", quote = FALSE, row.names = FALSE, sep="\t")

acc_other_clusters <- acc_clusters_info %>% 
  filter(Phylum != "Pseudomonadota") %>% 
  arrange(desc(alntmscore))

# write tables
write.table(structure_clusters_info, "results/PPK1-structure-clusters-info.tsv", quote=FALSE, row.names=FALSE, sep="\t")

# save figures
ggsave("figs/ppk1_pca_tsne_plot_full.jpg", pca_tsne_plot, width=30, height=25, units=c("cm"))

ggsave("figs/ppk1_pca_tsne_plot_full.pdf", pca_tsne_plot, width=11, height=8, units=c("in"))

ggsave("figs/ppk1_pca_umap_plot_full.jpg", pca_umap_plot, width=30, height=25, units=c("cm"))
