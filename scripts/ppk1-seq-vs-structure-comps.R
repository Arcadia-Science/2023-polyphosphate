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
ppk1_mmseqs <- read.table("results/A0A369XMZ4_CAP_ppk1_sequences_search.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))

# foldseek easy-search results
ppk1_foldseek <- read.table("results/A0A369XMZ4_CAP_ppk1_structures_search.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

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

ggsave("figs/all-CAP-ppk1-seq-structure-comps.jpg", all_CAP_ppk1_comps_plot, width=30, height=20, units=c("cm"))
ggsave("figs/all-CAP-ppk1-seq-structure-comps.pdf", all_CAP_ppk1_comps_plot, width=11, height=8, units=c("in"))

ppk1_results_metadata %>% 
  filter(alntmscore > 0.95) %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>% 
  ggplot(aes(x=seqid, y=alntmscore)) + 
  geom_point(aes(color=Phylum)) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr(legend = c("bottom"))

# exploring top hits by tmscore
top_hits_info <- ppk1_results_metadata %>% 
  filter(alntmscore > 0.95) %>% 
  arrange(desc(alntmscore)) %>% 
  select(accession, seqid, alntmscore, Taxonomic.lineage, Phylum)

stringent_top_hits_info <- ppk1_results_metadata %>% 
  filter(alntmscore > 0.98) %>% 
  arrange(desc(alntmscore))


#################################### 
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

leiden_clusters <- read.table("results/ppk1_clustering_results/leiden_features.tsv", header=TRUE) %>% 
  mutate(accession = protid) %>% 
  select(-protid)

structure_clusters_info <- left_join(structure_clusters, pca_tsne_info) %>% 
  left_join(ppk1_results)

leiden_clusters_info <- left_join(leiden_clusters, pca_tsne_info) %>% 
  left_join(ppk1_results)

structure_clusters_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>%
  ggplot(aes(x=tSNE1, y=tSNE2)) +
  geom_point(aes(color=StruCluster), alpha=0.5) +
  theme_pubr()

sc59 <- structure_clusters_info %>% 
  filter(StruCluster == 'SC59') %>% 
  filter(Phylum == 'Pseudomonadota') %>% 
  mutate(info = "SC59_cluster")

high_hits <- ppk1_results_metadata %>% 
  filter(alntmscore > 0.98) %>% 
  mutate(hit = 'high98')
  

# save figures
ggsave("figs/ppk1_pca_tsne_plot_full.png", pca_tsne_plot, width=30, height=25, units=c("cm"))

ggsave("figs/ppk1_pca_umap_plot_full.png", pca_umap_plot, width=30, height=25, units=c("cm"))


# pseudomonadota-specific information 
# retrieve all pseudomonadota hits
ppk1_all_pseud <- ppk1_results_metadata %>% 
  filter(Phylum == 'Pseudomonadota') %>% 
  select(accession)

ppk1_all_pseud_metadata <- ppk1_results_metadata %>% 
  filter(Phylum == 'Pseudomonadota')

# cluster information
ppk1_pseud_metadata_clusters <- left_join(ppk1_all_pseud_metadata, high_hits) %>% 
  left_join(structure_clusters_info)

ppk1_pseud_leiden_metadata <- left_join(ppk1_all_pseud_metadata, lc32)

write.table(ppk1_pseud_accession_list, "metadata/pseudomonadota-ppk1-list.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(ppk1_pseud_metadata_clusters, "metadata/pseudomonadota-ppk1-metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE)

write.table(ppk1_pseud_leiden_metadata, "metadata/pseudomonadota-leiden-ppk1-metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE)

write.table(ppk1_all_pseud, "metadata/all-pseudo-ppk1-list.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
