library(tidyverse)
library(ArcadiaColorBrewer)
library(ggpubr)

# metadata
ppk1_metadata <- all_filtered_ppk1_accessions %>% 
  mutate(accession = Entry) %>% 
  select(accession, Taxonomic.lineage, Organism, Phylum)

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
pca_tsne_plot <- pca_tsne_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>%
  ggplot(aes(x=tSNE1, y=tSNE2)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_pubr()

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

structure_clusters_info <- left_join(structure_clusters, ppk1_metadata)

structure_clusters_info %>% 
  group_by(StruCluster) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=66)


# save figures
ggsave("figs/ppk1_pca_tsne_plot_full.png", pca_tsne_plot, width=30, height=25, units=c("cm"))

ggsave("figs/ppk1_pca_umap_plot_full.png", pca_umap_plot, width=30, height=25, units=c("cm"))
