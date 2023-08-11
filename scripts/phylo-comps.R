library(ape)
library(adephylo)
library(tidyverse)

# read in tree and calculate patristic distance
ppk1_tree <- read.tree("results/ppk1_trees/pseud_reps_rerooted.tre")
ppk1_distances <- distTips(ppk1_tree, method = "patristic")
ppk1_distances_df <- ppk1_distances %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("query") %>% 
  gather(key = "target", value = "phylogenetic_distance", -query)

A0A369XMZ4_comparisons <- ppk1_distances_df %>% 
  filter(query == 'A0A369XMZ4') %>% 
  filter(target != 'Q9KZV6') # remove the outgroup comparison

# sequence similarity results for A0A369XMZ4 PPK1
acc_ppk1_mmseqs_table <- read.table("results/CAP_ppk1_sequences_search.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits")) %>% 
  select(mmseqs_query, mmseqs_target, seqid) %>% 
  mutate(query = mmseqs_query) %>% 
  mutate(target = mmseqs_target) %>% 
  select(-mmseqs_query, -mmseqs_target)

# structural similarity results for A0A369XMZ4 PPK1
acc_ppk1_foldseek_table <- read.table("results/CAP_ppk1_structures_search.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits")) %>% 
  select(foldseek_query, foldseek_target, alntmscore) %>% 
  mutate(query = gsub("-F1-model_v4.pdb", "", foldseek_query)) %>% 
  mutate(query = gsub("AF-", "", query)) %>% 
  mutate(target = gsub("-F1-model_v1.pdb", "", foldseek_target)) %>% 
  mutate(target = gsub("AF-", "", target)) %>% 
  select(query, target, alntmscore)

# structural clusters information
structure_clusters <- read.table("results/ppk1_clustering_results/struclusters_features.tsv", sep="\t", header = TRUE) %>% 
  mutate(target = protid) %>% 
  mutate(cluster = StruCluster) %>% 
  select(target, cluster)

# joined table of phylogenetic distance and seqid for all comparisons to A0A369XMZ4 PPK1 with structure clusters info
ppk1_phylo_seq_df <- left_join(acc_ppk1_mmseqs_table, A0A369XMZ4_comparisons) %>% 
  left_join(structure_clusters) %>% 
  drop_na() 

# joined table of phylogenetic distance and alntmscores for all comparisons to A0A369XMZ4 PPK1 with structure clusters info
ppk1_phylo_struc_df <- left_join(acc_ppk1_foldseek_table, A0A369XMZ4_comparisons) %>% 
  left_join(structure_clusters) %>% 
  drop_na() # only want comparisons in the tree

top_clusters <- ppk1_phylo_struc_df %>% 
  count(cluster) %>% 
  top_n(20, n) %>% 
  pull(cluster)

# plot comparison of phylogenetic distance and seqid to the A0A369XMZ4 Ppk1 query
ppk1_phylo_seq_df %>% 
  mutate(cluster = if_else(cluster %in% top_clusters, cluster, "other")) %>% 
  ggplot(aes(x=phylogenetic_distance, y=seqid)) +
  geom_point(aes(color=cluster)) +
  theme_classic()

# plot comparison of phylogenetic distance and Tm score to the A0A369XMZ4 Ppk1 query
ppk1_phylo_struc_df %>% 
  mutate(cluster = if_else(cluster %in% top_clusters, cluster, "other")) %>% 
  ggplot(aes(x=phylogenetic_distance, y=alntmscore)) +
  geom_point(aes(color=cluster)) +
  theme_classic()
