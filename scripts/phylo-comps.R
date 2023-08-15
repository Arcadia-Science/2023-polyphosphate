library(ape)
library(adephylo)
library(tidyverse)
library(plotly)
library(webshot)

# read in tree and calculate patristic distance
ppk1_tree <- read.tree("results/ppk1_trees/2023-08-15-itol.tre.txt")
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

# join with metadata information
ppk1_info <- read_tsv("metadata/all-filtered-ppk1-accessions.tsv") %>% 
  select(Entry, Organism, Taxonomic.lineage, Phylum) %>% 
  mutate(target = Entry) %>% 
  mutate(lineage = Taxonomic.lineage) %>% 
  select(-Entry, -Taxonomic.lineage)

ppk1_phylo_struc_info <- left_join(ppk1_phylo_struc_df, ppk1_info)
ppk1_phylo_seq_info <- left_join(ppk1_phylo_seq_df, ppk1_info)

# plot comparisons of phylogenetic distance and seqid to the A0A369XMZ4 Ppk1 query

# highlight specific pathogens
ppk1_phylo_seq_info$highlight <- ifelse(ppk1_phylo_seq_info$target == "Q5FAJ0", "Neisseria gonorrhoeae", ifelse(ppk1_phylo_seq_info$target == "P0DP44", "Pseudomonas aeruginosa", ifelse(ppk1_phylo_seq_info$target == "A0A5B7U1Z3", "Ralstonia solanacearum", ifelse(ppk1_phylo_seq_info$target == "A0A829RFS7", "Acinetobacter baumannii",
                                               "normal"))))
mycolours <- c("Pseudomonas aeruginosa" = "#5088C5", "Neisseria gonorrhoeae" = "#F28360", "Ralstonia solanacearum" = "#7A77AB", "Acinetobacter baumannii" = "#F898AE", "normal" = "grey")

# sequence and phylogenetic distance comparisons
seq_plot <- ppk1_phylo_seq_info %>% 
  ggplot(aes(x=phylogenetic_distance, y=seqid, text=paste("query:", query,
                                                           "\norganism:", Organism))) +
  geom_point(size = 2, aes(colour = highlight), alpha=0.5) +
  scale_color_manual("highlight", values = mycolours) +
  theme_classic()

seq_plot

int_seq_plot <- ggplotly(seq_plot) %>% 
  layout(xaxis=list(title = "Phylogenetic Distance", showgrid=FALSE), yaxis=list(title="Sequence Identity", showgrid=FALSE))

# plot comparison of phylogenetic distance and Tm score to the A0A369XMZ4 Ppk1 query
# highlight specific points
ppk1_phylo_struc_info$highlight <- ifelse(ppk1_phylo_struc_info$target == "Q5FAJ0", "Neisseria gonorrhoeae", ifelse(ppk1_phylo_struc_info$target == "P0DP44", "Pseudomonas aeruginosa", ifelse(ppk1_phylo_struc_info$target == "A0A5B7U1Z3", "Ralstonia solanacearum", ifelse(ppk1_phylo_struc_info$target == "A0A829RFS7", "Acinetobacter baumannii","normal"))))

# structure plot
struc_plot <- ppk1_phylo_struc_info %>% 
  ggplot(aes(x=phylogenetic_distance, y=alntmscore, text=paste("query:", query,
                                                                 "\norganism:", Organism))) +
  geom_point(size = 2, aes(colour = highlight), alpha=0.5) +
  scale_color_manual("highlight", values = mycolours) +
  theme_classic() +
  theme(legend.position = c("none"))

struc_plot

int_struc_plot <- ggplotly(struc_plot) %>% 
  layout(xaxis=list(title = "Phylogenetic Distance", showgrid=FALSE), yaxis=list(title="Structure Similarity (Tmscore)", showgrid=FALSE))

# grid of interactive plots
grid <- subplot(int_seq_plot, int_struc_plot, nrows=2, titleY=TRUE, titleX=TRUE, shareX=TRUE, margin = 0.05)
grid 

# export plots
plotly::export(p = grid,
               file = "figs/phylo-distance-comps.jpeg")

htmlwidgets::saveWidget(
  widget = grid, #the plotly object
  file = "figs/phylo-distance-comps.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)
