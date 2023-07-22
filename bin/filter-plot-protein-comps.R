library(tidyverse)

#################################################
# Script to join mmseqs and folseek results and plot comparisons
#################################################

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# inputs 
metadata_tsv             <- args[1]
mmseqs_tsv               <- args[2]
foldseek_tsv             <- args[3]
query                    <- args[4]
outdir                   <- args[5]

# metadata cleaning
metadata_table <- read_tsv(metadata_tsv, col_names = TRUE)

metadata_df <- metadata_table %>% 
  mutate(accession = Entry) %>% 
  select(accession, Taxonomic.lineage, Organism, Phylum) %>% 
  mutate(Phylum = if_else(is.na(Phylum), "Other", Phylum))

# mmseqs tsv
mmseqs_df <- read.table(mmseqs_tsv, sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))

# foldseek tsv
foldseek_df <- read.table(foldseek_tsv, sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

# reorganize tables and combine by accession 
mmseqs_filtered <- mmseqs_df %>% 
  select(mmseqs_target, seqid) %>% 
  mutate(accession = mmseqs_target) %>% 
  select(accession, seqid)

foldseek_filtered <- foldseek_df %>% 
  select(foldseek_target, alntmscore) %>% 
  mutate(accession = gsub(".pdb", "", foldseek_target)) %>% 
  select(accession, alntmscore)

result_df <- left_join(mmseqs_filtered, foldseek_filtered)

# combine with metadata information, get top phyla to assign colors to
result_df_info <- left_join(result_df, metadata_df) %>% 
  mutate(Phylum = if_else(is.na(Phylum), "Other", Phylum))

top_filtered_phyla <- result_df_info %>% 
  count(Phylum) %>% 
  top_n(10, n) %>% 
  pull(Phylum)

plot <- result_df_info %>% 
  mutate(Phylum = if_else(Phylum %in% top_filtered_phyla, Phylum, "Other")) %>% 
  ggplot(aes(x=seqid, y=alntmscore)) +
  geom_point(aes(color=Phylum), alpha=0.5) + 
  scale_color_manual(values = c("#5088C5", "#F28360", "#3B9886", "#F898AE", "#7A77AB", "#F7B846", "#97CD78", "#BAB0A8", "#C85152", "#8A99AD")) + 
  theme_classic() +
  labs(x="Protein Sequence Identity", y="Protein Strucutre Alignment (Tm score)", title = paste("Comparisons of Protein Sequence Identity and Structure Alignment to", query))

ggsave(outdir, plot, width=30, height=20, units=c("cm"))


