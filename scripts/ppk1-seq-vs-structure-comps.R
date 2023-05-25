library(tidyverse)

# metadata
ppk1_metadata <- all_filtered_ppk1_accessions %>% 
  mutate(accession = Entry) %>% 
  select(accession, Organism, Phylum)

# mmseqs easy-search results 
ppk1_mmseqs <- read.table("results/test_alignment.m8", sep="\t", col.names = c("mmseqs_query", "mmseqs_target", "seqid", "alnlen", "mismatch", "gaps", "qstart", "qend", "tstart", "tend", "evalue", "bits"))

# foldseek easy-search results
ppk1_foldseek <- read.table("results/test_structures.m8", sep="\t", col.names=c("foldseek_query","foldseek_target","fident","alnlen","alntmscore","qstart","qend","tstart","tend","evalue","bits"))

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
ppk1_results_metadata <- left_join(ppk1_results, ppk1_metadata)

# plot comparison of protein seqid to alntmscore against Candidatus Accumulibacter ppk1 query 
ppk1_results_metadata %>% 
  ggplot(aes(x=alntmscore, y=seqid)) +
  geom_point(aes(color=Phylum))
