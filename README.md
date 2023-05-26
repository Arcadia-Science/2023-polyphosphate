# Sequence and Structural Homology of Proteins involved in Microbial Polyphosphate Cycling
This repository contains notebooks, scripts, workflows, and results for exploring the sequence and structural homology of proteins that are involved in microbial polyphosphate cycling.

## Background
Inorganic polyphosphates (polyP) are polymers of orthophosphate and are ubiquitious across life from bacteria to higher eukaryotes. It spans numerous functions such as basic metabolism, sensing/responding to environmental changes, stress survival, sources of ATP, maintenace of cellular structure, etc.

All bacteria have the genetic repertoire for taking in inorganic phosphorus and forming chains of polyphosphate, however some bacteria store substantial amounts of intracellular polyphosphate in response to environmental conditions. This leads to the question: _Why are some microorganisms good at accumulating polyphosphate and others are not?_

**Hypothesis: The structures of proteins involved in polyphosphate metabolism and phosphorus transport are conserved across microorganisms with “enhanced” polyphosphate accumulation, regardless of sequence divergence.**

## Planned Approach
Using proteins from known polyphosphate-accumulating bacteria from diverse biomes, compare the protein sequences and structures of proteins involved in inorganic phosphorus transport and polyphosphate metabolism. We will first start with the Ppk1 protein because this protein is crucial for polyphosphate formation. The outline of steps:

1. Collect Uniprot ppk1 accessions of protein sequences and corresponding Alphafold structures for accessions that already have structures predicted
2. Using Ppk1 proteins from known polyphosphate-accumulating organisms (PAOs) or pathogens etc. known to be good at polyphosphate cycling, compare against all reference protein sequences and structures using `mmseqs` and `foldseek` respectively
3. Plot comparisons of sequence vs structure identity for the query proteins from known PAOs

### Prep and Download Accessions
Accessions on Uniprot were retrieved by searching "ppk1" and either filtered by taxonomy with "Bacteria" or "Archaea." The resulting accessions were then further filtered by length using the R script `scripts/ppk1-uniprot-accessions-filtering.R` where taxonomy information for each accession is also organized.

Protein FASTA accessions and corresponding Alphafold structures are downloaded with `scripts/download_alphafold_pdbs.py`. For Uniprot the Uniprot accession mostly matches the name of the Alphafold PDB file, but the script still checks against the alphafold accessions CSV that was downloaded from the Alphafold website.

### `mmseqs` and `folseek` Comparisons
First run `mmseqs easy-search` using the Accumulibacter ppk1 protein sequence as a query against all downloaded Uniprot ppk1 accessions with:
```
mmseqs easy-search ref_ppk1/A0A369XMZ4.fasta dbs/all_ppk1_protein_seqs.fasta ../results/CAP_ppk1_results.m8 ../results/tmp --exhaustive-search
```

Then run `foldseek easy-search` using the Accumulibacter ppk1 protein Alphafold structure as a query against all downloaded Alphafold ppk1 accessions with:
```
foldseek easy-search ref_ppk1/AF-A0A369XMZ4-F1-model_v4.pdb structures ../results/CAP_ppk1_structures_search.m8 ../results/tmp --format-output "query,target,fident,alnlen,alntmscore,qstart,qend,tstart,tend,evalue,bits" --exhaustive-search 1
```

### Exploration and Visualization
