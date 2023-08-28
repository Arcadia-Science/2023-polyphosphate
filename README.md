# Sequence and Structural Homology of Proteins involved in Microbial Polyphosphate Cycling
This repository contains notebooks, scripts, workflows, and results for exploring the sequence and structural homology of proteins that are involved in microbial polyphosphate cycling.

## Background
Inorganic polyphosphates (polyP) are polymers of orthophosphate and are ubiquitious across life from bacteria to higher eukaryotes. It spans numerous functions such as basic metabolism, sensing/responding to environmental changes, stress survival, sources of ATP, maintenace of cellular structure, etc.

All bacteria have the genetic repertoire for taking in inorganic phosphorus and forming chains of polyphosphate, however some bacteria store substantial amounts of intracellular polyphosphate in response to environmental conditions. This leads to the question: _Why are some microorganisms good at accumulating polyphosphate and others are not?_

**Hypothesis: The structures of proteins involved in polyphosphate metabolism and phosphorus transport are conserved across microorganisms with “enhanced” polyphosphate accumulation, regardless of sequence divergence.**

## Planned Approach
Within polyphosphate-accumulating bacteria from diverse biomes, compare the sequence and structural similarity of proteins involved in inorganic phosphorus transport and polyphosphate metabolism. We will first start with the Ppk1 protein because this protein is crucial for polyphosphate formation. The outline of steps:

1. From Uniprot, collect Ppk1 amino acid sequences and corresponding pre-computed Alphafold protein structural predictions
2. Cluster all Ppk1 proteins with `foldseek`
3. Using the Ppk1 protein from Accumulibacter as a reference, compare against all other Ppk1 protein sequences and structures using `mmseqs` and `foldseek` respectively
4. Plot comparisons of sequence vs structure identity for the query proteins from known PAOs
5. Infer a phylogeny of Ppk1 proteins within the _Pseudomonadota_ phylum, explore phylogenetic distance against protein sequence identity/structural homology

To install the software required for data processing and analysis, you can install with conda:
```
conda env create -n polyphosphate environment.yml
```
### Prep and Download Accessions
Accessions on Uniprot were retrieved by searching "ppk1" and either filtered by taxonomy with "Bacteria" or "Archaea." The resulting accessions were then further filtered by length using the R script `scripts/ppk1-uniprot-accessions-filtering.R` where taxonomy information for each accession is also organized.

Protein FASTA accessions and corresponding Alphafold structures are downloaded with `scripts/download_alphafold_pdbs.py`. For Uniprot the Uniprot accession mostly matches the name of the Alphafold PDB file, but the script still checks against the alphafold accessions CSV that was downloaded from the Alphafold website.

### Cluster all proteins
Using the [`gene-family-cartography`](https://github.com/Arcadia-Science/gene-family-cartography/blob/das/clustering/Cartography_explainer.ipynb) workflow with the `from-folder` configuration, I clustered all ppk1 PDB structure files. The workflow is a Snakemake pipeline that runs with:

```
snakemake --snakefile Snakefile_ff --configfile config_ff_ppk1.yml --cores n
```

And the config file looks like:

```
input_dir: "polyphosphate/protein_structures/structures/"
output_dir: "polyphosphate/results/ppk1_ff"
analysis_name: "ppk1"

features_file: "uniprot_features.tsv"

plotting_modes:
- "pca_tsne"
- "pca_umap"

taxon_focus: 'bac'
```

Where I manually made the `uniprot_features.tsv` file from the prior metadata cleaning I did from when I downloaded lists and metadata of the bacterial and archaeal ppk1 accessions and filtered down to a set I was confident in. This file is in the `polyphosphate/protein_structures/structures/` as the snakemake pipeline expects it to be there with all the PDB files. It is analogous to the `metadata/all-filtered-ppk1-accessions.tsv` file.

### Workflow for `mmseqs` and `foldseek` comparisons to a reference protein accession
The steps for running `mmseqs easy-search` and `foldseek easy-search` and plotting the comparison of protein sequence identity and Tm-score is automated with a Nextflow workflow.

To use the workflow, you will need to have Docker and Nextflow installed:
1. Install Docker [according to these instructions for your operating system](https://docs.docker.com/engine/install/).
2. The easiest way to install Nextflow without worrying about dependency issues on your machine is through a conda environment, and can [install according to the instructions for your operation system](https://docs.conda.io/en/latest/miniconda.html). This is included in the `environment.yml` file.

Then run the workflow with:

```
nextflow run main.nf --query A0A369XMZ4 \\
    --sequence_dir protein_structures/fastas \\
    --structure_dir protein_structures/structures \\
    --all_proteins protein_structures/dbs/all_ppk1_protein_seqs.fasta \\
    --metadata metadata/all-filtered-ppk1-accessions.tsv
    --outdir results
```

Where all protein sequences and folded structures are downloaded and in separate directories, including the query accession. Additionally for this example, all the proteins are appended together into a single FASTA file, which you can point to a directory to do so with `scripts/reorganize_fastas.py`. You can see an example of the figure that is produced by the pipeline in `figs`, and resulting TSV files in `results`. The R script `scripts/ppk1-seq-vs-structure-comps.R` combines the outputs of `mmseqs2` and `foldseek` and plots the comparison of protein sequence identity to Tm score to the provided query protein, which is provided as a command-line R script for this workflow.

### Phylogeny of _Pseudomonadota_ Ppk1 Sequences and Phylogenetic Distance Comparisons
To investigate how phylogenetic distance is related to pairwise sequence identity/structural homology compared to Accumulibacter, I created a phylogenetic tree of Ppk1 proteins within the _Pseudomonadota_ phylum, which is what Accumulibacter as classified within.

1. Cluster sequences at 80% identity using `mmseqs` to get the number of sequences down to a manageable number to make a tree (from ~20,000 to around 1,500)
2. Align sequences with `muscle -super5`
3. Create the phylogenetic tree with `FastTree`
4. Root the tree in iTOL using the outgroup _S. coleicolor_ Ppk1 sequence that was propagated in
5. Visualize in `Empress` coloring by Tmscore, highlighting individuals with > 0.98 Tmscore
6. Use the Rscript `scripts/phylo-comps.R` to calculate pairwise Patrisian distance, and compare to pairwise Seqid and Tmscore to Accumulibacter for the clustered representatives within the _Pseudomonadota_ phylum, which is output as two interactive `plotly` plots for exploration
