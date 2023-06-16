#! /usr/bin/env nextflow

// Description
// small workflow to take input directories of protein fastas and Alphafold PDBs and run mmseqs and foldseek against a query accession, plot the results
// requires downloading directories of protein sequences and structures that have common accession numbers

nextflow.enable.dsl=2

params.threads=10
params.outdir=null

log.info """\
COMPARE PROTEIN SEQUENCE IDENTITY AND STRUCTURE HOMOLOGY TO QUERY PROTEIN
=========================================
query               : $params.query
metadata            : $params.metadata
sequence_dir        : $params.sequence_dir
all_proteins        : $params.all_proteins
structure_dir       : $params.structure_dir
outdir              : $params.outdir
"""

sequence_dir = channel.fromPath(params.sequence_dir, checkIfExists: true)
structure_dir = channel.fromPath(params.structure_dir, checkIfExists: true)
all_proteins = channel.fromPath(params.all_proteins, checkIfExists: true)
metadata_ch = channel.fromPath(params.metadata, checkIfExists: true)

workflow {
    // define channels and check parameters
    sequence_dir = channel.fromPath(params.sequence_dir, checkIfExists: true)
    structure_dir = channel.fromPath(params.structure_dir, checkIfExists: true)
    all_proteins = channel.fromPath(params.all_proteins, checkIfExists: true)
    metadata_ch = channel.fromPath(params.metadata, checkIfExists: true)
    query = channel.from([params.query])

    // workflow
    mmseqs_easy_search(sequence_dir, query, all_proteins)
    mmseqs_result_ch = mmseqs_easy_search.out.mmseqs_tsv
    foldseek_easy_search(query, structure_dir)
    foldseek_result_ch = foldseek_easy_search.out.foldseek_tsv
    combine_visualize(query, metadata_ch, mmseqs_result_ch, foldseek_result_ch)

}

process mmseqs_easy_search {
    // run mmseqs easy-search
    tag "${query}_mmseqs"
    publishDir "${params.outdir}/mmseqs", mode: 'copy', pattern:"*.m8"

    container "quay.io/biocontainers/mmseqs2:14.7e284--pl5321h6a68c12_2"

    input:
    path(sequence_dir)
    val(query)
    path(all_proteins)

    output:
    path("*.m8"), emit: mmseqs_tsv

    script:
    def args = task.ext.args ?: "--exhaustive-search"
    """
    mmseqs easy-search ${sequence_dir}/${query}.fasta ${all_proteins} ${query}.mmseqs_result.m8 tmp --threads 5 $args
    """
}

process foldseek_easy_search {
    // run foldseek easy-search

    tag "${query}_foldseek"
    publishDir "${params.outdir}/foldseek", mode: 'copy', pattern:"*.m8"

    container "quay.io/biocontainers/foldseek:6.29e2557--pl5321h6f8c7b7_1"

    input:
    val(query)
    path(structure_dir)

    output:
    path("*.m8"), emit: foldseek_tsv

    script:
    def args = task.ext.args ?: "--format-output query,target,fident,alnlen,alntmscore,qstart,qend,tstart,tend,evalue,bits --exhaustive-search 1"
    """
    foldseek easy-search ${structure_dir}/${query}-F1-model_v1.pdb ${structure_dir} ${query}.foldseek_result.m8 tmp --threads 5 $args
    """

}

process combine_visualize {
    // Rscript to combine mmseqs and foldseek results files and visualize

    tag "${query}_viz"
    publishDir "${params.outdir}/figures", mode: 'copy', pattern:"*.png"

    container "rocker/tidyverse"

    input:
    val(query)
    path(metadata)
    path(mmseqs_tsv)
    path(foldseek_tsv)

    output:
    path("*.png"), emit: result_png

    script:
    """
    Rscript bin/filter-plot-protein-comps.R ${metadata} ${mmseqs_tsv} ${foldseek_tsv} ${query}
    """
}
