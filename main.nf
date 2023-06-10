#! /usr/bin/env nextflow

// Description
// small workflow to take input directories of protein fastas and Alphafold PDBs and run mmseqs and foldseek against a query accession, plot the results
// requires downloading directories of protein sequences and structures that have common accession numbers

nextflow.enable.dsl=2

params.threads=4
params.outdir=null

log.info """\
COMPARE PROTEIN SEQUENCE IDENTITY AND STRUCTURE HOMOLOGY TO QUERY PROTEIN
=========================================
query               : $params.query
sequence_dir        : $params.sequence_dir
all_proteins        : $params.all_proteins
structure_dir       : $params.structure_dir
outdir              : $params.outdir
"""

workflow {
    // define channels and check parameters
    sequence_dir = channel.fromPath(params.sequence_dir, checkIfExists: true)
    structure_dir = channel.fromPath(params.structure_dir, checkIfExists: true)
    all_proteins = channel.fromPath(params.all_proteins, checkIfExists: true)

    Channel.from([params.query])
        .map { "${it}.fasta"}
        .set {query_fasta}

    Channel.from([params.query])
        .map { "${it}-F1-model_v1.pdb"}
        .set {query_structure}

    // workflow
    mmseqs_easy_search(sequence_dir, query_fasta, all_proteins)
    foldseek_easy_search(query_structure, structure_dir)

}

process mmseqs_easy_search {
    // run mmseqs easy-search
    tag "${query_fasta}_mmseqs"
    publishDir "${params.outdir}/mmseqs", mode: 'copy', pattern:"*.m8"

    container "quay.io/biocontainers/mmseqs2:14.7e284--pl5321h6a68c12_2"

    input:
    path(sequence_dir)
    val(query_fasta)
    path(all_proteins)

    output:
    path("*.m8"), emit: mmseqs_tsv

    script:
    def args = task.ext.args ?: "--exhaustive-search"
    """
    mmseqs easy-search ${sequence_dir}/${query_fasta} ${all_proteins} mmseqs_result.m8 tmp $args
    """
}

process foldseek_easy_search {
    // run foldseek easy-search

    tag "${query_structure}_foldseek"
    publishDir "${params.outdir}/foldseek", mode: 'copy', pattern:"*.m8"

    container "quay.io/biocontainers/foldseek:6.29e2557--pl5321h6f8c7b7_1"

    input:
    val(query_structure)
    path(structure_dir)

    output:
    path("*.m8"), emit: foldseek_tsv

    script:
    def args = task.ext.args ?: "--format-output query,target,fident,alnlen,alntmscore,qstart,qend,tstart,tend,evalue,bits --exhaustive-search 1"
    """
    foldseek easy-search ${structure_dir}/${query_structure} ${structure_dir} foldseek_result.m8 tmp $args
    """

}
