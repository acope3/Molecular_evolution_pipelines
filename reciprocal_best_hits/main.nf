#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.query   = null
params.subject = null
params.outdir  = 'results'
params.evalue  = 1e-5
params.max_target_seqs = 5
params.outfmt_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
params.threads = 4

if (!params.query || !params.subject) {
    error "Please provide --query and --subject protein FASTA files"
}

process MAKEBLASTDB {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/tmp/blastdbs", mode: 'copy'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta.baseName}_db*"), emit: db

    script:
    """
    makeblastdb -in ${fasta} -dbtype prot -out ${fasta.baseName}_db
    """
}

process BLASTP_FORWARD {
    tag "query_vs_subject"
    publishDir "${params.outdir}/tmp/blast_results", mode: 'copy'

    input:
    tuple path(subject_fasta), path(db_files)
    path query_fasta

    output:
    path 'forward_blast.tsv'

    script:
    """
    blastp \
        -query ${query_fasta} \
        -db ${subject_fasta.baseName}_db \
        -evalue ${params.evalue} \
        -max_target_seqs ${params.max_target_seqs} \
        -outfmt "6 ${params.outfmt_cols}" \
        -num_threads ${params.threads} \
        -out forward_blast.tsv
    """
}

process BLASTP_REVERSE {
    tag "subject_vs_query"
    publishDir "${params.outdir}/tmp/blast_results", mode: 'copy'

    input:
    tuple path(query_fasta), path(db_files)
    path subject_fasta

    output:
    path 'reverse_blast.tsv'

    script:
    """
    blastp \
        -query ${subject_fasta} \
        -db ${query_fasta.baseName}_db \
        -evalue ${params.evalue} \
        -max_target_seqs ${params.max_target_seqs} \
        -outfmt "6 ${params.outfmt_cols}" \
        -num_threads ${params.threads} \
        -out reverse_blast.tsv
    """
}

process RECIPROCAL_BEST_HITS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path forward_hits
    path reverse_hits

    output:
    path 'reciprocal_best_hits.tsv'

    script:
    """
    find_reciprocal_best_hits.py \
        --forward ${forward_hits} \
        --reverse ${reverse_hits} \
        --outfmt_cols "${params.outfmt_cols}" \
        --output reciprocal_best_hits.tsv
    """
}

workflow {
    query_ch   = Channel.fromPath(params.query, checkIfExists: true)
    subject_ch = Channel.fromPath(params.subject, checkIfExists: true)

    // Build BLAST databases for both
    MAKEBLASTDB(query_ch.mix(subject_ch))

    // Separate outputs by matching to original filenames
    query_db = MAKEBLASTDB.out.db.filter { tuple ->
        tuple[0].name == file(params.query).name
    }
    subject_db = MAKEBLASTDB.out.db.filter { tuple ->
        tuple[0].name == file(params.subject).name
    }

    // Forward search: query sequences against subject database
    BLASTP_FORWARD(subject_db, query_ch)

    // Reverse search: subject sequences against query database
    BLASTP_REVERSE(query_db, subject_ch)

    // Identify reciprocal best hits
    RECIPROCAL_BEST_HITS(BLASTP_FORWARD.out, BLASTP_REVERSE.out)
}
