#!/usr/bin/env nextflow

process alignment {
	 label 'process_hpc'
   tag "${og_id}"
   publishDir "${params.pep_alignment_dir}/", mode: "${params.mode}", overwrite: true
   errorStrategy 'ignore'
   input:
       tuple val(og_id), path(pep_file)
   output:
       tuple val(og_id), path("${pep_file.baseName}.pep_aln"), emit: pep_align
   script:
      if ("${params.alignment_tool}" == "mafft")    	
      """
      mafft --version
      mafft ${params.alignment_options} --thread ${task.cpus} ${pep_file} > ${pep_file.baseName}.pep_aln
      """
      else if ("${params.alignment_tool}" == "muscle") //TODO: test
      """
      muscle -version
      muscle ${params.alignment_options} -in {pep_file} -out ${pep_file.baseName}.pep_aln
      """
      else if ("${params.alignment_tool}" == "prank") //TODO: test
      """
      prank -v
      prank ${params.alignment_options} -f=fasta -protein -o=${pep_file.baseName}.pep_aln -d=${pep_file}
      """
      else
				error "Alignment tool ${params.alignment_tool} is not currently implemented."
}


process trimAl {
	label 'process_local'
	tag "${og_id}"
	publishDir "${params.trimmed_cds_alignment_dir}/", mode: "${params.mode}", overwrite: true
  errorStrategy 'ignore'
  input:
  	tuple val(og_id), path(pep_align), path(cds_file)
  output:
  	path("${cds_file.baseName}.cds_aln")
  script:
    """
  	trimal --version
  	trimal -in ${pep_align} -out ${cds_file.baseName}.cds_aln -backtrans ${cds_file} ${params.trim_options}
  	"""
}


process pal2nal {
	label 'process_local'
	tag "${og_id}"
	publishDir "${params.untrimmed_cds_alignment_dir}/", mode: "${params.mode}", overwrite: true
  errorStrategy 'ignore'
  input:
  	tuple val(og_id), path(pep_align), path(cds_file)
  output:
  	tuple val(og_id), path("${cds_file.baseName}.untrimmed_cds_aln"), emit: untrimmed_cds_align
  script:
  	"""
  	pal2nal.pl ${pep_align} ${cds_file} ${params.pal2nal_options} > ${cds_file.baseName}.untrimmed_cds_aln
  	"""
}

process trimClipKit {
	label 'process_local'
	tag "${og_id}"
	publishDir "${params.trimmed_cds_alignment_dir}/", mode: "${params.mode}", overwrite: true
  errorStrategy 'ignore'
  input:
  	tuple val(og_id), path(untrimmed_cds_align)
  output:
  	path("${untrimmed_cds_align.baseName}.cds_aln")
  script:
    """
    clipkit --version
  	clipkit ${untrimmed_cds_align} --codon -s nt -if fasta -of fasta ${params.trim_options} -o ${untrimmed_cds_align.baseName}.cds_aln
  	"""
}
