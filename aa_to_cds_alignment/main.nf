#!/usr/bin/env nextflow
import org.yaml.snakeyaml.Yaml

// Currently implemented for Nextflow v25.10.2
/**
 * Log a help message.
 */
def helpMessage() {
    help = """
    Usage
    -----

    Run:

    \$ nextflow run main.nf -params-file <YAML_CONFIG_FILE> [--help -c <NEXTFLOW_CONFIG_FILE]

    where:

    * '<YAML_CONFIG_FILE>' is a YAML configuration file. The YAML
      configuration parameters are described below).
    * '--help' displays this help information and exits.
    * 'NEXTFLOW_CONFIG_FILE' is a file that contains many configuration parameters for Nextflow (e.g., queueSize, cpus). This is helpful when running on an HPC system. See https://docs.seqera.io/nextflow/config. 
    * Configuration parameters can also be provided via the
      command-line in the form '--<PARAMETER>=<VALUE>' (for example
      '--min_sequences=10').

   
    Configuration
    -------------

    Input:

    * 'pep_dir': Input directory that contains protein sequences
    * 'cds_dir': Input directory that contains corresponding CDS sequences to pair with files in 'pep_dir'
    
    Output:
    
    * 'pep_alignment_dir': Output directory for aligned protein sequences will go
    * 'untrimmed_cds_alignment_dir': Output directory for aligned CDS sequences before trimming. Note that this will only be used if using pal2nal. Using trimAl will both trim and backtranslate. 
		* 'trimmed_cds_alignment_dir': Output directory for aligned CDS sequences after trimming. Will not be used if no trimming (trim_alignment=false) is specified.
    
    Data filtering:
    
    * min_sequences: Mimimum number of sequences needed in the orthogroup to perform analysis. 
    
    Multi-sequence Alignment:

    * 'alignment_tool': Specify the alignment tool to use. Options are 'mafft', 'muscle', or 'prank'. Default is 'mafft'. NOTE: Have only tested 'mafft' so far.
    * 'alignment_options': Alignment options to pass to aligner. Default '--auto --quiet' is for the default aligner Mafft. This will need to be provided if the user uses any aligner other than Mafft. Can be the empty character "". NOTE: some options are forced, e.g., '-f=fasta' is always specifed for Prank.

    Trim alignment:

    * 'trim_alignment': true or false. Default is true. If true, trim the alignment. If not, backtranslation will occur without trimming using pal2nal
    * 'trim_tool': Trimming software to use (options: 'trimal' or 'clipkit'). The default is 'trimal' which will also provide backtranslation. If ClipKit is used, sequences will first be backtranslated with pal2Nal and then trimmed at the codon level using the '--codon' option. 
    * 'trim_options': Alignment options to pass to trimmer. Default is "-fasta -ignorestopcodon -gappyout" for trimAl. Will need to be modified to if using ClipKit. Can be the empyty character "".
    
    Back translated protein to CDS :
    
    * pal2nal_options: Any options to pass to pal2Nal. Default is '-output fasta -nomismatch'.
    
    Nextflow options:
    
    * mode: Parameter to determine how output files are handled. See https://docs.seqera.io/nextflow/reference/process. Default is 'copy'    """.stripIndent()
    print(help)
}
// Help message implementation, following
// https://github.com/nf-core/rnaseq/blob/master/main.nf (MIT License)
params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

params.min_sequences = 10

params.cds_dir = "test_data/cds/"
params.pep_dir = "test_data/aa/"
params.pep_alignment_dir = "test_data/pep_align/"
params.untrimmed_cds_alignment_dir = "test_data/untrimmed_cds_align/"
params.trimmed_cds_alignment_dir = "test_data/trimmed_cds_align/"


params.alignment_tool = "mafft" //options: mafft, muscle, prank
params.alignment_options = "--auto --quiet" //default for mafft

params.trim_alignment = true
params.trim_tool = "trimal" //options: trimal, clipkit
params.trim_options = "-fasta -ignorestopcodon -gappyout" //default for trimal
params.pal2nal_options = "-output fasta -nomismatch"

params.mode = 'copy' //default for Nextflow is 'symlink'

orthogroup_aa_files = Channel.fromPath(params.pep_dir+"/*.fa")
 .map { file -> tuple(file.baseName, file) }

orthogroup_cds_files = Channel.fromPath(params.cds_dir+"/*.fa")
 .map { file -> tuple(file.baseName, file) }
 
 
//orthogroup_aa_files = orthogroup_aa_files
//	.filter { v -> !(file(params.untrimmed_cds_alignment_dir+"/"+v[0]+".untrimmed_cds_aln").exists()) }

	
//orthogroup_cds_files = orthogroup_cds_files
//	.filter { v -> !(file(params.untrimmed_cds_alignment_dir+"/"+v[0]+".untrimmed_cds_aln").exists()) }


orthogroup_aa_cds = orthogroup_aa_files.combine(orthogroup_cds_files,by: 0)
 
include { alignment;trimAl;trimClipKit;pal2nal } from './modules/alignment'
 
 
workflow pepAlignment {
	take:
	og_aa_file
	
	main:
	pep_align = alignment(og_aa_file)
	
	emit:
	pep_align
} 



workflow {
	
 seq_include = orthogroup_aa_files
		.filter {
			og_id, pep_fasta -> 
        num_seq = 0
        enough_seq = false
        // make sure that the VCF has at least 1 variant, then stop counting
        pep_fasta.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.startsWith(">")) num_seq++
                if (num_seq >= params.min_sequences) {
                    enough_seq = true
                    break
                }
            }
            enough_seq
        }
	}
	pep_align = pepAlignment(seq_include)
	pep_align_cds = orthogroup_cds_files
								 .concat(pep_align)
								 .groupTuple(by: 0, size: 2)
	if (params.trim_alignment & params.trim_tool == "trimal") //trimal performs both trimming and backtranslation
	{
		trimAl(pep_align_cds.map{id,files -> tuple(id,files[1],files[0])})
	} else {
		untrimmed_cds_align = pal2nal(pep_align_cds.map{id,files -> tuple(id,files[1],files[0])})
		if (params.trim_alignment & params.trim_tool == "clipkit")
		{
			trimClipKit(untrimmed_cds_align)
		}
	} 
	
}