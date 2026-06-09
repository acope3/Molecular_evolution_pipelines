# Molecular_evolution_pipelines
Nextflow pipelines to assist is various aspects of molecular evolution, phylogenomics, etc.
Nextflow pipelines are currently being implemented and tested with v25.10.2.
Each pipeline will have a self-contained README.md containing instructions



## aa_to_cds_alignment

Many common software tools rely on codon-level alignments, but these can be difficult to generate using standard nucleotide alignment tools due to, e.g., insertion/deletions in the protein-coding sequences. 
As an alternative approach, it is common to align at the amino-acid level and then translate this back to the nucleotide level, generating codon-level alignments.
These can be used as input into tools such as HyPhy, PAML, IQTree, etc. 

### Requirements

