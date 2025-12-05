# bulkRNAseq_snakemake_pipeline
Pipeline for bulkRNAseq analysis using several modules including star/hisat alignment, feature counts quantification, edgeR/DEseq2 differential expression, and gene ontology.

To do:
- Increase number of organisms available for Gene Ontology (currently only drosophila is configured)
- Alternative input types (currently only fastq/fastq.gz), e.g. bam
- multiple comparisons deseq2
- gene ontology for the multiple comparisons
- FastQC implementation
- VT EdgeR tricks

####################################################
#### Snakemake pipeline for bulkRNAseq analysis ####
####            Author: Ethan Baird             ####
####################################################

This pipeline is designed to be reproducible and modular. It can be run from any point in the pipeline,
automatically handling dependencies, inputs, outputs...

Contact ethan_baird@imbb.forth.gr with issues or features you would like to see integrated

Minimum required prerequisites:
- Fastq files of reads 
- Miniconda or anaconda installed
- .gtf file of annotations
- Fasta file of reference genome
- Samples metadata .csv

Pipeline modules:
- Star or Hisat2 alignment
- featureCounts quantification
- MultiQC quality control
- DESeq2 or edgeR normalisation and differential expression
- Gene ontology

How to run:

1. Run the system_setup.sh script (Only needs to be run once per system)

chmod +x setup_rnaseq_pipeline.sh # Make it executable

./setup_rnaseq_pipeline.sh # Run script

2. Duplicate the pipeline directory into a project specific directory

cp -r /path/to/pipeline_directory /path/to/new_project_directory/

cd /path/to/project_directory   ### Navigate to the newly created project directory where the Snakefile is located

3. Provide a sample metadata .csv in the config folder

e.g.

sample,condition
SRR18430863,Mi2_RNAi
SRR18430862,Mi2_RNAi
SRR18430861,Mi2_RNAi
SRR18430860,mcherry_RNAi
SRR18430859,mcherry_RNAi
SRR18430858,mcherry_RNAi

4. Configure the config.yaml file with appropriate directories and parameters

5. Start the pipeline with the following command

snakemake --profile slurm
or
snakemake --profile slurm <rule_name> #Replace rule_name with the specific rule to start the pipeline from that rule
or
snakemake --profile slurm --until <rule_name> #Run pipeline until specified rule

see snakemake documentation for other run methods
