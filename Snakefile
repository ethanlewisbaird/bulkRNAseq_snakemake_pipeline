configfile: "config/config.yaml"

include: "rules/alignment.smk"
include: "rules/quantification.smk"
include: "rules/differential_expression.smk"
include: "rules/qc.smk"
include: "rules/gene_ontology.smk"

rule all:
    input:
        "results/gene_ontology/upregulated_GO.csv",
        "results/multiqc_report/multiqc_report.html"
