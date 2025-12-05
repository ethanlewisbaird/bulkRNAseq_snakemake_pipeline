rule gene_ontology:
    input:
        de_file = "results/differential_expression/differential_expression_results.csv"
    output:
        upregulated_dotplot = "results/gene_ontology/upregulated_dotplot.png",
        upregulated_GO = "results/gene_ontology/upregulated_GO.csv",
        downregulated_dotplot = "results/gene_ontology/downregulated_dotplot.png",
        downregulated_GO = "results/gene_ontology/downregulated_GO.csv"
    resources:
        time_min=120,
        mem_mb=32000
    threads: 8
    conda: "../envs/gene_ontology.yaml"
    script: "../scripts/gene_ontology.R"