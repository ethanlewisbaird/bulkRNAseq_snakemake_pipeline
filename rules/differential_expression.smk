diffexp_tool = config.get("diffexp_tool", "DESeq2").lower()

if diffexp_tool == "deseq2":
    rule differential_expression:
        input:
            counts = "results/quantification/counts.txt",
            sample_info = "config/samples.csv"
        output:
            normalized_counts = "results/differential_expression/normalized_counts.csv",
            de_results = "results/differential_expression/differential_expression_results.csv"
        resources:
            time_min=120,
            mem_mb=16000
        threads: 4
        conda: "../envs/deseq2.yaml"
        script: "../scripts/deseq2.R"

elif diffexp_tool == "edger":
    rule differential_expression:
        input:
            counts = "results/quantification/counts.txt",
            sample_info = "config/samples.csv"
        output:
            normalized_counts = "results/differential_expression/normalized_counts.csv",
            de_results = "results/differential_expression/differential_expression_results.csv"
        resources:
            time_min=120,
            mem_mb=16000
        threads: 4
        conda: "../envs/edger.yaml"
        script: "../scripts/edger.R"