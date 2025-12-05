rule feature_counts:
    input:
        bams = expand(f"results/alignment/{BAM_PATTERN}", sample=config["samples"]),
        gtf = config["gtf_file"]
    output:
        counts = "results/quantification/counts.txt",
        summary = "results/quantification/counts.txt.summary"
    resources: time_min=240, mem_mb=32000
    threads: 16
    conda: "../envs/subread.yaml"  # Ensure subread is installed here
    shell:
        "featureCounts -T {threads} -p -t exon -g gene_id "
        "-a {input.gtf} -o {output.counts} {input.bams}"