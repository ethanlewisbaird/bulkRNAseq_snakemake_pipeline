aligner = config["aligner"]
samples = config["samples"]

if aligner == "star":
    qc_inputs = expand("results/alignment/{sample}_Log.final.out", sample=samples)
elif aligner == "hisat2":
    qc_inputs = expand("results/alignment/{sample}_hisat2.log", sample=samples)
else:
    raise ValueError(f"Unsupported aligner: {aligner}")

qc_inputs.append("results/quantification/counts.txt.summary")

rule multiqc:
    input:
        qc_inputs
    output:
        html = "results/multiqc_report/multiqc_report.html",
        data = directory("results/multiqc_report/multiqc_data")
    resources:
        time_min=60,
        mem_mb=8000
    threads: 2
    conda: "../envs/multiqc.yaml"
    shell:
        """
        multiqc results/ -o results/multiqc_report/
        """
