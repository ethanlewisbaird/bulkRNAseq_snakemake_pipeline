import os

# Get configuration values
aligner = config.get("aligner", "hisat2").lower()
BAM_PATTERN = "{sample}_Aligned.sortedByCoord.out.bam" if aligner == "star" else "{sample}_hisat2.sorted.bam"
single_end = config.get("single_end", False)

def get_r1(wc):
    # Try paired-end naming first
    base1 = f"{config['samples_dir']}{wc.sample}_1"
    # Try single-end naming (no _1)
    base_no1 = f"{config['samples_dir']}{wc.sample}"
    # Check all possible extensions
    for base in [base1, base_no1]:
        for ext in [".fastq.gz", ".fastq"]:
            path = base + ext
            if os.path.exists(path):
                return path
    raise FileNotFoundError(
        f"Could not find FASTQ file for sample {wc.sample} with or without _1 suffix."
    )

def get_r2(wc):
    if single_end:
        return []  # Return empty list for single-end
    base = f"{config['samples_dir']}{wc.sample}_2"
    for ext in [".fastq.gz", ".fastq"]:
        path = base + ext
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Could not find paired-end FASTQ file for sample {wc.sample} (_2)."
    )

# Function to format read files for STAR
def star_read_files(wildcards, input):
    if single_end:
        return f"{input.r1}"
    else:
        return f"{input.r1} {input.r2}"

# Function to format read files for HISAT2
def hisat2_read_files(wildcards, input):
    if single_end:
        return f"-U {input.r1}"
    else:
        return f"-1 {input.r1} -2 {input.r2}"

# Function to determine read files command
def get_read_files_command(wildcards):
    # Check if R1 file is gzipped
    r1_file = get_r1(wildcards)
    return "zcat" if r1_file.endswith(".gz") else "cat"

# STAR alignment rules (conditional)
if aligner == "star":
    star_index = config.get("star_index", None)
    if star_index:
        # Use prebuilt STAR index, skip building
        rule star_align:
            input:
                index = star_index,
                r1 = get_r1,
                r2 = get_r2 if not single_end else []  # Empty list for single-end
            output:
                bam = f"results/alignment/{BAM_PATTERN}",
                log = "results/alignment/{sample}_Log.final.out"
            params:
                read_files = star_read_files,
                read_files_command = get_read_files_command
            resources: 
                time_min=1440, 
                mem_mb=64000
            threads: 16
            conda: "../envs/star.yaml"
            shell:
                "STAR --runThreadN {threads} "
                "--genomeDir {input.index} "
                "--readFilesIn {params.read_files} "
                "--readFilesCommand {params.read_files_command} "
                "--outFileNamePrefix results/alignment/{wildcards.sample}_ "
                "--outSAMtype BAM SortedByCoordinate"
    else:
        # Build STAR index, then align
        rule star_index:
            input:
                fasta = config["genome_fasta"],
                gtf = config["gtf_file"]
            output:
                directory("results/star_index")
            resources: 
                time_min=1440, 
                mem_mb=64000
            threads: 16
            conda: "../envs/star.yaml"
            shell:
                "mkdir -p {output} && "
                "STAR --runThreadN {threads} "
                "--runMode genomeGenerate "
                "--genomeDir {output} "
                "--genomeFastaFiles {input.fasta} "
                "--sjdbGTFfile {input.gtf} "
                "--sjdbOverhang 199"

        rule star_align:
            input:
                index = "results/star_index",
                r1 = get_r1,
                r2 = get_r2 if not single_end else []  # Empty list for single-end
            output:
                bam = f"results/alignment/{BAM_PATTERN}",
                log = "results/alignment/{sample}_Log.final.out"
            params:
                read_files = star_read_files,
                read_files_command = get_read_files_command
            resources: 
                time_min=1440, 
                mem_mb=64000
            threads: 16
            conda: "../envs/star.yaml"
            shell:
                "STAR --runThreadN {threads} "
                "--genomeDir {input.index} "
                "--readFilesIn {params.read_files} "
                "--readFilesCommand {params.read_files_command} "
                "--outFileNamePrefix results/alignment/{wildcards.sample}_ "
                "--outSAMtype BAM SortedByCoordinate"

# HISAT2 alignment rules (conditional)
elif aligner == "hisat2":
    hisat2_index_prefix = config.get("hisat2_index", None)
    if hisat2_index_prefix:
        # Use prebuilt HISAT2 index, skip building
        hisat2_index_files = [f"{hisat2_index_prefix}.{i}.ht2" for i in range(1,9)]
        rule hisat2_align:
            input:
                hisat2_index = hisat2_index_files,
                r1 = get_r1,
                r2 = get_r2 if not single_end else []  # Empty list for single-end
            output:
                bam = f"results/alignment/{BAM_PATTERN}",
                log = "results/alignment/{sample}_hisat2.log"
            params:
                index_prefix = hisat2_index_prefix,
                read_files = hisat2_read_files
            resources: 
                time_min=1440, 
                mem_mb=32000
            threads: 16
            conda: "../envs/hisat2.yaml"
            shell:
                "hisat2 -p {threads} -x {params.index_prefix} "
                "{params.read_files} "
                "2> {output.log} | "
                "samtools sort -o {output.bam} -"
    else:
        # Build HISAT2 index, then align
        hisat2_index_prefix = "results/hisat2_index/genome"
        hisat2_index_files = expand("results/hisat2_index/genome.{i}.ht2", i=range(1,9))
        
        rule hisat2_index:
            input:
                fa = config["genome_fasta"],
                gtf = config["gtf_file"]
            output:
                hisat2_index = hisat2_index_files
            resources: 
                time_min=1440, 
                mem_mb=64000
            threads: 8
            conda: "../envs/hisat2.yaml"
            shell:
                """
                mkdir -p results/hisat2_index
                hisat2_extract_splice_sites.py {input.gtf} > results/hisat2_index/splice_sites.txt
                hisat2_extract_exons.py {input.gtf} > results/hisat2_index/exons.txt
                hisat2-build --threads {threads} \
                    --ss results/hisat2_index/splice_sites.txt --exon results/hisat2_index/exons.txt \
                    {input.fa} results/hisat2_index/genome
                """

        rule hisat2_align:
            input:
                hisat2_index = hisat2_index_files,
                r1 = get_r1,
                r2 = get_r2 if not single_end else []  # Empty list for single-end
            output:
                bam = f"results/alignment/{BAM_PATTERN}",
                log = "results/alignment/{sample}_hisat2.log"
            params:
                index_prefix = hisat2_index_prefix,
                read_files = hisat2_read_files
            resources: 
                time_min=1440, 
                mem_mb=32000
            threads: 16
            conda: "../envs/hisat2.yaml"
            shell:
                "hisat2 -p {threads} -x {params.index_prefix} "
                "{params.read_files} "
                "2> {output.log} | "
                "samtools sort -o {output.bam} -"
