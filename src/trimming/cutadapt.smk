def getFastqs(wildcards):
    fastq1 = config["samples"][wildcards.sample]
    fastq2 = fastq1.replace("_R1_","_R2_")
    return [fastq1, fastq2]

rule cutadapt:
    input:
        getFastqs
    output:
        fastq1="fastq/trimmed/{sample}_R1_trimmed.fastq",
        fastq2="fastq/trimmed/{sample}_R2_trimmed.fastq",
        qc="fastq/trimmed/{sample}/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters_r1 = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  #"-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT",
        adapters_r2 = "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", #"-A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = "--minimum-length 2 -q 20"
    log:
        "logs/cutadapt/{sample}.log"
    threads:    8
    singularity:
        "cutadaptv2.5-0.simg"
    wrapper:
        "0.38.0/bio/cutadapt/pe"
