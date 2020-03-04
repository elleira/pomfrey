
rule bwa_mem:
    input:
        reads=["data_processing/{sample}/{sample}_R1_trimmed.fastq", "data_processing/{sample}/{sample}_R2_trimmed.fastq"]
    output:
        "data_processing/{sample}/{sample}.bam"
    log:
        "logs/map/bwa_mem/{sample}.log"
    params:   #-M
        index="/medstore/External_References/hs37d5/hs37d5.fa",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    singularity:
        config["singularitys"]["bwa"] #bwa 0.7.17, samtools 1.9, picard 2.20.11
    wrapper:
        "0.38.0/bio/bwa/mem"

rule samtools_index:
    input:
        "data_processing/{sample}/{sample}.bam"
    output:
        "data_processing/{sample}/{sample}.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}.log" # optional params string
    singularity:
        config["singularitys"]["bwa"]
    wrapper:
        "0.38.0/bio/samtools/index"
