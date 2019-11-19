
rule bwa_mem:
    input:
        reads=["fastq/trimmed/{sample}_R1_trimmed.fastq", "fastq/trimmed/{sample}_R2_trimmed.fastq"]
    output:
        "Results/{sample}/Data/{sample}.bam"
    log:
        "logs/map/bwa_mem/{sample}.log"
    params:   #-M
        index="/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    singularity:
        "bwa0.7.17-samtools-1.9.simg" #bwa 0.7.17, samtools 1.9, picard 2.20.11
    wrapper:
        "0.38.0/bio/bwa/mem"

rule samtools_index:
    input:
        "Results/{sample}/Data/{sample}.bam"
    output:
        "Results/{sample}/Data/{sample}.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}.log" # optional params string
    singularity:
        "bwa0.7.17-samtools-1.9.simg"
    wrapper:
        "0.38.0/bio/samtools/index"
