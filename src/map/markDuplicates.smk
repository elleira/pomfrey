rule markDuplicates:
    input:
        bam = "data_processing/{sample}/{sample}.bam",
        bai =  "data_processing/{sample}/{sample}.bam.bai"
    output:
        "Results/{sample}/Data/{sample}-dedup.bam"
    params:
        metric = "qc/{sample}_DuplicationMetrics.txt"
    log:
        "logs/map/{sample}-dedup.log"
    threads:
        4 ##2??
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output} METRICS_FILE={params.metric}) &> {log}"

rule samtools_index_dedup:
    input:
        "Results/{sample}/Data/{sample}-dedup.bam"
    output:
        "Results/{sample}/Data/{sample}-dedup.bam.bai"
    params:
        "" # optional params string
    log:
        "logs/map/samtools_index/{sample}-dedup.log" # optional params string
    singularity:
        config["singularitys"]["bwa"]
    wrapper:
        "0.38.0/bio/samtools/index"
