rule samtools_stats:
    input:
        "Results/{sample}/Data/{sample}-dedup.bam"
    output:
        "qc/{sample}/{sample}.samtools-stats.txt"
    params:
        extra = "-t "+config["bed"]["bedfile"],                       # Optional: extra arguments.
        # region="1:1000000-2000000"      # Optional: region string.
    log:
        "logs/qc/samtools_stats/{sample}.log"
    singularity:
        config["singularitys"]["bwa"]
    wrapper:
        "0.38.0/bio/samtools/stats"

rule picardHsMetrics:
    input:
        bam = "Results/{sample}/Data/{sample}-dedup.bam",
        intervals = config["bed"]["intervals"] ##Create with gatk ..
    output:
        "qc/{sample}/{sample}.HsMetrics.txt"
    log:
        "logs/qc/picardHsMetrics/{sample}.log"
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectHsMetrics BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals} INPUT={input.bam} OUTPUT={output}) &> {log}"
