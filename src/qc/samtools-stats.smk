rule samtools_stats:
    input:
        "Results/{sample}/Data/{sample}.bam"
    output:
        "qc/{sample}/{sample}.samtools-stats.txt"
    params:
        extra="-t "+config["bed"]["bedfile"],                       # Optional: extra arguments.
        # region="1:1000000-2000000"      # Optional: region string.
    log:
        "logs/samtools_stats/{sample}.log"
    singularity:
        config["singularitys"]["bwa"]
    wrapper:
        "0.38.0/bio/samtools/stats"
