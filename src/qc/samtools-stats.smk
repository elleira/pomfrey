rule samtools_stats:
    input:
        "mapped/{sample}.bam"
    output:
        "qc/{sample}/{sample}.samtools-stats.txt"
    params:
        extra="-t "+config["bed"]["bedfile"],                       # Optional: extra arguments.
        # region="1:1000000-2000000"      # Optional: region string.
    log:
        "logs/samtools_stats/{sample}.log"
    singularity:
        "bwa0.7.17-samtools-1.9.simg"
    wrapper:
        "0.38.0/bio/samtools/stats"
