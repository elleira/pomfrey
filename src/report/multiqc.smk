rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}/{sample}.samtools-stats.txt", "qc/{sample}/{sample}_fastqc.zip", "fastq/trimmed/{sample}/{sample}.qc.txt","qc/{sample}/{sample}_cartool_mqc.csv"]
    output:
        "Results/{sample}/Reports/{sample}.html"
    params:
        "-c src/report/multiqc_config.yaml"  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/{sample}.log"
    singularity:
        "multiqc-1.7.simg"
    wrapper:
        "0.38.0/bio/multiqc"
