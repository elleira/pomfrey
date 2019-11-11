rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}/{sample}.samtools-stats.txt", "qc/{sample}/{sample}_fastqc.zip", "fastq/trimmed/{sample}/{sample}.qc.txt","qc/{sample}/{sample}_cartool_mqc.csv"]
    output:
        "reports/{sample}/{sample}.html"
    params:
        "-c /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle/somaticpipeline/src/multiqc_config.yaml"  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/{sample}.log"
    singularity:
        "multiqc-1.7.simg"
    wrapper:
        "0.38.0/bio/multiqc"
