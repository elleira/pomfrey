rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}/{sample}.samtools-stats.txt", "qc/{sample}/{sample}_fastqc.zip", "fastq/trimmed/{sample}/{sample}.qc.txt","qc/{sample}/{sample}_cartool_mqc.csv"]
    output:
        "Results/{sample}/Reports/{sample}.html"
    params:
        "-c /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/multiqc_config.yaml"  # Optional: extra parameters for multiqc.
    log:
        "logs/report/multiqc/{sample}.log"
    singularity:
        config["singularitys"]["multiqc"]
    wrapper:
        "0.38.0/bio/multiqc"
