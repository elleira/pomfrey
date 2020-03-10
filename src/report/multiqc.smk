rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip", "qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip", "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt", "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt"] #, "qc/{sample}/{sample}_cartool_mqc.csv"
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.html"
    params:
        "-c "+config["configCache"]["multiqc"]  # Optional: extra parameters for multiqc.
    log:
        "logs/report/multiqc/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["multiqc"]
    wrapper:
        "0.38.0/bio/multiqc"
