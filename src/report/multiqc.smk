rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip", "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt", "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt"] #, "qc/{sample}/{sample}_cartool_mqc.csv" "qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip",
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}_MultiQC.html"
    params:
        "-c "+config["configCache"]["multiqc"]  # Optional: extra parameters for multiqc.
    log:
        "logs/report/multiqc/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["multiqc"]
    wrapper:
        "0.38.0/bio/multiqc"


rule multiqcBatch:
    input:
        expand("qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        # expand("qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        expand("data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        "Results/batchQC_{seqID}/{seqID}_stats_mqc.csv",
        expand("qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]) #Wait until all in table
    output:
        "Results/batchQC_{seqID}/{seqID}_MultiQC.html"
    params:
        "-c "+config["configCache"]["multiqc"]+" --ignore *_{seqID}_stats_mqc.csv --ignore *HsMetrics.txt --ignore *samtools-stats.txt"
    log:
        "logs/report/multiqc/{seqID}.log"
    singularity:
        config["singularitys"]["multiqc"]
    wrapper:
        "0.38.0/bio/multiqc"
