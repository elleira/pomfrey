rule multiqc:
    input: #Add bcl2fastq folder?
        ["qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip", "data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt", "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt", "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt"] #, "qc/{sample}/{sample}_cartool_mqc.csv" "qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip",
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}_MultiQC.html"
    params:
        extra = "-c "+config["configCache"]["multiqc"] , # Optional: extra parameters for multiqc.
        output_dir = "Results/{sample}_{seqID}/Reports/", #Fix nicer
        output_name = "{sample}_{seqID}_MultiQC.html"
    log:
        "logs/report/multiqc/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["multiqc"]
    shell:
        "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input} ) &> {log}"
    # run: #directly copied from 0.38 wrapper
    #     from os import path
    #     from snakemake.shell import shell
    #
    #     input_dirs = set(path.dirname(fp) for fp in snakemake.input)
    #     output_dir = path.dirname(snakemake.output[0])
    #     output_name = path.basename(snakemake.output[0])
    #     log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    #
    #     shell(
    #         "multiqc"
    #         " {snakemake.params}"
    #         " --force"
    #         " -o {output_dir}"
    #         " -n {output_name}"
    #         " {input_dirs}"
    #         " {log}"
    #     )



rule multiqcBatch:
    input:
        expand("qc/{sample}_{seqID}/{sample}_{seqID}_R1_trimmed_fastqc.zip", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        # expand("qc/{sample}_{seqID}/{sample}_{seqID}_R2_trimmed_fastqc.zip", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        expand("data_processing/{sample}_{seqID}/{sample}_{seqID}.qc.txt", sample=config["samples"], seqID=config["seqID"]["sequencerun"]),
        "Results/batchQC_{seqID}/{seqID}_stats_mqc.json",
        expand("qc/{sample}_{seqID}/{sample}_batchStats.done", sample=config["samples"], seqID=config["seqID"]["sequencerun"]) #Wait until all in table
    output:
        "Results/batchQC_{seqID}/{seqID}_MultiQC.html"
    params:
        extra = "-c "+config["configCache"]["multiqc"]+" --ignore *_{seqID}_stats_mqc.csv", # --ignore *HsMetrics.txt --ignore *samtools-stats.txt",
        output_dir = "Results/batchQC_{seqID}",
        output_name = "{seqID}_MultiQC.html"
    log:
        "logs/report/multiqc/{seqID}.log"
    singularity:
        config["singularitys"]["multiqc"]
    # run: #directly copied from 0.38 wrapper
    #     from os import path
    #     from snakemake.shell import shell
    #
    #     input_dirs = set(path.dirname(fp) for fp in snakemake.input)
    #     output_dir = path.dirname(snakemake.output[0])
    #     output_name = path.basename(snakemake.output[0])
    #     log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    shell:
        "( multiqc {params.extra} --force -o {params.output_dir} -n {params.output_name} {input} ) &> {log}"
