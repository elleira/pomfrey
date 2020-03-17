localrules: headerBatchMqc, getStatsforMqc

rule samtools_stats:
    input:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam"
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt"
    params:
        extra = "-t "+config["bed"]["bedfile"],                       # Optional: extra arguments.
        # region="1:1000000-2000000"      # Optional: region string.
    log:
        "logs/qc/samtools_stats/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["bwa"]
    wrapper:
        "0.38.0/bio/samtools/stats"

rule picardHsMetrics:
    input:
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        intervals = config["bed"]["intervals"] ##Create with gatk ..
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt"
    log:
        "logs/qc/picardHsMetrics/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(java -Xmx4g -jar /opt/conda/share/picard-2.20.1-0/picard.jar CollectHsMetrics BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals} INPUT={input.bam} OUTPUT={output}) &> {log}"


rule headerBatchMqc:
    input:
        header = config["programdir"]["dir"]+"src/qc/multiqc-header.txt"
    output:
        batch = "Results/batchQC_{seqID}/{seqID}_stats_mqc.tsv"
    log:
        "logs/qc/batchHeader_{seqID}.log"
    # singularity:
    #     config["singularitys"]["python"]
    run:
        import csv
        header = ['Sample','Tot seq','Reads mapped','Avg Coverage','Breadth 500x','Reads paired [%]','Insert size','Insert size s.d.','Average Quality','Duplicates [%]','Target bases 50x','Target bases 100x','Bases on target']
        with open(input.header, 'r') as f:
            with open(output.batch, "w") as file:
                for mqcline in f:
                    file.write(mqcline)
                writer = csv.writer(file, delimiter='\t')
                writer.writerow(header)

rule getStatsforMqc:
    input:
        picardDup = "qc/{sample}_{seqID}/{sample}_{seqID}_DuplicationMetrics.txt",
        picardMet = "qc/{sample}_{seqID}/{sample}_{seqID}.HsMetrics.txt",
        samtools = "qc/{sample}_{seqID}/{sample}_{seqID}.samtools-stats.txt",
        multiQCheader = config["programdir"]["dir"]+"src/qc/multiqc-header.txt",
        cartool = "qc/{sample}_{seqID}/{sample}_{seqID}_Log.csv",
        batch =  "Results/batchQC_{seqID}/{seqID}_stats_mqc.tsv"
    output:
        batchTmp = temp("qc/{sample}_{seqID}/{sample}_batchStats.done"),
        # batch = "qc/{seqID}_stats_mqc.tsv",
        sample = "qc/{sample}_{seqID}/{sample}_{seqID}_stats_mqc.tsv"
    params:
        dir = config["programdir"]["dir"]
    log:
        "logs/qc/{sample}_{seqID}_stats.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/qc/get_stats.py {input.picardDup} {input.picardMet} {input.samtools} {input.multiQCheader} {input.cartool} {output.sample} {input.batch} && touch {output.batchTmp}) &> {log}"
