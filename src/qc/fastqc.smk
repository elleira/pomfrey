rule fastqc:
    input:
        "fastq/trimmed/{sample}/{sample}_R1_trimmed.fastq" ##one for each R1 and one for R2 should be from a samples.yaml file
    output:
        html="qc/{sample}/{sample}_fastqc.html",
        zip="qc/{sample}/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    singularity:
        config["singularitys"]["fastqc"]
    wrapper:
        "0.38.0/bio/fastqc"
