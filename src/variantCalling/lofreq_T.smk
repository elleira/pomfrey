rule lofreq:
    input:
        "Results/{sample}/Data/{sample}.bam"
    output:
        temp("variantCalls/callers/lofreq/{sample}.lofreq.vcf")
    log:
        "logs/lofreq_call/{sample}.log"
    params:
        ref="/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        extra="-l "+ config["bed"]["bedfile"]
    singularity:
        config["singularitys"]["lofreq"]
    threads: 8
    wrapper:
        "0.38.0/bio/lofreq/call" ##When version 40, add .bai file
