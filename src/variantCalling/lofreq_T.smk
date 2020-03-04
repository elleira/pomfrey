rule lofreq:
    input:
        "Results/{sample}/Data/{sample}-dedup.bam"
    output:
        temp("variantCalls/callers/lofreq/{sample}.lofreq.vcf")
    log:
        "logs/variantCalling/lofreq_call/{sample}.log"
    params:
        ref = config["reference"]["ref"],
        extra = "-l "+ config["bed"]["bedfile"]
    singularity:
        config["singularitys"]["lofreq"]
    threads: 8
    wrapper:
        "0.38.0/bio/lofreq/call" ##When version 40, add .bai file
