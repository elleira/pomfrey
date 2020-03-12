localrules: indexSnver, concatSnver
rule snver:
    input:
        bam = "Results/{sample}/Data/{sample}-dedup.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
        index = "Results/{sample}/Data/{sample}-dedup.bam.bai",
        ref = "/medstore/External_References/hs37d5/hs37d5.fa",
        bed = lambda wildcards: config["bed"]["bedfile"]
    output:
        temp("variantCalls/callers/snver/{sample}.snver.raw.vcf"),
        temp("variantCalls/callers/snver/{sample}.snver.filter.vcf"),
        temp("variantCalls/callers/snver/{sample}.snver.indel.raw.vcf"), ##Borde kanske use this instead?
        temp("variantCalls/callers/snver/{sample}.snver.indel.filter.vcf")
    params:
        outfolder = "variantCalls/callers/snver/{sample}.snver"
    log:
        "logs/variantCalling/snver/{sample}.log"
    singularity:
        config["singularitys"]["snver"]
    shell:
        "(snver -i {input.bam} -r {input.ref} -l {input.bed} -o {params.outfolder}) &> {log}"

rule indexSnver:
    input:
        "variantCalls/callers/snver/{sample}.filter.vcf"
    output:
        temp("variantCalls/callers/snver/{sample}.filter.vcf.gz.tbi"),
        "variantCalls/callers/snver/{sample}.filter.vcf.gz"
    log:
        "logs/variantCalling/snver/{sample}.index.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"

rule concatSnver:
    input:
        snver = "variantCalls/callers/snver/{sample}.snver.filter.vcf.gz",
        indel = "variantCalls/callers/snver/{sample}.snver.indel.filter.vcf.gz",
        index = "variantCalls/callers/snver/{sample}.snver.indel.filter.vcf.gz.tbi",
        index2 = "variantCalls/callers/snver/{sample}.snver.filter.vcf.gz.tbi"
    output:
        temp("variantCalls/callers/snver/{sample}.snver.weirdAF.vcf")
    log:
        "logs/variantCalling/snver/concat_{sample}.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat -a -Ou {input.snver} {input.indel} | bcftools sort -Ov -o {output} -) &> {log}"
