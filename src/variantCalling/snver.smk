rule snver:
    input:
        bam = "Results/{sample}/Data/{sample}.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
        index = "Results/{sample}/Data/{sample}.bam.bai",
        ref = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        bed = lambda wildcards: config["bed"]["bedfile"]
    output:
        temp("variantCalls/callers/snver/{sample}.snver.raw.vcf"),
        temp("variantCalls/callers/snver/{sample}.snver.filter.vcf"),
        temp("variantCalls/callers/snver/{sample}.snver.indel.raw.vcf"), ##Borde kanske use this instead?
        temp("variantCalls/callers/snver/{sample}.snver.indel.filter.vcf")
    params:
        outfolder = "variantCalls/callers/snver/{sample}.snver"
    log:
        "logs/snver/{sample}.log"
    singularity:
        "snver-0.5.3-0.simg"
    shell:
        "(snver -i {input.bam} -r {input.ref} -l {input.bed} -o {params.outfolder}) &> {log}"

rule indexSnver:
    input:
        "variantCalls/callers/snver/{sample}.filter.vcf"
    output:
        temp("variantCalls/callers/snver/{sample}.filter.vcf.gz.tbi"),
        "variantCalls/callers/snver/{sample}.filter.vcf.gz"
    log:
        "logs/snver/{sample}.index.log"
    singularity:
        "bcftools-1.9--8.simg"
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
        "logs/snver/concat_{sample}.log"
    singularity:
        "bcftools-1.9--8.simg"
    shell:
        "(bcftools concat -a -Ou {input.snver} {input.indel} | bcftools sort -Ov -o {output} -) &> {log}"
