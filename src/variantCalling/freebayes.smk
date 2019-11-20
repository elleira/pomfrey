rule freebayes:
    input:
        ref = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        samples = "Results/{sample}/Data/{sample}.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
        index = "Results/{sample}/Data/{sample}.bam.bai"
    output:
        temp("variantCalls/callers/freebayes/{sample}.freebayes.unsort.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes/{sample}.log"
    singularity:
        config["singularitys"]["freebayes"]  ##Not including bcftools and parallel
    params:
        extra = " --min-alternate-fraction 0.01 --allele-balance-priors-off --pooled-discrete --pooled-continuous --report-genotype-likelihood-max -t " +config["bed"]["bedfile"],         # optional parameters. Add regions file, bed-format.
        chunksize = 100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 1
    wrapper:
        "0.34.0/bio/freebayes"

rule sortFreebayes:
    input:
        "variantCalls/callers/freebayes/{sample}.freebayes.unsort.vcf"
    output:
        temp("variantCalls/callers/freebayes/{sample}.freebayes.weirdAF.vcf")
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/freebayes/{sample}.sort.log"
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"
