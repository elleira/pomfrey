rule freebayes:
    input:
        ref = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta",
        samples = "mapped/{sample}.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
        index = "mapped/{sample}.bam.bai"
    output:
        temp("variantCalls/callers/freebayes/{sample}.freebayes.unsort.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes/{sample}.log"
    singularity:
        "freebayes-1.3.1-0.simg"  ##Not including bcftools and parallel
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
        "bcftools-1.9--8.simg"
    log:
        "logs/freebayes/{sample}.sort.log"
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"
