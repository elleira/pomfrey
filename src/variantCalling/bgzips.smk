
rule bgzipCallers:
    input:
        vcf = "variantCalls/callers/{method}/{sample}.{method}.vcf" #[m+"/{sample}."+m+".vcf" for m in config["methods"]]
    output:
        vcf = temp("variantCalls/callers/{method}/{sample}.{method}.vcf.gz"),
        tabix = temp("variantCalls/callers/{method}/{sample}.{method}.vcf.gz.tbi")
    log:
        "logs/variantCalling/bgzip/{method}/{sample}.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input.vcf} && tabix {output.vcf}) 2> {log}"
