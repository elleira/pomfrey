def get_minCov(wildcards):
    allCov = config["cartool"]["cov"]
    minCov = allCov.split(' ')[0]
    return minCov

rule makeContainersList:  ##Depends on not removing any slurms!
    input:
        expand("reports/{sample}/{sample}.html", sample=config["samples"]),
        expand("variantCalls/pindel/{sample}.pindel.vcf.gz", sample=config["samples"]),
        expand("variantCalls/annotation/{sample}.3.filt.vcf.gz", sample=config["samples"]),
	expand("reports/{sample}/done.txt", sample=config["samples"])
    output:
        temp("containers.txt")
    log:
        "logs/report/containersLog.log"
    shell:
        "(cat slurm-*out | grep singularity | sort | uniq | cut -d' ' -f4 > {output}) &> {log}"

rule vcf2excel:
    input:
        snv =  "variantCalls/annotation/{sample}.3.filt.vcf.gz",
        indel = "variantCalls/pindel/{sample}.pindel.vcf.gz",
        cart =  "qc/{sample}/{sample}_MeanCoverageShortList.csv",
        sing = "containers.txt",
        bed = lambda wildcards: config["bed"]["pindel"]
    output:
        "reports/{sample}/{sample}.xlsx"
    params:
        get_minCov
    log:
        "logs/report/{sample}.vcf2excel.log"
    singularity:
        "python3.6.0-pysam-xlsxwriter.simg"
    shell:
        "(python /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle/somaticpipeline/src/report/vcf2excel.py {input.snv} {input.indel} {input.cart} {params} {input.bed} {output} ) &> {log}"
