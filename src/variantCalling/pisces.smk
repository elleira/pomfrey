rule pisces:
    input:
      bam = "mapped/{sample}.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
      reffolder = "/data/ref_genomes/hg19/genome_fasta/",
      index = "mapped/{sample}.bam.bai"
    output:
        vcf = "variantCalls/callers/pisces/{sample}/{sample}.vcf"
    params:
        outfolder = "variantCalls/callers/pisces/{sample}",
        bed = lambda wildcards: config["bed"]["bedfile"]
    threads: 1
    log:
        "logs/pisces/{sample}.log"
    singularity:
        "Pisces-5.2.11.simg"
    shell:
        "(dotnet /app/Pisces/Pisces.dll -b {input.bam} -g {input.reffolder} -i {params.bed} -gVCF FALSE -t 1 --filterduplicates TRUE --outfolder {params.outfolder} ) 2> {log}"
##Bed file?

rule piscesFix: ##Might not be needed with -gVCF FALSE
    input:
        "variantCalls/callers/pisces/{sample}/{sample}.vcf"
    output:
        "variantCalls/callers/pisces/{sample}.pisces.unsorted.vcf"
    log:
        "logs/pisces/{sample}.2.log"
    shell:
        """( awk '{{if($5 != "." || $1 ~ /^"#"/)print $0}}' {input} >{output} ) 2> {log}"""

rule sortPisces:
    input:
        "variantCalls/callers/pisces/{sample}.pisces.unsorted.vcf"
    output:
        temp("variantCalls/callers/pisces/{sample}.pisces.weirdAF.vcf")
    singularity:
        "bcftools-1.9--8.simg"
    log:
        "logs/pisces/{sample}.sort.log"
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"
