localrules: piscesFix, sortPisces, gVCFfinalIndex
rule pisces:
    input:
      bam = "Results/{sample}/Data/{sample}-dedup.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
      reffolder = "/medstore/External_References/hs37d5/",
      index = "Results/{sample}/Data/{sample}-dedup.bam.bai"
    output:
        vcf = temp("variantCalls/callers/pisces/{sample}/{sample}-dedup.genome.vcf")
    params:
        outfolder = "variantCalls/callers/pisces/{sample}",
        bed = config["bed"]["bedfile"]
    threads: 1
    log:
        "logs/variantCalling/pisces/{sample}.log"
    singularity:
        config["singularitys"]["pisces"]
    shell:  #Remove gVCF False for genome vcf and save for db, and artifacts? -gVCF FALSE
        "(dotnet /app/Pisces/Pisces.dll -b {input.bam} -g {input.reffolder} -i {params.bed} -t {threads} --filterduplicates TRUE --outfolder {params.outfolder} ) &> {log}"
##Bed file?

rule piscesFix: ## use bcftools view --minalleles 2 {input} instead?
    input:
        "variantCalls/callers/pisces/{sample}/{sample}-dedup.genome.vcf"
    output:
        temp("variantCalls/callers/pisces/{sample}/{sample}.pisces.unsorted.vcf")
    log:
        "logs/variantCalling/pisces/{sample}.2.log"
    shell:
        """( awk '{{if($5 != "." || $1 ~ /^"#"/)print $0}}' {input} >{output} ) 2> {log}"""

rule sortPisces:
    input:
        "variantCalls/callers/pisces/{sample}/{sample}.pisces.unsorted.vcf"
    output:
        temp("variantCalls/callers/pisces/{sample}.pisces.weirdAF.vcf")
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/variantCalling/pisces/{sample}.sort.log"
    shell:
        "(bcftools sort -o {output} -O v {input}) &> {log}"

# rule gVCFindex: #Remove? Not Needed?
#     input:
#         vcf = "variantCalls/callers/pisces/{sample}/{sample}.genome.vcf",
#         wait = "variantCalls/callers/pisces/{sample}.pisces.weirdAF.vcf"
#     output:
#         temp("variantCalls/callers/pisces/{sample}/{sample}.genome.vcf.tbi")
#     log:
#         "logs/variantCalling/pisces/{sample}.index.gVCF.log"
#     singularity:
#         config["singularitys"]["bcftools"]
#     shell:
#         "(tabix {input.vcf} )&>{log}"

rule gVCFdecompose:
    input:
        vcf = "variantCalls/callers/pisces/{sample}/{sample}-dedup.genome.vcf",
        # tbi = "variantCalls/callers/pisces/{sample}/{sample}.genome.vcf.tbi"
    output:
        temp("variantCalls/callers/pisces/{sample}/{sample}.decomp.genome.vcf")
    log:
        "logs/variantCalling/pisces/{sample}.genome.decomp.log"
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output} -) &> {log}"

rule gVCFnormalize:
    input:
        vcf = "variantCalls/callers/pisces/{sample}/{sample}.decomp.genome.vcf",
        fasta = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta"
    output:
        "Results/{sample}/Data/{sample}.normalized.genome.vcf.gz"
    log:
        "logs/variantCalling/pisces/{sample}.normalized.gVCF.log"
    singularity:
        config["singularitys"]["vt"]
    shell:
        "(vt normalize -n -r {input.fasta} -o {output} {input.vcf} ) &> {log}"

rule gVCFfinalIndex:
    input:
        vcf = "Results/{sample}/Data/{sample}.normalized.genome.vcf.gz"
    output:
        "Results/{sample}/Data/{sample}.normalized.genome.vcf.gz.tbi"
    singularity:
        config["singularitys"]["bcftools"]
    log:
        "logs/variantCalling/pisces/{sample}.gz.log"
    shell:
        "(tabix {input.vcf}) 2> {log}"
