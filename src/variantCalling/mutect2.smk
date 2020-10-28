localrules: split_bedfile, bgzipCallers, fixSB
chrom_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

rule Split_bam:
    input:
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bai = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam.bai"
        # vcf = "Results/DNA/{sample}/vcf/{sample}-ensemble.final.no.introns.vcf.gz"
    output:
        bam = temp("variantCalls/callers/mutect2/infiles/{sample,[A-Za-z0-9_-]+}_{seqID}.{chr}.bam"),
        bai = temp("variantCalls/callers/mutect2/infiles/{sample}_{seqID}.{chr}.bam.bai")
    log:
        "logs/split_bam_{sample}_{seqID}-{chr}.log"
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(samtools view -b {input.bam} {wildcards.chr} > {output.bam} && samtools index {output.bam}) &> {log}"

rule split_bedfile:
    input:
        config["bed"]["bedfile"]
    output:
        temp("variantCalls/callers/mutect2/infiles/bedfile_{seqID}.{chr}.bed")
    log:
        "logs/variantCalling/split_bed_{seqID}.{chr}.log"
    shell:
        "(grep -w {wildcards.chr} {input} > {output}) &> {log}"


rule Mutect2:
    input:
        bam = "variantCalls/callers/mutect2/infiles/{sample}_{seqID}.{chr}.bam",
        bai = "variantCalls/callers/mutect2/infiles/{sample}_{seqID}.{chr}.bam.bai",
        fasta = config["reference"]["ref"],
        bed = "variantCalls/callers/mutect2/infiles/bedfile_{seqID}.{chr}.bed"
    output:
        bam = temp("variantCalls/callers/mutect2/perChr/{sample,[A-Za-z0-9_-]+}_{seqID}.{chr}.indel.bam"),
        bai = temp("variantCalls/callers/mutect2/perChr/{sample,[A-Za-z0-9_-]+}_{seqID}.{chr}.indel.bai"),
        stats = temp("variantCalls/callers/mutect2/perChr/{sample}_{seqID}.{chr}.mutect2.unfilt.vcf.gz.stats"),
        vcf = temp("variantCalls/callers/mutect2/perChr/{sample}_{seqID}.{chr}.mutect2.unfilt.vcf.gz")
    log:
        "logs/variantCalling/mutect2_{sample}_{seqID}.{chr}.log"
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' Mutect2 -R {input.fasta} -I {input.bam} -L {input.bed} --bam-output {output.bam} -O {output.vcf}) &> {log}"
# #        "(java -Xmx4g -jar /gatk/gatk-package-4.1.4.1-local.jar Mutect2 -R {input.fasta} -I {input.bam} -L {input.bed} --bam-output {output.bam} -O {output.vcf} ) &> {log}"

rule filterMutect2:
    input:
        vcf = "variantCalls/callers/mutect2/perChr/{sample}_{seqID}.{chr}.mutect2.unfilt.vcf.gz",
        stats = "variantCalls/callers/mutect2/perChr/{sample}_{seqID}.{chr}.mutect2.unfilt.vcf.gz.stats",
        fasta = config["reference"]["ref"]
    output:
        temp("variantCalls/callers/mutect2/perChr/{sample}_{seqID}.{chr}.mutect2.vcf.gz")
    log:
        "logs/variantCalling/mutect2/filter_{sample}_{seqID}.{chr}.log"
    singularity:
        config["singularitys"]["gatk4"]
    shell:
        "(gatk --java-options '-Xmx4g' FilterMutectCalls --max-alt-allele-count 3 --max-events-in-region 8 -R {input.fasta} -V {input.vcf} -O {output}) &> {log}"

rule Merge_vcf:
    input:
        vcf = expand("variantCalls/callers/mutect2/perChr/{{sample}}_{seqID}.{chr}.mutect2.vcf.gz", chr=chrom_list, seqID=config["seqID"]["sequencerun"])
    output:
        "variantCalls/callers/mutect2/{sample,[A-Za-z0-9_-]+}_{seqID}.mutect2.SB.vcf"
    log:
        "logs/variantCalling/mutect2/merge_vcf_{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat -o {output} -O v {input} ) &> {log}"

rule fixSB:
    input:
        "variantCalls/callers/mutect2/{sample}_{seqID}.mutect2.SB.vcf"
    output:
        temp(touch("variantCalls/callers/mutect2/{sample,[A-Za-z0-9_-]+}_{seqID}.SB.done"))
    log:
        "logs/variantCalling/mutect2/{sample}_{seqID}.fixSB.log"
    shell:
        "(sed -i 's/=SB/=SB_mutect2/g' {input}  && sed -i 's/:SB/:SB_mutect2/g' {input}) &> {log}"

rule mutect2HardFilter:
    input:
        vcf = "variantCalls/callers/mutect2/{sample}_{seqID}.mutect2.SB.vcf",
        wait = "variantCalls/callers/mutect2/{sample}_{seqID}.SB.done"
    output:
        temp("variantCalls/callers/mutect2/{sample}_{seqID}.mutect2.weirdAF.vcf")
    params:
        config["programdir"]["dir"]
    log:
        "logs/variantCalling/mutect2/{sample}_{seqID}.hardFilt.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params}/src/variantCalling/hardFilter_mutect2.py {input.vcf} {output}) &> {log}"
        # "(python3.6 {params}/src/variantCalling/hardFilter_mutect2.py {input.vcf} {output}) &> {log}"

rule Merge_bam:
    input:
        bams = expand("variantCalls/callers/mutect2/perChr/{{sample}}_{seqID}.{chr}.indel.bam", chr=chrom_list, seqID=config["seqID"]["sequencerun"])#["mutect2/perChr/{sample}" + str(c) + ".indel.bam" for c in chrom_list]
    output:
        bam = "Results/{sample,[A-Za-z0-9_-]+}_{seqID}/Data/{sample}_{seqID}.indel.bam",
        bai = "Results/{sample}_{seqID}/Data/{sample}_{seqID}.indel.bam.bai"
    log:
        "logs/variantCalling/merge_bam_{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["bwa"]
    shell:
        "(samtools merge {output.bam} {input.bams} && samtools index {output.bam}) &> {log}"
