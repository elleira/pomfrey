rule mergeSNVPindel:
    input:
        snv = "variantCalls/annotation/{sample}.3.filt.vcf.gz",
        pindel = "variantCalls/pindel/{sample}.pindel.filt.vcf.gz"
    output:
        "Results/{sample}/Data/{sample}.SNV-pindel.vcf"
    params:
        "-a -D -O v" #Allow overlap, Remove duplicates, output format vcf
    log:
        "logs/report/{sample}.mergeVcf.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.snv} {input.pindel} ) &> {log}"

rule makePassVCF:
    input:
        vcf = "Results/{sample}/Data/{sample}.SNV-pindel.vcf",
        artefact = config["bed"]["artefact"],
        germline = config["bed"]["germline"]
    output:
        temp("Results/{sample}/Reports/{sample}.PASS.vcf")
    log:
        "logs/report/{sample}.PASS.vcf.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "( python3.6 /apps/bio/repos/somatic-twist/src/report/makePASSvcf.py {input.vcf} {input.artefact} {input.germline} {output} ) &>{log}"

rule createBatFile:
    input:
        vcf = "Results/{sample}/Reports/{sample}.PASS.vcf",
        bam = "Results/{sample}/Data/{sample}-dedup.bam",
        bed = config["bed"]["cartool"],
        ref = "/seqstore/webfolders/igv/genomes/hg19.genome"
    output:
        "Results/{sample}/Reports/{sample}-igv.bat"
    params:
        outfolder = "Results/{sample}/Reports/",
        padding = "40",
        sort = "base", #Type of sorting: base, position, strand, quality, sample or readgroup. Could add pos after, but always uses middle.
        view = "squish",  #Type of view, collaps, squished...
        format = "svg" #svg, jpg
    log:
        "logs/report/{sample}-makeBat.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /apps/bio/repos/somatic-twist/src/report/makeBatfile.py {output} {input.vcf} {input.bam} {input.ref} {input.bed} {params.outfolder} {params.padding} {params.sort} {params.view} {params.format}) &> {log}"

rule igv:
    input:
        bat = "Results/{sample}/Reports/{sample}-igv.bat"
    output:
        touch("Results/{sample}/Reports/done-igv.txt")##Several files, add a done.txt
    log:
        "logs/report/{sample}.igv.log"
    threads:
        2
    singularity:
        config["singularitys"]["igv"]
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} ) &> {log}"
