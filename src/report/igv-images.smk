rule mergeSNVPindel:
    input:
        snv = "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        pindel = "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz"
    output:
        "Results/{sample}_{seqID}/Data/{sample}_{seqID}.SNV-pindel.vcf"
    params:
        "-a -D -O v" #Allow overlap, Remove duplicates, output format vcf
    log:
        "logs/report/{sample}_{seqID}.mergeVcf.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.snv} {input.pindel} ) &> {log}"

rule makePassVCF:
    input:
        vcf = "Results/{sample}_{seqID}/Data/{sample}_{seqID}.SNV-pindel.vcf",
        artefact = config["bed"]["artefact"],
        germline = config["bed"]["germline"]
    output:
        temp("Results/{sample}_{seqID}/Reports/{sample}_{seqID}.PASS.vcf")
    params:
        config["programdir"]["dir"]
    log:
        "logs/report/{sample}_{seqID}.PASS.vcf.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "( python3.6 {params}/src/report/makePASSvcf.py {input.vcf} {input.artefact} {input.germline} {output} ) &>{log}"

rule createBatFile:
    input:
        vcf = "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.PASS.vcf",
        bam = "Results/{sample}_{seqID}/Data/{sample}_{seqID}-dedup.bam",
        bed = config["bed"]["cartool"],
        ref = config["configCache"]["igv"] #cache
    output:
        "Results/{sample}_{seqID}/Reports/IGV/{sample}_{seqID}-igv.bat"
    params:
        outfolder = "Results/{sample}_{seqID}/Reports/IGV/",
        padding = "40",
        sort = "base", #Type of sorting: base, position, strand, quality, sample or readgroup. Could add pos after, but always uses middle.
        view = "squish",  #Type of view, collaps, squished...
        format = "svg", #svg, jpg
        dir = config["programdir"]["dir"]
    log:
        "logs/report/{sample}_{seqID}-makeBat.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/report/makeBatfile.py {output} {input.vcf} {input.bam} {input.ref} {input.bed} {params.outfolder} {params.padding} {params.sort} {params.view} {params.format}) &> {log}"

rule igv:
    input:
        bat = "Results/{sample}_{seqID}/Reports/IGV/{sample}_{seqID}-igv.bat"
    output:
        touch("Results/{sample}_{seqID}/Reports/IGV/done-igv.txt")##Several files, add a done.txt
    log:
        "logs/report/{sample}_{seqID}.igv.log"
    threads:
        2
    singularity:
        config["singularitys"]["igv"]
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} ) &> {log}"
