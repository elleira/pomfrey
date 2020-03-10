localrules: fixCoverageHotspot
# def get_minCov(wildcards):
#     allCov = config["cartool"]["cov"]
#     minCov = allCov.split(' ')[0]
#     return minCov

rule fixCoverageHotspot:
    input:
        tsv = "qc/{sample}_{seqID}/{sample}_{seqID}_coverage.tsv",
        bed = config["bed"]["hotspot"]
    output:
        "qc/{sample}_{seqID}/{sample}_{seqID}_coverageShortHotspot.tsv"
    log:
        "logs/report/{sample}_{seqID}.covShortHotspot.log"
    shell:
        """ ( while read line; do chr=$(echo $line | awk '{{print $1}}');pos=$(echo $line | awk '{{print $2}}');cat {input.tsv} | grep ${{chr}} | grep ${{pos}} >>{output} ; done < {input.bed} ) &> {log} """


rule vcf2excel:
    input:
        snv =  "variantCalls/annotation/{sample}_{seqID}.filt.vcf.gz",
        indel = "variantCalls/pindel/{sample}_{seqID}.pindel.filt.vcf.gz",
        cart =  "qc/{sample}_{seqID}/{sample}_{seqID}_MeanCoverageShortList.csv",
        sing = "containers.txt",
        bed = config["bed"]["pindel"],
        hotspot = config["bed"]["hotspot"],
        artefact = config["bed"]["artefact"],
        germline = config["bed"]["germline"],
        shortCov = "qc/{sample}_{seqID}/{sample}_{seqID}_coverageShortHotspot.tsv"
    output:
        "Results/{sample}_{seqID}/Reports/{sample}_{seqID}.xlsx"
    params:
        coverage = config["cartool"]["cov"], #All coverage, goes in as three sys.argv[], get_minCov,
        seqID = config["seqID"]["sequencerun"],
        dir = config["programdir"]["dir"]
    log:
        "logs/report/{sample}_{seqID}.vcf2excel.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 {params.dir}/src/report/vcf2excel.py {input.snv} {input.indel} {params.seqID} {input.cart} {params.coverage} {input.bed} {input.hotspot} {input.artefact} {input.germline} {output}) &> {log}"
