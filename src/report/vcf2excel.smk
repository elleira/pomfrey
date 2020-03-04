localrules: fixCoverageHotspot
# def get_minCov(wildcards):
#     allCov = config["cartool"]["cov"]
#     minCov = allCov.split(' ')[0]
#     return minCov

rule fixCoverageHotspot:
    input:
        tsv = "qc/{sample}/{sample}_coverage.tsv",
        bed = config["bed"]["hotspot"]
    output:
        "qc/{sample}/{sample}_coverageShortHotspot.tsv"
    log:
        "logs/report/{sample}.covShortHotspot.log"
    shell:
        """ ( while read line; do chr=$(echo $line | awk '{{print $1}}');pos=$(echo $line | awk '{{print $2}}');cat {input.tsv} | grep ${{chr}} | grep ${{pos}} >>{output} ; done < {input.bed} ) &> {log} """


rule vcf2excel:
    input:
        snv =  "variantCalls/annotation/{sample}.{support}.filt.vcf.gz",
        indel = "variantCalls/pindel/{sample}.pindel.filt.vcf.gz",
        cart =  "qc/{sample}/{sample}_MeanCoverageShortList.csv",
        sing = "containers.txt",
        bed = config["bed"]["pindel"],
        hotspot = config["bed"]["hotspot"],
        artefact = config["bed"]["artefact"],
        germline = config["bed"]["germline"],
        shortCov = "qc/{sample}/{sample}_coverageShortHotspot.tsv"
    output:
        "Results/{sample}/Reports/{sample}.{support}.xlsx"
    params:
        coverage = config["cartool"]["cov"], #All coverage, goes in as three sys.argv[], get_minCov,
        seqID = config["seqID"]["sequencerun"]
    log:
        "logs/report/{sample}.{support}.vcf2excel.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/vcf2excel.py {input.snv} {input.indel} {params.seqID} {input.cart} {params.coverage} {input.bed} {input.hotspot} {input.artefact} {input.germline} {output}) &> {log}"
