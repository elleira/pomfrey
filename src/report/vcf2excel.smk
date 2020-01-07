def get_minCov(wildcards):
    allCov = config["cartool"]["cov"]
    minCov = allCov.split(' ')[0]
    return minCov

rule makeContainersList:  ##From bedfile, not really dependent on sample
    input:
        expand("Results/{sample}/Reports/{sample}.html", sample=config["samples"])
    output:
        temp("containers.txt")
    log:
        "logs/report/containersLog.log"
    run:
        for k,v in config["singularitys"].items():
            shell("echo {v} >> containers.txt")
        # "(cat slurm-*out | grep singularity | sort | uniq | cut -d' ' -f4 > {output}) &> {log}"


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
        snv =  "Results/{sample}/Data/{sample}.{support}.filt.vcf",
        indel = "variantCalls/pindel/{sample}.pindel.ann.vcf",
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
        get_minCov
    log:
        "logs/report/{sample}.{support}.vcf2excel.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/vcf2excel.py {input.snv} {input.indel} {input.cart} {params} {input.bed} {input.hotspot} {input.artefact} {input.germline} {output}) &> {log}"
