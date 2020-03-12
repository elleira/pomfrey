localrules: fixoutput

rule cartool:
    input:
        bam = "Results/{sample}/Data/{sample}-dedup.bam",
        bed = config["bed"]["cartool"],
        bai = "Results/{sample}/Data/{sample}-dedup.bam.bai"
    output:
        statstable = "qc/{sample}/{sample}_Stat_table.csv",
        cartoollog =  "qc/{sample}/{sample}_Log.csv",
        coverage = "qc/{sample}/{sample}_coverage.tsv",
        full = "qc/{sample}/{sample}_MeanCoverageFullList.csv",
        short = "qc/{sample}/{sample}_MeanCoverageShortList.csv"
    params:
        user = "arielle",
        coverage = config["cartool"]["cov"],
        extra = "-k" #k = combine, p= mapQ
    log:
        "logs/qc/CARTool/{sample}.cartool.log"
    singularity:
        config["singularitys"]["cartool"]
    shell: #Need to fix -o so no space is needed.
        "( python3.6 /opt/CARtool/ProgramLancher.py -a {input.bed} -b {input.bam} -c {params.coverage} -e {params.user} -o qc/{wildcards.sample}/ {wildcards.sample} {params.extra} )&> {log}"


rule fixoutput:
    input:
        header = "/apps/bio/repos/somatic-twist/src/report/multiqc-header.txt", #Always change to correct full path
        stat = "qc/{sample}/{sample}_Stat_table.csv"
    output:
        "qc/{sample}/{sample}_cartool_mqc.csv"
    log:
        "logs/qc/CARTool/{sample}.fix.log"
    shell:
        """ (cat {input.header} > {output} && cut -d ',' -f2- {input.stat} | tr -d "\15\32" >>{output} ) &> {log}"""
