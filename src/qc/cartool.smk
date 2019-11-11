rule cartool:
    input:
        bam = "mapped/{sample}.bam",
        bed = lambda wilcards: config["bed"]["cartool"],
        bai = "mapped/{sample}.bam.bai"
    output:
        statstable = "qc/{sample}/{sample}_Stat_table.csv",
        cartoollog =  "qc/{sample}/{sample}_Log.csv",
        full = "qc/{sample}/{sample}_MeanCoverageFullList.csv",
        short = "qc/{sample}/{sample}_MeanCoverageShortList.csv"
    params:
        user = "arielle",
        coverage = lambda wildcards: config["cartool"]["cov"],
        extra = "-k" #k = combine, p= mapQ
    log:
        "logs/qc/CARTool/{sample}.cartool.log"
    singularity:
        "CARTool.simg"
    shell: #Need to fix -o so no space is needed.
        "( python /opt/CARtool/ProgramLancher.py -a {input.bed} -b {input.bam} -c {params.coverage} -e {params.user} -o qc/{wildcards.sample}/ {wildcards.sample} {params.extra} )&> {log}"


rule fixoutput:
    input:
        header = "../multiqc-header.txt",
        stat = "qc/{sample}/{sample}_Stat_table.csv"
    output:
        "qc/{sample}/{sample}_cartool_mqc.csv"
    log:
        "logs/qc/cartool/{sample}.fix.log"
    shell:
        """ (cat {input.header} > {output} && cut -d ',' -f2- {input.stat} | tr -d "\15\32" >>{output} ) &> {log}"""
