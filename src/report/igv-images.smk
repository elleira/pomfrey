rule makePassVCF:
    input:
        snv = "Results/{sample}/Data/{sample}.{support}.filt.vcf",
        pindel = "variantCalls/pindel/{sample}.pindel.ann.vcf",
        artefact = config["bed"]["artefact"],
        germline = config["bed"]["germline"]
    output:
        "Results/{sample}/Reports/{sample}.{support}.PASS.vcf"
    log:
        "logs/report/{sample}.{support}.PASS.vcf.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        """( python3.6 /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/makePASSvcf.py {input.snv} {input.artefact} {input.germline} {output}  && cat {input.pindel} | grep -v '^#' | grep PASS  >> {output} || true ) &>{log}"""

rule createBatFile:
    input:
        vcf = "Results/{sample}/Reports/{sample}.{support}.PASS.vcf",
        indel = "variantCalls/pindel/{sample}.pindel.ann.vcf",
        bam = "Results/{sample}/Data/{sample}.bam",
        bed = config["bed"]["cartool"],
        ref = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/caches/igv/genomes/hg19.genome"
    output:
        "Results/{sample}/Reports/{sample}.{support}-igv.bat"
    params:
        outfolder = "Results/{sample}/Reports/",
        padding = "40",
        sort = "base", #Type of sorting: base, position, strand, quality, sample or readgroup. Could add pos after, but always uses middle.
        view = "squish",  #Type of view, collaps, squished...
        format = "svg" #svg, jpg
    log:
        "logs/report/{sample}.{support}-makeBat.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/report/makeBatfile.py {output} {input.vcf} {input.indel} {input.bam} {input.ref} {input.bed} {params.outfolder} {params.padding} {params.sort} {params.view} {params.format}) &> {log}"

rule igv:
    input:
        bat = "Results/{sample}/Reports/{sample}.{support}-igv.bat"
    output:
        touch("Results/{sample}/Reports/done.{support}-igv.txt")##Several files, add a done.txt
    log:
        "logs/report/{sample}.{support}.igv.log"
    threads:
        2
    singularity:
        config["singularitys"]["igv"]
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} ) &> {log}"
