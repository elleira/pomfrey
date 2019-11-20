rule makePassVCF:
    input:
        snv = "Results/{sample}/Data/{sample}.3.filt.vcf",
        pindel = "variantCalls/pindel/{sample}.pindel.ann.vcf"
    output:
        temp("Results/{sample}/Reports/{sample}.PASS.vcf")
    log:
        "logs/report/{sample}.PASS.vcf.log"
    shell:
        """( zcat {input.snv} | grep -v '^#' | grep PASS >> {output}  && cat {input.pindel} | grep -v '^#' | grep PASS  >> {output} || true ) &>{log}"""

rule createBatFile:
    input:
        vcf = "Results/{sample}/Reports/{sample}.PASS.vcf",
        indel = "variantCalls/pindel/{sample}.pindel.ann.vcf",
        bam = "Results/{sample}/Data/{sample}.bam",
        bed = lambda wildcards: config["bed"]["cartool"],
        ref = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/caches/igv/genomes/hg19.genome"
    output:
        "Results/{sample}/Reports/{sample}-igv.bat"
    params:
        outfolder = "Results/{sample}/Reports/",
        padding = "60",
        sort = "base", #Type of sorting: base, position, strand, quality, sample or readgroup. Could add pos after, but always uses middle.
        view = "squish",  #Type of view, collaps, squished...
        format = "svg" #svg, jpg
    log:
        "logs/report/{sample}-makeBat.log"
    shell:
        "(python src/report/makeBatfile.py {output} {input.vcf} {input.indel} {input.bam} {input.ref} {input.bed} {params.outfolder} {params.padding} {params.sort} {params.view} {params.format}) &> {log}"

rule igv:
    input:
        bat = "Results/{sample}/Reports/{sample}-igv.bat"
    output:
        "Results/{sample}/Reports/done-igv.txt"##Several files, add a done?
    log:
        "logs/report/{sample}.igv.log"
    threads:
        2
    singularity:
        "igv-2.4.10-0.simg"
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} && touch {output}) &> {log}"
