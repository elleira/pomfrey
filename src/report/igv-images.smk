rule makePassVCF:
    input:
        "variantCalls/annotation/{sample}.3.filt.vcf.gz"
    output:
        temp("reports/{sample}/{sample}.PASS.vcf")
    log:
        "logs/report/{sample}.PASS.vcf.log"
    shell:
        "( zcat {input} | grep -v '^#' | grep PASS >> {output} ) &>{log}"  #( zcat {input} | grep '^#' >{output} &&

rule createBatFile:
    input:
        vcf = "reports/{sample}/{sample}.PASS.vcf",
        bam = "mapped/{sample}.bam",
        bed = lambda wildcards: config["bed"]["cartool"],
        ref = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/caches/igv/genomes/hg19.genome"
    output:
        "reports/{sample}/{sample}-igv.bat"
    params:
        outfolder = "reports/{sample}/",
        padding = "60",
        sort = "base", #Type of sorting: base, position, strand, quality, sample or readgroup. Could add pos after, but always uses middle.
        view = "squish",  #Type of view, collaps, squished...
        format = "svg" #svg, jpg
    log:
        "logs/report/{sample}-makeBat.log"
    shell:
        "(python src/report/makeBatfile.py {output} {input.vcf} {input.bam} {input.ref} {input.bed} {params.outfolder} {params.padding} {params.sort} {params.view} {params.format}) &> {log}"

rule igv:
    input:
        bat = "reports/{sample}/{sample}-igv.bat"
    output:
        "reports/{sample}/done.txt"##Several files, add a done?
    log:
        "logs/report/{sample}.igv.log"
    threads:
        2
    singularity:
        "igv-2.4.10-0.simg"
    shell:
        "(xvfb-run --server-args='-screen 0 3200x2400x24' --auto-servernum --server-num=1 igv.sh -b {input.bat} && touch {output}) &> {log}"
