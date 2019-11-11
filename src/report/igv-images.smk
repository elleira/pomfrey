rule makePassVCF:
    input:
        "variantCalls/annotation/{sample}.3.filt.vcf.gz"
    output:
        temp("reports/{sample}/{sample}.PASS.vcf")
    log:
        "logs/report/{sample}.PASS.vcf.log"
    shell:
        "( zcat {input} | grep '^#' >{output} && zcat {input} | grep -v '^#' | grep PASS >> {output} ) &>{log}"

rule makeBatFile:
    input:
        vcf = "reports/{sample}/{sample}.PASS.vcf"##Needs to be PASS only vcf!
    output:
        temp("reports/{sample}/igv.temp")
    params:
        outdir = "reports/{sample}/",
        padding = "60",
        sort = "base", #ort: base, position, strand, quality, sample, and readGroup ## Hava to add pos after?? Why not always working?
        view = "-clps", # Need to manually change all collaps to squshi after bedtools,
        format = "svg" #svg, png, eps
    log:
        "logs/report/{sample}.batFile.log"
    singularity:
        "bedtools-2.29.0--3.simg"
    shell:
        "( bedtools igv -path {params.outdir} -slop {params.padding} -sort {params.sort} {params.view} -i {input.vcf} -img {params.format} >> {output}) &> {log}"

def get_bed(wildcards):
    bedFil=config["bed"]["cartool"]
    bed="\/".join(bedFil.split('/'))
    return bed

def get_bam(wildcards):
    bamFil="mapped/"+wildcards.sample+".bam"
    bam="\/".join(bamFil.split('/'))
    return bam

rule fixBatFile:
    input:
        bam = "mapped/{sample}.bam",
        batTmp = "reports/{sample}/igv.temp",
        bedfile =  config["bed"]["cartool"]
    output:
        bat = "reports/{sample}/{sample}-igv.bat"
    params:
        ref = "hg19", ## where does it get the genomefile from, could you specify? Could run with fasta?,
        bamfile = get_bam,
        bedfile = get_bed
    log:
        log = "logs/report/{sample}-fixbat.log"
    shell:
        "( sed '1s/^/new\\ngenome {params.ref}\\nload {params.bedfile}\\nload {params.bamfile}\\n/ ; s/collapse/squish/ ; $aexit' {input.batTmp} >> {output.bat}) &> {log}"
    # run:
    #     import subprocess
    #     bam = "\/".join(input.bam.split('/'))
    #     bedA = "\/".join(input.bedfile.split('/'))
    #     subprocess.call("( sed '1s/^/new\\ngenome "+params.ref+"\\nload "+bedA+"\\nload "+bam+"\\n/ ; s/collapse/squish/ ; $aexit' "+input.batTmp+" >> "+output.bat+") &> "+log.log, shell=True)

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
