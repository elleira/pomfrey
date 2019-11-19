#configfile: "/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/samples.yaml"
##Specify configfile and singularity folder in snakemake command.

rule all:
    input:
        expand("Results/{sample}/Reports/{sample}.html", sample=config["samples"]), ##Borde vi addera which run?
        expand("Results/{sample}/Reports/{sample}.xlsx", sample=config["samples"]),
        expand("Results/{sample}/Reports/done-igv.txt", sample=config["samples"]),  ## For the igv images
        expand("Results/{sample}/Data/{sample}.3.filt.vcf", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}.genome.vcf.gz", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}.bam", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}.bam.bai", sample=config["samples"]),
        expand("variantCalls/pindel/{sample}.pindel.ann.vcf", sample=config["samples"]),
        # expand("variantCalls/annotation/{sample}.3.filt.vcf", sample=config["samples"]),
        # expand("variantCalls/recall/{sample}.3.vcf.gz", sample=config["samples"]) ## Reports, final vcf, bam, fastqs..

### QC modules
include:    "qc/fastqc.smk" #fastq in html/text out
include:    "qc/samtools-stats.smk" #bam in txt out
include:    "qc/cartool.smk" #bam in tables out

## Demultiplexing out runfolder/{sample}_S[0-9]_R[12]_001.fastq.g
# include:    "demultiplexing/bcl2fastq.smk"

## Trimming in runfolder/{sample}_S[0-9]_R[12]_001.fastq.gz out trimming/{sample}_R[12]_trimmed.fastq.gz
include:    "trimming/cutadapt.smk"

## Map in trimming/{sample}_R[12]_trimmed.fastq.gz out bam/{sample}.bam
include:    "map/bwa-mem.smk" #fastq R1 R2 from trimming in, bam out.

# include:    "map/indel-realign.smk" #not really better for AF...


## Variant callers
##bamfiles in! then a annotation/{sample}.3.filt.vcf.gz and indel/{sample}.pindel.vcf.gz
include:    "variantCalling/tumor_only.smk"
include:    "variantCalling/pindel.smk"

## CNV?

## Rapportering
include:    "report/multiqc.smk" # per sample, add per batch as well but only certain results?
include:    "report/vcf2excel.smk"
include:    "report/igv-images.smk" #per sample

#Make a results folder:
# Results/Sample/Reports/  and Results/Sample/data/

rule resultFiles:
    input:

    output:
        # "Results/{sample}/Reports/{sample}.html", #Same as reports/{sample}/
        # "Results/{sample}/Reports/{sample}.xlsx",
        # #All svg files? How?
        "Results/{sample}/Data/{sample}.3.filt.vcf",
        "Results/{sample}/Reports/{sample}.genome.vcf.gz",
        "Results/{sample}/Reports/{sample}.bam",
        "Results/{sample}/Reports/{sample}.bam.bai"
    params:

    log:

    shell:
