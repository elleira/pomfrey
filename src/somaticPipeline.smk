#configfile: "/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/samples.yaml"
##Specify configfile and singularity folder in snakemake command.

rule all:
    input:
        expand("Results/{sample}/Reports/{sample}.html", sample=config["samples"]), ##Borde vi addera which run?
        expand("Results/{sample}/Reports/{sample}.3.xlsx", sample=config["samples"]),
        expand("Results/{sample}/Reports/done-igv.txt", sample=config["samples"]),  ## For the igv images
        expand("Results/{sample}/Data/{sample}.SNV-pindel.vcf", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}.normalized.genome.vcf.gz", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}.normalized.genome.vcf.gz.tbi", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}-dedup.bam", sample=config["samples"]),
        expand("Results/{sample}/Data/{sample}-dedup.bam.bai", sample=config["samples"]),
        expand("variantCalls/pindel/{sample}.pindel.filt.vcf.gz", sample=config["samples"]),
        expand("variantCalls/annotation/{sample}.3.filt.vcf.gz", sample=config["samples"])
        # expand("variantCalls/recall/{sample}.3.vcf.gz", sample=config["samples"]) ## Reports, final vcf, bam, fastqs..

wildcard_constraints:
    sample = "[a-zA-Z0-9-_\.]+",
    support = "3" #"\.[0-9]+\."


### QC modules
include:    "qc/fastqc.smk" #fastq in html/text out
include:    "qc/samtools-picard-stats.smk" #bam in txt out
include:    "qc/cartool.smk" #bam in tables out

## Demultiplexing out runfolder/{sample}_S[0-9]_R[12]_001.fastq.g
# include:    "demultiplexing/bcl2fastq.smk"

## Trimming in runfolder/{sample}_S[0-9]_R[12]_001.fastq.gz out trimming/{sample}_R[12]_trimmed.fastq.gz
include:    "trimming/cutadapt.smk"

## Map in trimming/{sample}_R[12]_trimmed.fastq.gz out bam/{sample}.bam
include:    "map/bwa-mem.smk" #fastq R1 R2 from trimming in, bam out.
include:    "map/markDuplicates.smk"

## Variant callers
##bamfiles in! then a annotation/{sample}.3.filt.vcf.gz and indel/{sample}.pindel.vcf.gz
include:    "variantCalling/tumor_only.smk"
include:    "variantCalling/pindel.smk"

## CNV?

## Rapportering
rule makeContainersList:  ##From bedfile, not really dependent on sample
    input:
        expand("data_processing/{sample}/{sample}_R1_trimmed.fastq", sample=config["samples"])
    output:
        temp("containers.txt")
    log:
        "logs/report/containersLog.log"
    run:
        for k,v in config["singularitys"].items():
            shell("echo {v} >> containers.txt")
        "(cat slurm-*out | grep singularity | sort | uniq | cut -d' ' -f4 > {output}) &> {log}"

include:    "report/multiqc.smk" # per sample, add per batch as well but only certain results?
include:    "report/vcf2excel.smk"
include:    "report/igv-images.smk" #per sample
