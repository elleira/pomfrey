rule vep:
    input:
        vcf = "variantCalls/recall/{sample}.{support}.vcf.gz",
        cache = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/caches/vep",#"/opt/vep/.vep", ## always remeber the --bind vep-data:/opt/vep/.vep command in singularity args
        fasta = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta"
    output:
        vcf = temp("variantCalls/annotation/raw/{sample}.{support}.raw.vcf")
    params:
        "--everything --check_existing --pick"  #--exclude_null_alleles
    log:
        "logs/variantCalling/vep/{sample}.{support}.log"
    singularity:
        config["singularitys"]["vep"]
    threads:    8
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --merged --offline --fasta {input.fasta} {params} ) 2> {log}"

rule bgzipVep:
    input:
        "variantCalls/annotation/raw/{sample}.{support}.raw.vcf"
    output:
        "variantCalls/annotation/raw/{sample}.{support}.raw.vcf.gz"
    log:
        "logs/variantCalling/vep/{sample}.{support}.bgzip.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input}) 2> {log}"

rule filterVep:
    input:
        vcf="variantCalls/annotation/raw/{sample}.{support}.raw.vcf.gz"
    output:
        "Results/{sample}/Data/{sample}.{support}.filt.vcf"
    log:
        "logs/variantCalling/vep/filter/{sample}.{support}.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/variantCalling/filter_vcf.py {input.vcf} {output}) 2> {log}"
