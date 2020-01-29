rule vep:
    input:
        vcf = "variantCalls/recall/{sample}.{support}.vcf.gz",
        cache = "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/caches/vep",#"/opt/vep/.vep", ## always remeber the --bind vep-data:/opt/vep/.vep command in singularity args
        fasta = "/data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta"
    output:
        vcf = temp("variantCalls/annotation/raw/{sample}.{support}.raw.vcf")
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class "
        # "--everything --check_existing --pick"  #--exclude_null_alleles
    log:
        "logs/variantCalling/vep/{sample}.{support}.log"
    singularity:
        config["singularitys"]["vep"]
    threads:    8
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --merged --offline --fasta {input.fasta} {params} ) &> {log}"

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
        "(bgzip {input}) &> {log}"

rule filterVep:
    input:
        vcf="variantCalls/annotation/raw/{sample}.{support}.raw.vcf.gz"
    output:
        temp("variantCalls/annotation/{sample}.{support}.filt.vcf")
    log:
        "logs/variantCalling/vep/filter/{sample}.{support}.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/variantCalling/filter_vcf.py {input.vcf} {output}) &> {log}"

rule bgzipSNV:
    input:
        "variantCalls/annotation/{sample}.{support}.filt.vcf"
    output:
        "variantCalls/annotation/{sample}.{support}.filt.vcf.gz",
        "variantCalls/annotation/{sample}.{support}.filt.vcf.gz.tbi"
    log:
        "logs/variantCalling/{sample}.{support}.bgzip.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bgzip {input} && tabix {input}.gz) &> {log}"
