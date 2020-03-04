localrules: fixContigPindel, pindelConf, fixPindelDPoAF, filterPindel, bgzipPindel

rule pindelConf: ##Add in excel file what genes were used.
    input:
        bam = "Results/{sample}/Data/{sample}-dedup.bam",
        bai = "Results/{sample}/Data/{sample}-dedup.bam.bai"
    output:
        "variantCalls/pindel/{sample}/{sample}-config.txt"
    log:
        "logs/variantCalling/pindel/{sample}.config.log"
    shell:
        "( echo -e '{input.bam}\t250\t{wildcards.sample}'>{output} ) &> {log}"

rule pindel:
    input:
        bed = lambda wildcards: config["bed"]["pindel"],
        ref = "/medstore/External_References/hs37d5/hs37d5.fa",
        bamconfig = "variantCalls/pindel/{sample}/{sample}-config.txt" #path to bam \t insert size \t sample name
    output:
        "variantCalls/pindel/{sample}/{sample}_BP",
        "variantCalls/pindel/{sample}/{sample}_CloseEndMapped",
        "variantCalls/pindel/{sample}/{sample}_D",
        "variantCalls/pindel/{sample}/{sample}_INT_final",
        "variantCalls/pindel/{sample}/{sample}_INV",
        "variantCalls/pindel/{sample}/{sample}_LI",
        "variantCalls/pindel/{sample}/{sample}_RP",
        "variantCalls/pindel/{sample}/{sample}_SI",
        "variantCalls/pindel/{sample}/{sample}_TD"
    params:
        x = 2,
        B = 60
    log:
        "logs/variantCalling/pindel/{sample}.pindel.log"
    singularity:
        config["singularitys"]["pindel"]
    threads:    4
    shell:
        " (pindel -f {input.ref} -i {input.bamconfig} -T {threads} -x {params.x} -B {params.B} -j {input.bed} -o variantCalls/pindel/{wildcards.sample}/{wildcards.sample} ) &> {log}"

rule pindel2vcf:
    input:
        ref = "/medstore/External_References/hs37d5/hs37d5.fa",
        bp = "variantCalls/pindel/{sample}/{sample}_BP",
        closeend = "variantCalls/pindel/{sample}/{sample}_CloseEndMapped",
        d = "variantCalls/pindel/{sample}/{sample}_D",
        final = "variantCalls/pindel/{sample}/{sample}_INT_final",
        inv = "variantCalls/pindel/{sample}/{sample}_INV",
        li = "variantCalls/pindel/{sample}/{sample}_LI",
        rp = "variantCalls/pindel/{sample}/{sample}_RP",
        si = "variantCalls/pindel/{sample}/{sample}_SI",
        td = "variantCalls/pindel/{sample}/{sample}_TD"
    output:
        temp("variantCalls/pindel/{sample}.pindel.noDP.noContig.vcf")
    params:
        e = 10, #min supporting reads 35
        mc = 10, #min coverage
        minsize = 10, #min size of reported 5
        refname = "hg19",
        refdate = 000000  #Can I add seqID instead? config["seqID"]["sequencerun"]
    log:
        "logs/variantCalling/pindel/{sample}.pindel2vcf.log"
    singularity:
        config["singularitys"]["pindel"]
    threads:    1
    shell:
        "(pindel2vcf -P variantCalls/pindel/{wildcards.sample}/{wildcards.sample} -r {input.ref} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize} ) &> {log}"

rule fixContigPindel:
    input:
        "variantCalls/pindel/{sample}.pindel.noDP.noContig.vcf"
    output:
        temp("variantCalls/pindel/{sample}.pindel.noDP.vcf")
    params: ## awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
        "\"##contig=<ID=1,length=249250621>\\n##contig=<ID=2,length=243199373>\\n##contig=<ID=3,length=198022430>\\n##contig=<ID=4,length=191154276>\\n##contig=<ID=5,length=180915260>\\n##contig=<ID=6,length=171115067>\\n##contig=<ID=7,length=159138663>\\n##contig=<ID=8,length=146364022>\\n##contig=<ID=9,length=141213431>\\n##contig=<ID=10,length=135534747>\\n##contig=<ID=11,length=135006516>\\n##contig=<ID=12,length=133851895>\\n##contig=<ID=13,length=115169878>\\n##contig=<ID=14,length=107349540>\\n##contig=<ID=15,length=102531392>\\n##contig=<ID=16,length=90354753>\\n##contig=<ID=17,length=81195210>\\n##contig=<ID=18,length=78077248>\\n##contig=<ID=19,length=59128983>\\n##contig=<ID=20,length=63025520>\\n##contig=<ID=21,length=48129895>\\n##contig=<ID=22,length=51304566>\\n##contig=<ID=X,length=155270560>\\n##contig=<ID=Y,length=59373566>\\n##contig=<ID=MT,length=16569>\\n\""
    log:
        "logs/variantCalling/pindel/{sample}.fixContig.log"
    shell:
        """(cat {input} | grep -v "^##contig" | awk '/^#CHROM/ {{ printf({params});}} {{print;}}' > {output} )&> {log}  """

rule fixPindelDPoAF:
    input:
        "variantCalls/pindel/{sample}.pindel.noDP.vcf"
    output:
        "variantCalls/pindel/{sample}.pindel.vcf"
    log:
        "logs/variantCalling/{sample}.fixDP.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /apps/bio/repos/somatic-twist/src/variantCalling/fix_pindelDPoAF.py {input} {output}) &> {log}"

rule annotatePindel:
    input:
        vcf = "variantCalls/pindel/{sample}.pindel.vcf",
        fasta = "/medstore/External_References/hs37d5/hs37d5.fa",
        cache = "/medstore/External_References/VEP/vep-data-99.0"
    output:
        temp("variantCalls/pindel/{sample}.pindel.ann.vcf")
    params:
        "--check_existing --pick --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --af --af_1kg --af_gnomad --max_af --pubmed --variant_class "
        # "--everything --check_existing --pick"
    log:
        "logs/variantCalling/pindel/{sample}.ann.log"
    threads:    8
    singularity:
        config["singularitys"]["vep"]
    shell:
        """(if [[ $(cat {input.vcf} | grep -v '^#' | wc -l) -eq 0 ]]; then mv {input.vcf} {output}
        else vep --vcf --no_stats -o {output} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --cache --merged --offline --fasta {input.fasta} {params} ; fi) &> {log}"""

rule filterPindel:
    input:
        vcf = "variantCalls/pindel/{sample}.pindel.ann.vcf"
    output:
        temp("variantCalls/pindel/{sample}.pindel.filt.vcf")
    log:
        "logs/variantCalling/pindel.{sample}.filt.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python3.6 /apps/bio/repos/somatic-twist/src/variantCalling/filter_vcf.py {input.vcf} {output}) &> {log}"

rule bgzipPindel:
    input:
        "variantCalls/pindel/{sample}.pindel.filt.vcf"
    output:
        "variantCalls/pindel/{sample}.pindel.filt.vcf.gz",
        "variantCalls/pindel/{sample}.pindel.filt.vcf.gz.tbi"
    log:
        "logs/variantCalling/pindel/{sample}.bgzip-index.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "( bgzip {input} && tabix {input}.gz ) &> {log}"
