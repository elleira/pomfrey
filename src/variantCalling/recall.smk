localrules: fixContig, sort_recall

rule recall:
    input:
        vcfs = expand("variantCalls/callers/{method}/{{sample}}.{method}.normalized.vcf.gz",  method=config["methods"]) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("variantCalls/callers/{method}/{{sample}}.{method}.normalized.vcf.gz.tbi",  method=config["methods"]),
        ref = "/medstore/External_References/hs37d5/hs37d5.fa"
    output:
        vcf = "variantCalls/recall/{sample}.{support}.unsorted.vcf.gz"
    params:
        support =  "{support}" ,
        order = "Vardict,Pisces,Freebayes,SNVer" #,Manta" #Make sure that the order is correct! Order of methods in configfile
    log:
        "logs/variantCalling/recall/{sample}.{support}.log"
    singularity:
        config["singularitys"]["recall"]
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle/somaticpipeline/src/singularity/bcbio-variation-recall-0.2.6-0.simg" Fungerar inte beh;ver bcftool!!
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/somatic_dev/bcbio-variation-recall.simg" #Dev
    shell: ##Remove filtered?? if so --nofiltered
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.order} {output.vcf} {input.ref} {input.vcfs}) 2> {log}"

rule fixContig:
    input:
        "variantCalls/recall/{sample}.{support}.unsorted.vcf.gz"
    output:
        temp("variantCalls/recall/{sample}.{support}.unsorted2.vcf")
    params: ## awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
        "\"##contig=<ID=1,length=249250621>\\n##contig=<ID=2,length=243199373>\\n##contig=<ID=3,length=198022430>\\n##contig=<ID=4,length=191154276>\\n##contig=<ID=5,length=180915260>\\n##contig=<ID=6,length=171115067>\\n##contig=<ID=7,length=159138663>\\n##contig=<ID=8,length=146364022>\\n##contig=<ID=9,length=141213431>\\n##contig=<ID=10,length=135534747>\\n##contig=<ID=11,length=135006516>\\n##contig=<ID=12,length=133851895>\\n##contig=<ID=13,length=115169878>\\n##contig=<ID=14,length=107349540>\\n##contig=<ID=15,length=102531392>\\n##contig=<ID=16,length=90354753>\\n##contig=<ID=17,length=81195210>\\n##contig=<ID=18,length=78077248>\\n##contig=<ID=19,length=59128983>\\n##contig=<ID=20,length=63025520>\\n##contig=<ID=21,length=48129895>\\n##contig=<ID=22,length=51304566>\\n##contig=<ID=X,length=155270560>\\n##contig=<ID=Y,length=59373566>\\n##contig=<ID=MT,length=16569>\\n\""
    log:
        "logs/variantCalling/recall/{sample}.{support}.fixContig.log"
    shell:
        """(zcat {input} | grep -v "^##contig" | awk '/^#CHROM/ {{ printf({params});}} {{print;}}' > {output} )&> {log}  """

rule sort_recall:
    input:
        "variantCalls/recall/{sample}.{support}.unsorted2.vcf"
    output:
        "variantCalls/recall/{sample}.{support}.vcf.gz"
    log:
        "logs/variantCalling/recall/{sample}.{support}.sort.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "( bcftools sort -o {output} -O z {input} )&> {log}"
