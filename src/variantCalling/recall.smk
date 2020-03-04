localrules: fixContig, sort_recall,createMultiVcf

rule recall:
    input:
        vcfs = expand("variantCalls/callers/{method}/{{sample}}.{method}.normalized.vcf.gz",  method=config["methods"]) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("variantCalls/callers/{method}/{{sample}}.{method}.normalized.vcf.gz.tbi",  method=config["methods"]),
        ref = config["reference"]["ref"]
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
        "variantCalls/recall/{sample}.{support}.unsorted2.vcf"
    params: ## awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
        "\"##contig=<ID=chr1,length=249250621>\\n##contig=<ID=chr2,length=243199373>\\n##contig=<ID=chr3,length=198022430>\\n##contig=<ID=chr4,length=191154276>\\n##contig=<ID=chr5,length=180915260>\\n##contig=<ID=chr6,length=171115067>\\n##contig=<ID=chr7,length=159138663>\\n##contig=<ID=chr8,length=146364022>\\n##contig=<ID=chr9,length=141213431>\\n##contig=<ID=chr10,length=135534747>\\n##contig=<ID=chr11,length=135006516>\\n##contig=<ID=chr12,length=133851895>\\n##contig=<ID=chr13,length=115169878>\\n##contig=<ID=chr14,length=107349540>\\n##contig=<ID=chr15,length=102531392>\\n##contig=<ID=chr16,length=90354753>\\n##contig=<ID=chr17,length=81195210>\\n##contig=<ID=chr18,length=78077248>\\n##contig=<ID=chr19,length=59128983>\\n##contig=<ID=chr20,length=63025520>\\n##contig=<ID=chr21,length=48129895>\\n##contig=<ID=chr22,length=51304566>\\n##contig=<ID=chrX,length=155270560>\\n##contig=<ID=chrY,length=59373566>\\n##contig=<ID=chrM,length=16571>\\n\""
    log:
        "logs/variantCalling/recall/{sample}.{support}.fixContig.log"
    shell:
        """(zcat {input} | grep -v "^##contig" | awk '/^#CHROM/ {{ printf({params});}} {{print;}}' > {output} ) &> {log}  """
rule sort_recall:
    input:
        "variantCalls/recall/{sample}.{support}.unsorted2.vcf" #multiAllelic.vcf"
    output:
        vcf = "variantCalls/recall/{sample}.{support}.notMulti.vcf.gz",
        tbi = "variantCalls/recall/{sample}.{support}.notMulti.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}.{support}.sort.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "( bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf} ) &> {log}"

# ##Add in multiallelic Variants
rule createMultiVcf:
     input:
         "variantCalls/callers/pisces/{sample}.pisces.normalized.vcf.gz" ##based on pisces!!!
     output:
         "variantCalls/recall/{sample}.{support}.multiPASS.vcf"
     log:
         "logs/recall/{sample}.{support}.multiPASS.log"
     shell:
         """(zcat {input} | grep '^#' >{output} &&
         for pos in $(zcat {input} | grep -v '^#' | awk '{{print $2}}' | sort |uniq -d);do
         zcat {input}|grep $pos | grep PASS  >> {output} || true ; done) &> {log}"""

rule sort_multiPASS:
    input:
        "variantCalls/recall/{sample}.{support}.multiPASS.vcf"
    output:
        vcf = "variantCalls/recall/{sample}.{support}.multiPASS.sort.vcf.gz", #temp
        tbi = "variantCalls/recall/{sample}.{support}.multiPASS.sort.vcf.gz.tbi" #temp
    log:
        "logs/recall/{sample}.{support}.multiPASS.sort.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf}) &> {log}"

rule concatMulti:
    input:
        vcf = "variantCalls/recall/{sample}.{support}.notMulti.vcf.gz",
        multi = "variantCalls/recall/{sample}.{support}.multiPASS.sort.vcf.gz"
    output:
        "variantCalls/recall/{sample}.{support}.vcf.gz"
    params:
        "--allow-overlaps -d all -O z"
    log:
        "logs/recall/{sample}.{support}.concat.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.vcf} {input.multi}) &> {log}"
