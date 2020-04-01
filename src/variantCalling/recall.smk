localrules: fixContig, sort_recall,createMultiVcf

rule recall:
    input:
        vcfs = expand("variantCalls/callers/{method}/{{sample}}_{{seqID}}.{method}.normalized.vcf.gz",  method=config["methods"]) ,  # same order as methods in config!! make sure that is correct
        tabix = expand("variantCalls/callers/{method}/{{sample}}_{{seqID}}.{method}.normalized.vcf.gz.tbi",  method=config["methods"]),
        ref = config["reference"]["ref"]
    output:
        vcf = "variantCalls/recall/{sample}_{seqID}.unsorted.vcf.gz"
    params:
        support =  "3", #"{support}" ,
        order = "Vardict,Pisces,Freebayes,SNVer" #,Manta" #Make sure that the order is correct! Order of methods in configfile
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.log"
    singularity:
        config["singularitys"]["recall"]
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle/somaticpipeline/src/singularity/bcbio-variation-recall-0.2.6-0.simg" Fungerar inte beh;ver bcftool!!
        # "/gluster-storage-volume/projects/wp4/nobackup/workspace/somatic_dev/bcbio-variation-recall.simg" #Dev
    shell: ##Remove filtered?? if so --nofiltered
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.order} {output.vcf} {input.ref} {input.vcfs}) 2> {log}"

rule fixContig:
    input:
        "variantCalls/recall/{sample}_{seqID}.unsorted.vcf.gz"
    output:
        "variantCalls/recall/{sample}_{seqID}.unsorted2.vcf"
    params: ## awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai
        "\"##contig=<ID=1,length=249250621>\\n##contig=<ID=2,length=243199373>\\n##contig=<ID=3,length=198022430>\\n##contig=<ID=4,length=191154276>\\n##contig=<ID=5,length=180915260>\\n##contig=<ID=6,length=171115067>\\n##contig=<ID=7,length=159138663>\\n##contig=<ID=8,length=146364022>\\n##contig=<ID=9,length=141213431>\\n##contig=<ID=10,length=135534747>\\n##contig=<ID=11,length=135006516>\\n##contig=<ID=12,length=133851895>\\n##contig=<ID=13,length=115169878>\\n##contig=<ID=14,length=107349540>\\n##contig=<ID=15,length=102531392>\\n##contig=<ID=16,length=90354753>\\n##contig=<ID=17,length=81195210>\\n##contig=<ID=18,length=78077248>\\n##contig=<ID=19,length=59128983>\\n##contig=<ID=20,length=63025520>\\n##contig=<ID=21,length=48129895>\\n##contig=<ID=22,length=51304566>\\n##contig=<ID=X,length=155270560>\\n##contig=<ID=Y,length=59373566>\\n##contig=<ID=MT,length=16569>\\n\""
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.fixContig.log"
    shell:
        """(zcat {input} | grep -v "^##contig" | awk '/^#CHROM/ {{ printf({params});}} {{print;}}' > {output} ) &> {log}  """
rule sort_recall:
    input:
        "variantCalls/recall/{sample}_{seqID}.unsorted2.vcf" #multiAllelic.vcf"
    output:
        vcf = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz",
        tbi = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz.tbi"
    log:
        "logs/variantCalling/recall/{sample}_{seqID}.sort.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "( bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf} ) &> {log}"

# ##Add in multiallelic Variants
rule createMultiVcf:
     input:
         "variantCalls/callers/pisces/{sample}_{seqID}.pisces.normalized.vcf.gz" ##based on pisces!!!
     output:
         "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
     log:
         "logs/recall/{sample}_{seqID}.multiPASS.log"
     shell:
         """(zcat {input} | grep '^#' >{output} &&
         for pos in $(zcat {input} | grep -v '^#' | awk '{{print $2}}' | sort |uniq -d);do
         zcat {input}|grep $pos | grep PASS  >> {output} || true ; done) &> {log}"""

rule sort_multiPASS:
    input:
        "variantCalls/recall/{sample}_{seqID}.multiPASS.vcf"
    output:
        vcf = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz", #temp
        tbi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz.tbi" #temp
    log:
        "logs/recall/{sample}_{seqID}.multiPASS.sort.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools sort -o {output.vcf} -O z {input} && tabix {output.vcf}) &> {log}"

rule concatMulti:
    input:
        vcf = "variantCalls/recall/{sample}_{seqID}.notMulti.vcf.gz",
        multi = "variantCalls/recall/{sample}_{seqID}.multiPASS.sort.vcf.gz"
    output:
        "variantCalls/recall/{sample}_{seqID}.vcf.gz"
    params:
        "--allow-overlaps -d all -O z"
    log:
        "logs/recall/{sample}_{seqID}.concat.log"
    singularity:
        config["singularitys"]["bcftools"]
    shell:
        "(bcftools concat {params} -o {output} {input.vcf} {input.multi}) &> {log}"
