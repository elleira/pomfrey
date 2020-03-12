
rule vardict:
    input:
        bam = "Results/{sample}/Data/{sample}-dedup.bam",  # differnet path sort of like: "{delivery}/bam/{sample}.bam"
        index = "Results/{sample}/Data/{sample}-dedup.bam.bai",
        ref = "/medstore/External_References/hs37d5/hs37d5.fa",
        bed =  config["bed"]["bedfile"]#"/gluster-storage-volume/projects/wp4/nobackup/workspace/somatic_dev/bedfiles/TST500C_manifest.bed"
    output:
        temp("variantCalls/callers/vardict/{sample}.vardict.vcf")
    params:
        af = "0.01"
    log:
        "logs/variantCalling/vardict/{sample}.log"
    threads:
        4
    singularity:
        config["singularitys"]["vardict"]
    shell:
        "(vardict-java -G {input.ref} -f {params.af} -th {threads} -N '{wildcards.sample}' -z -c 1 -S 2 -E 3 -b {input.bam} {input.bed} | "
        "teststrandbias.R | var2vcf_valid.pl -N '{wildcards.sample}' -E -f {params.af} > {output}) 2> {log}"
