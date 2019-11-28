
include:    "freebayes.smk"

include:    "lofreq_T.smk"

include:    "snver.smk" #since based on samtools cannot handle high base qualities, everything becomes N

include:    "vardict_T.smk"

include:    "pisces.smk"

#include:    "manta_T.smk"
rule fixAF:
    input:
        "variantCalls/callers/{method}/{sample}.{method}.weirdAF.vcf"
    output:
        temp("variantCalls/callers/{method}/{sample}.{method}.vcf")
    log:
        "logs/fixAF/{method}/{sample}.log"
    singularity:
        config["singularitys"]["python"]
    shell:
        "(python /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/variantCalling/fix_af.py {input} {output}) &> {log}"


include:    "bgzips.smk"

include:    "normalize.smk"


include:    "recall.smk"

include:    "vep.smk"
