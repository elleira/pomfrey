
rule indelReali: #ca 1h for twist sample
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample], #Outfile from bwa-mem mapped/bwa/{sample}.bam
        ref = "refs/hg19.with.mt.fasta",
        intervals = "realigBam/TargetCreator-Mills.intervals" #Full better since not mills not in order! Just skipp targetIntervals then
    output:
        bam = "mapped/indelReali/{sample}-realign.bam" #{bam}/{sample}.realign.bam
    log:
        "logs/indel/{sample}-realin.log"
    singularity:
        "/projects/wp1/nobackup/ngs/utveckling/software/Pipeline/singularity/swift.simg"
    shell:
        "(gatk3 -T IndelRealigner -R {input.ref} -I {input.bam} -targetIntervals {input.intervals} -o {output}) 2> {log}"


#only needed once as long as known indels does not change..? Maybe redo it with different panels? Also need to sort bam if using Mills!
# rule targetCreator:
#     input:
#         bam = lambda wildcards: config["samples"][wildcards.sample],
#         ref = "refs/hg19.with.mt.fasta"
#     output:
#         intervals = "{sample}-forIndelRealigner.intervals"
#     log:
#         "logs/indel/{sample}-target.log"
#     threads: 4
#     singularity:
#         "/projects/wp1/nobackup/ngs/utveckling/software/Pipeline/singularity/swift.simg"
#     shell:
#         "(gatk3 -T RealignerTargetCreator -nt {threads} -R {input.ref} -I {input.bam} -o {output}) 2> {log}"
