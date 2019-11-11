# somatic-twist

To run the pipeline start with snakemake using:
1. create config file from sample sheet configfile.yaml
2. download pipeline from git
3. module load slurm-drmaa
4. ` snakemake -p -j 32 -drmaa "-A wp4 -s -p core -t {cluster.time} -n {cluster.n} " -s src/somaticPipeline.smk --use-singularity --singularity-prefix /projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/src/singularity/ --singularity-args "--bind /data --bind /gluster-storage-volume " --latency-wait 5 --cluster-config src/cluster-som.json --rerun-incomplete --configfile configfile.yaml `
