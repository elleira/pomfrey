#!/bin/bash -l
#module load slurm-drmaa/1.0.7
#module load singularity
shopt -s extglob

rawdatafolder=$1  #Always with the last /
sequencerun=$2
#Stand in folder where you want to run and get the results.
#somaticFolder=/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline
somaticFolder=/apps/bio/repos/somatic-twist
singularityFolder=/apps/bio/singularities/gms_hematology

date=$(date +'%y%m%d') #Until runfolder or should it always be rundate? And then have the runfolder with seq date
cp ${somaticFolder}/src/configfileSetup/configdefaults200217-gbg.yaml ${date}.config.yaml ## Copy configfile default to runfolder
for i in $(ls ${rawdatafolder}*R1_001.fastq.gz); do ##Add each sample in seq run to the config yaml file What happens to undet.?
  filename=$(echo $i | awk -F/ '{print $NF}')
  sample=${filename%%_S+([[:digit:]])_R1*}
  echo -e "  \"${sample}_${sequencerun}\": \"${i}\"" >> ${date}.config.yaml
done
echo "dates:">> ${date}.config.yaml
echo -e "  sequencerun: \"${sequencerun}\"">> ${date}.config.yaml

export DRMAA_LIBRARY_PATH=/apps/univa-gridengine/lib/lx-amd64/libdrmaa.so.1.0
echo "DRMAA_LIBRARY_PATH: $DRMAA_LIBRARY_PATH"
export SINGULARITY_BIND="/medstore,/seqstore,/apps"
echo "SINGULARITY_BIND: $SINGULARITY_BIND"
export SGE_ROOT=/apps/univa-gridengine
echo "SGE_ROOT: $SGE_ROOT"
export SGE_CELL=default
echo "SGE_CELL: $SGE_CELL"
set +u;  PATH=$PATH:/apps/bio/software/singularity/bin; set -u
echo "PATH: $PATH"
#Actually starting the pipeline with snakemake
#snakemake -p -j 999 --cluster "qsub -V -S /bin/bash -q {cluster.queue} -pe mpi {cluster.n} -e ./logs/stderr.log -o ./logs/stderr.log " \
snakemake -p -j 999 --drmaa " -V -q {cluster.queue} -pe mpi {cluster.n} -e {cluster.error} -o {cluster.output} -l excl=1 --jobname {cluster.name} " \
-s ${somaticFolder}/src/somaticPipeline.smk \
--configfile ${date}.config.yaml \
--use-singularity --singularity-prefix ${singularityFolder} \
--singularity-args "--bind /medstore --bind /seqstore --bind /apps " --latency-wait 60 \
--cluster-config ${somaticFolder}/cluster-config.json \
--rerun-incomplete # \ # Uncomment this if running for DAG
#--dag | dot -Tsvg > dag.svg

#--director /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/ \  #Working dir
 #Configfile just created.
  #Snakemakefile
