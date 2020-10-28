#!/bin/bash -l
shopt -s extglob

rawdatafolder=$1  #Always with the last /
sequencerun=$2
sspath=$3
#Stand in folder where you want to run and get the results.
#somaticFolder=/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline
somaticFolder=/apps/bio/repos/somatic-twist
singularityFolder=/apps/bio/singularities/gms_hematology

today=$(date +'%y%m%d') #Until runfolder or should it always be rundate? And then have the runfolder with seq date
cp ${somaticFolder}/src/configfileSetup/configdefaults200317-gbg.yaml ${today}.config.yaml ## Copy configfile default to runfolder
for i in $(ls ${rawdatafolder}*R1_001.fastq.gz); do ##Add each sample in seq run to the config yaml file What happens to undet.?
  filename=$(echo $i | awk -F/ '{print $NF}')
  sample=${filename%%_S+([[:digit:]])_R1*}
  echo -e "  \"${sample}\": \"${i}\"" >> ${today}.config.yaml
done
echo "seqID:">> ${today}.config.yaml
echo -e "  sequencerun: \"${sequencerun}\"">> ${today}.config.yaml

export DRMAA_LIBRARY_PATH=/apps/univa-gridengine/lib/lx-amd64/libdrmaa.so.1.0
echo "DRMAA_LIBRARY_PATH: $DRMAA_LIBRARY_PATH"
#export SINGULARITY_BIND="/medstore,/seqstore,/apps"
#echo "SINGULARITY_BIND: $SINGULARITY_BIND"
export SGE_ROOT=/apps/univa-gridengine
echo "SGE_ROOT: $SGE_ROOT"
export SGE_CELL=default
echo "SGE_CELL: $SGE_CELL"
set +u;  PATH=$PATH:/apps/bio/software/singularity/bin; set -u
echo "PATH: $PATH"
unset LD_PRELOAD # förhindrar spam i singularities

mkdir -p ./fastq
mkdir -p ./logs/stdlogs
cp $sspath ./fastq/SampleSheetUsed.csv
my_cwd=$(pwd)

#unset PETASUITE_REFPATH
#Actually starting the pipeline with snakemake
#snakemake -p -j 999 --cluster "qsub -V -S /bin/bash -q {cluster.queue} -pe mpi {cluster.n} -e ./logs/stderr.log -o ./logs/stderr.log " \
snakemake -p -j 999 --drmaa " -V -q {cluster.queue} -pe mpi {cluster.n} -l excl=1 -o ${my_cwd}/logs/stdlogs/{cluster.stdout} -e ${my_cwd}/logs/stdlogs/{cluster.stderr}" \
--jobname "{cluster.name}_{jobid}" \
-s ${somaticFolder}/src/somaticPipeline.smk \
--configfile ${today}.config.yaml \
--use-singularity --singularity-prefix ${singularityFolder} \
--singularity-args "-e --bind /medstore --bind /seqstore --bind /apps " \
--latency-wait 60 \
--cluster-config ${somaticFolder}/cluster-config-gbg.json \
--rerun-incomplete # \ # Uncomment this backslash if running for DAG
#--dag | dot -Tsvg > dag.svg

#use this -e flag is for clean env in singularity
#--singularity-args "-e --bind /medstore --bind /seqstore --bind /apps " \

#not this
# --singularity-args "--bind /medstore --bind /seqstore --bind /apps " --latency-wait 60 \

#--director /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/ \  #Working dir
 #Configfile just created.
  #Snakemakefile
