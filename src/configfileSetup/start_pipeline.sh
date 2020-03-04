#!/usr/bin/env bash
module load slurm-drmaa/1.0.7
shopt -s extglob

rawdatafolder=$1  #Always with the last /
sequnecerun=$2
#Stand in folder where you want to run and get the results.
somaticFolder=/gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline
today=$(date +'%y%m%d') #Until runfolder or should it always be rundate? And then have the runfolder with seq date
cp ${somaticFolder}/src/configfileSetup/configdefaults200331.yaml ${today}.config.yaml ## Copy configfile default to runfolder
for i in $(ls ${rawdatafolder}*R1_001.fastq.gz); do ##Should use sample sheet instead!
  filename=$(echo $i | awk -F/ '{print $NF}')
  sample=${filename%%_S+([[:digit:]])_R1*}
  echo -e "  \"${sample}\": \"${i}\"" >> ${today}.config.yaml
done
echo "seqID:">> ${today}.config.yaml
echo -e "  sequencerun: \"${sequnecerun}\"">> ${today}.config.yaml


#Actually starting the pipeline with snakemake
snakemake -p -j 32 --drmaa "-A wp4 -s -p core -n {cluster.n} " \
-s ${somaticFolder}/src/somaticPipeline.smk \
--configfile ${today}.config.yaml \
--use-singularity --singularity-prefix ${somaticFolder}/src/singularity/ \
--singularity-args "--bind /data --bind /gluster-storage-volume " --latency-wait 5 \
--cluster-config ${somaticFolder}/cluster-config.json \
--rerun-incomplete

#--director /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/ \  #Working dir
 #Configfile just created.
  #Snakemakefile
