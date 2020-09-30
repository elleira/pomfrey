#!/usr/bin/env bash
module load slurm-drmaa/1.0.7
shopt -s extglob
#Stand in folder where you want to run and get the results!
rawdataFolder=$1  #To the demultiplexed lanes merged fastq.gz files
sequencerun=$2    #Sequence ID

somaticFolder=/projects/wp2/nobackup/Twist_Myeloid/Bin/pomfrey/
today=$(date +'%y%m%d') #Until runfolder or should it always be rundate? And then have the runfolder with seq date
# Merge samples perlane into same.
for i in $(cat fastq/SampleSheetUsed.csv | grep 'Sample_ID' -A20 |awk 'FS=","{print $1}' | grep -v 'Sample_ID'); do ## If more than 20 samples will not work!!
   sample=$(ls fastq-perLane/${i}*L001_R1_001.fastq.gz|cut -d '/' -f 2 |cut -d '_' -f1-2); #To get the SX part of samplename, do we rally need it?
   echo ${sample}
   cat fastq-perLane/${i}*R1_001.fastq.gz >${rawdataFolder}/${sample}_R1_001.fastq.gz
   cat fastq-perLane/${i}*R2_001.fastq.gz >${rawdataFolder}/${sample}_R2_001.fastq.gz
 done
# cp fastq-perlane/SampleSheetUsed.csv ${rawdataFolder}/SampleSheetUsed.csv
# Add samples to snakemake configfile
python3 ${somaticFolder}/src/configfileSetup/addSamples2config.py ${rawdataFolder} ${somaticFolder}/src/configfileSetup/configdefaults200417.yaml ${sequencerun}
# Collect statistics to KIBANA, needs outfolder
python3 ${somaticFolder}/src/configfileSetup/collectStatistics.py ${rawdataFolder}/SampleSheetUsed.csv ${sequencerun}.json

#Actually starting the pipeline with snakemake
snakemake -p -j 32 --drmaa "-A wp4 -s -p core -t {cluster.time} -n {cluster.n} " \
-s ${somaticFolder}/src/somaticPipeline.smk \
--configfile ${sequencerun}_config.yaml \
--use-singularity  \
--singularity-args "--bind /data --bind /projects " --latency-wait 5 \
--cluster-config ${somaticFolder}/cluster-config.json \
--rerun-incomplete

#--director /gluster-storage-volume/projects/wp4/nobackup/workspace/arielle_test/somaticpipeline/ \  #Working dir
 #Configfile just created.
  #Snakemakefile
