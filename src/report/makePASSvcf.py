#!/bin/python3.6
import sys
import subprocess
from pysam import VariantFile

## Define sys.argvs
vcf_in = VariantFile(sys.argv[1])
artefactFile = sys.argv[2]
germlineFile = sys.argv[3]
vcf_out = VariantFile(sys.argv[4],'w', header=vcf_in.header)

for record in vcf_in.fetch():
    synoCosmicN = 0
    if record.filter.keys()==["Syno"]: #If only syno! No popAF.    any(x in "Syno" for x in record.filter.keys()):
        csq = record.info["CSQ"][0]
        synoCosmicVepList = [cosmic for cosmic in csq.split("|")[17].split("&") if cosmic.startswith('CO')] #Get all cosmicID in list
        if len(synoCosmicVepList) != 0:
            for synoCosmicId in synoCosmicVepList:
                cmdCosmic = 'grep -w '+synoCosmicId+' /gluster-storage-volume/data/ref_data/COSMIC/COSMIC_v90_hemato_counts.txt | cut -f 16 '
                synoCosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                if len(synoCosmicNew) == 0:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)

    if record.filter.keys()==["PASS"] or synoCosmicN != 0 :
        # if len(record.ref) > len(record.alts[0]): Deletion but need to remove first base as well
        #     record.pos = record.pos-1
        cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
        artLine = subprocess.run(cmdArt, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip() ##What happens if two hits?
        if artLine and record.ref == artLine.split("\t")[2] and record.alts[0] == artLine.split("\t")[3]: #if pos exists and match in artefact file.
            #Artefact do not print
            continue
        else:
            # Germline /gluster-storage-volume/projects/wp2/nobackup/Twist_Myeloid/Artefact_files/Low_VAF_SNVs.txt
            cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
            germLine = subprocess.run(cmdGerm, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
            if germLine and record.ref == germLine.split("\t")[2] and record.alts[0] == germLine.split("\t")[3]: #if exists in germline file
                #Germline match, do nothing
                continue
            else:
                vcf_out.write(record)
