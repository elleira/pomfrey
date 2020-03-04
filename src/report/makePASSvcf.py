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
                cmdCosmic = 'grep -w '+synoCosmicId+' /apps/bio/singularities/gms_hematology/bedfiles/COSMIC_v90_hemato_counts.txt | cut -f 16 '
                synoCosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                if len(synoCosmicNew) == 0:
                    synoCosmicNew = 0
                synoCosmicN += int(synoCosmicNew)

    if record.filter.keys()==["PASS"] or synoCosmicN != 0 :
        # if len(record.ref) > len(record.alts[0]): Deletion but need to remove first base as well
        #     record.pos = record.pos-1
        if record.info["AF"][0] >= 0.03:
            cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
            artLines = subprocess.run(cmdArt, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip() ##What happens if two hits?
            artefact_variant = 0

            for artLine in artLines.split("\n"):
                if artLine and record.ref == artLine.split()[2] and record.alts[0] == artLine.split()[3]: #if pos exists and match in artefact file.
                    #Artefact do not print
                    artefact_variant = 1
                    continue

            if artefact_variant == 0:
                cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
                germLines = subprocess.run(cmdGerm, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                germline_variant = 0
                for germLine in germLines.split("\n"):
                    if germLine and record.ref == germLine.split("\t")[2] and record.alts[0] == germLine.split("\t")[3]: #if exists in germline file
                        #Germline match, do nothing
                        germline_variant = 1
                        continue
                if germline_variant==0:
                    vcf_out.write(record)
