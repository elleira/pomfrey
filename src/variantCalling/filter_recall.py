#!/bin/python3.6
import sys
from pysam import VariantFile
import subprocess

vcf_in = VariantFile(sys.argv[1])
new_header = new_header=vcf_in.header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)
sv_out = sys.argv[2]+'.svtypeDEL.txt'
indelArteFile = sys.argv[3]

for record in vcf_in.fetch():
    # import pdb; pdb.set_trace()
    try:
        if record.info["SVTYPE"] == 'DEL':
            with open(sv_out, 'a+') as svtype_out:
                svtype_out.write(str(record))
    except KeyError:
        if len(record.ref) != len(record.alts[0]): #if InDel
            if "Mutect2" in record.info["CALLERS"] or "Vardict" in record.info["CALLERS"]: #Support by either Vardict or Manta, ok.
                ##Check if indel artefact
                # import pdb; pdb.set_trace()
                write = 1
                cmdIndelArte = 'grep -w '+str(record.pos)+' '+indelArteFile
                artefactLines = subprocess.run(cmdIndelArte, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                for artefactLine in artefactLines.split("\n"):
                    if artefactLine and record.ref == artefactLine.split()[2] and record.alts[0] == artefactLine.split()[3]: #if ref and alt is same as artefacts
                        write = 0 #Could exit for loop?
                if write == 1:
                    vcf_out.write(record)
        elif len(record.info["CALLERS"]) > 2: #if substitution demand 3/4 support
            vcf_out.write(record)
