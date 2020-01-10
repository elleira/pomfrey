#!/bin/python3.6
import sys
from pysam import VariantFile

batfile = open(sys.argv[1], 'w')
vcf_in=VariantFile(sys.argv[2],'r')
#indel_in=sys.argv[3]
bamfile=sys.argv[4]
reffile = sys.argv[5]
bedfile = sys.argv[6]
outfolder = sys.argv[7]
padding = int(sys.argv[8])
sort = sys.argv[9]
view = sys.argv[10]
format = sys.argv[11]

batfile.write("new")
batfile.write("\ngenome "+reffile)
batfile.write("\nload "+bedfile+"\nload "+bamfile)
batfile.write("\nsnapshotDirectory "+outfolder)

for record in vcf_in.fetch():
        # vcfRow = line.split()
        chr = record.contig
        pos = record.pos-1
        low = pos-padding
        high = pos+len(record.alts[0])+padding

        batfile.write("\ngoto "+chr+":"+str(low)+"-"+str(high))
        if len(record.ref) > len(record.alts[0]):  #Deletions
            sortPos = record.pos+1
        else:
            sortPos = record.pos
        batfile.write("\nsort "+sort+" "+str(sortPos)+"\n"+view+" "+bamfile.split("/")[-1])
        batfile.write("\nsnapshot "+chr+"_"+str(pos)+"_"+str(pos+len(record.alts[0]))+"."+format)

batfile.write("\nexit")
batfile.close()
