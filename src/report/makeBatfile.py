#!/bin/python
import sys

batfile = open(sys.argv[1], 'w')
# vcf_in=sys.argv[2]
bamfile=sys.argv[3]
reffile = sys.argv[4]
bedfile = sys.argv[5]
outfolder = sys.argv[6]
padding = int(sys.argv[7])
sort = sys.argv[8]
view = sys.argv[9]
format = sys.argv[10]

batfile.write("new")
batfile.write("\ngenome "+reffile)
batfile.write("\nload "+bedfile+"\nload "+bamfile)
batfile.write("\nsnapshotDirectory "+outfolder)

with open(sys.argv[2], 'r') as vcfPASS:
    for line in vcfPASS:
        vcfRow = line.split()
        chr = vcfRow[0]
        pos = int(vcfRow[1])-1
        low = pos-padding
        high = pos+len(vcfRow[3])+padding

        batfile.write("\ngoto "+chr+":"+str(low)+"-"+str(high))
        batfile.write("\nsort "+sort+"\n"+view+" "+bamfile)
        batfile.write("\nsnapshot "+chr+"_"+str(pos)+"_"+str(pos+len(vcfRow[3]))+"."+format)

batfile.write("\nexit")
batfile.close()
