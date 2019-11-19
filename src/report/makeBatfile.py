#!/bin/python
import sys

batfile = open(sys.argv[1], 'w')
# vcf_in=sys.argv[2]
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

with open(sys.argv[2], 'r') as vcfPASS:
    for line in vcfPASS:
        vcfRow = line.split()
        chr = vcfRow[0]
        pos = int(vcfRow[1])-1
        low = pos-padding
        high = pos+len(vcfRow[3])+padding

        batfile.write("\ngoto "+chr+":"+str(low)+"-"+str(high))
        if len(vcfRow[3]) > len(vcfRow[4]):  #Deletions
            sortPos = int(vcfRow[1])+1
        else:
            sortPos = int(vcfRow[1])
        batfile.write("\nsort "+sort+" "+str(sortPos)+"\n"+view+" "+bamfile.split("/")[-1])
        batfile.write("\nsnapshot "+chr+"_"+str(pos)+"_"+str(pos+len(vcfRow[3]))+"."+format)

batfile.write("\nexit")
batfile.close()
