#!/bin/python3.6
import sys
import re
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  #dosen't matter if bgziped or not. Automatically recognizes
method = re.search('callers/(.+?)/',sys.argv[1]).group(1)  ##The folder after callers/
# Add new filter descriptions to new header
new_header = vcf_in.header
# import pdb; pdb.set_trace()

if method == "freebayes" or method == "pisces": #Byta description on AF for freebayes?
    sample = vcf_in.header.samples[0]

if method == "pisces":
    new_header.info.add("AF","A","Float","DescriptionDescription")

if  method == "snver":
    new_header.info.add("AF","A","Float","Allel count divided on depth, crude")

#start new vcf with the new_header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)

for record in vcf_in.fetch():
    if method == "freebayes":
        ad = record.samples[sample].get("AD")
        af=[]
        af = [ad[1]/(ad[0]+ad[1])]
        if len(ad) > 2:
            for item in ad[2:]:
                af.append(item/sum(ad))
    if method == "pisces":
        af = record.samples[sample].get("VF")
    if method == "snver":
        dp=record.info["DP"]
        ac=record.info["AC"]
        af=ac/dp

    record.info["AF"]=af

    vcf_out.write(record)
