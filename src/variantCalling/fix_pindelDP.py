#!/bin/python3.6
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  #dosen't matter if bgziped or not. Automatically recognizes

new_header = vcf_in.header
# import pdb; pdb.set_trace()
new_header.info.add("DP","1","Integer","Sum of AD fields")

#start new vcf with the new_header
vcf_out = VariantFile(sys.argv[2], 'w', header=new_header)

for record in vcf_in.fetch():
    dp = sum(record.samples[0].get("AD"))
    record.info["DP"]=dp
    vcf_out.write(record)
