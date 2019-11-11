#!/bin/bash
import sys
from pysam import VariantFile

## Define sys.argvs
vcf_in=VariantFile(sys.argv[1])
vcf_out=VariantFile(sys.argv[2],'w', header=vcf_in.header)

for record in vcf_in.fetch():
    if record.filter.keys()==["PASS"]:
        if len(record.ref) > len(record.alts[0]):
            record.pos = record.pos-1 ##Blir inte ratt
        vcf_out.write(record)
