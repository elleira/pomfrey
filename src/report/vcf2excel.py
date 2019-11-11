#!/bin/bash
import sys
import csv
from pysam import VariantFile
import xlsxwriter ## need to download? Can I use same as pysam singularity?

## Define sys.argvs
vcf_snv=VariantFile(sys.argv[1])
vcf_indel=VariantFile(sys.argv[2])
cartool=sys.argv[3]
minCov=sys.argv[4]  ##grep thresholds ../somaticpipeline/qc/cartool/10855-17_Log.csv | cut -d ',' -f2
bedfile=sys.argv[5]
output = sys.argv[6]
 ## Create execl file and sheets.
workbook = xlsxwriter.Workbook(output)
worksheetSNV = workbook.add_worksheet('SNVs')
worksheetIndel = workbook.add_worksheet('InDel') #.... sys.argv[2]
worksheetCov = workbook.add_worksheet('Coverage') #... sys.argv[3]
worksheetVersions = workbook.add_worksheet('Version Log')
## Define formats to be used.
boldFormat = workbook.add_format({'bold': 1})
lowCovFormat = workbook.add_format({'bg_color': '#f7a19a'})
# headerFormat = workbook.add_format()

## Define sample based on annotated vcf
sample = list(vcf_snv.header.samples)[0]

######## Create low cov dict #########

with open(cartool) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV) ##Skip header
    lowCovDict={}
    for line in readCSV:
        posList=list(range(int(line[2]),int(line[3])))
        if line[1] in lowCovDict.keys():
            lowCovDict[line[1]] = lowCovDict[line[1]] + posList
        else:
            lowCovDict[line[1]]=posList
######################################

######### SNV sheet ##################

widthC = max(len(str(record.pos)) for record in vcf_snv.fetch())
worksheetSNV.set_column(1,3, widthC+2) #Width for position column
worksheetSNV.set_column(1,11,9) #Width for MaxPopAf column.

for x in vcf_snv.header.records:
    if (x.key=='reference'):
        refV = x.value
    if (x.key=='VEP'):
        vepline = x.value

## Headers before variants
worksheetSNV.write('A1', 'Sample: '+str(sample), boldFormat)  #headerFormat)
worksheetSNV.write('A3', 'Reference used: '+str(refV))
# worksheetSNV.write('A4', 'Callers used: VardictJava v.1,6 ? Dubbelkolla!, Pisces 5.2.11, Freebayes v1.1.0, LoFreq v.2.1.3.1, SNVer v.2.1.3.1')
worksheetSNV.write('A5', 'VEP: '+vepline)
# Location of raw vcf?
# Todays date?
# Other info...

## Variant table
tableheading = ['Gene','Chr','Pos','Ref','Alt', 'AF', 'DP', 'COSMIC ids', 'Clinical significance', 'dbSNP','Max popAF','Max Pop']
worksheetSNV.write_row('A7', tableheading, boldFormat) #1 index
row = 7 #0 index
col=0

for record in vcf_snv.fetch():
    if record.filter.keys()==["PASS"]:

        if len(record.info["AF"]) == 1:
            af=record.info["AF"][0]
        else:
            print(record.info["AF"])
            sys.exit() ##Fix!!

        if len(record.alts) == 1:
            alt=record.alts[0]
        else:
            print(record.alts)
            sys.exit() ##Fix!!

        cosmicList = []
        geneList = []
        maxPopAfList = []
        maxPopList = []
        cliSigList = []
        rsList = []
        # import pdb; pdb.set_trace()
        for csq in record.info["CSQ"]:
            geneList.append(csq.split("|")[3])
            existing = csq.split("|")[17] ##What about the rest?? What is tmp_esp??
            cliSigList.append(csq.split("|")[62] )
            cosmicList += [cosm for cosm in existing.split("&") if cosm.startswith('COSM')]
            rsList += [rs for rs in existing.split("&") if rs.startswith('rs')]

        geneList = list(dict.fromkeys(geneList))
        gene = ', '.join(geneList)
        cosmicList = list(dict.fromkeys(cosmicList))
        cosmic = ', '.join(cosmicList)
        cliSigList = list(dict.fromkeys(cliSigList))
        clinical = ', '.join(cliSigList)
        rsList = list(dict.fromkeys(rsList))
        rs = ', '.join(rsList)

        maxPopAf = record.info["CSQ"][0].split("|")[60]
        if len(maxPopAf) > 1:
            maxPopAf = float(maxPopAf)
        maxPop = record.info["CSQ"][0].split("|")[61]

        # Should pos be shifted if del??
        snv = [gene, record.contig, record.pos, record.ref, alt, af, record.info["DP"], cosmic, clinical, rs, maxPopAf, maxPop]
        if record.contig in lowCovDict.keys() and record.pos in lowCovDict[record.contig]:
            worksheetSNV.write_row(row, col, snv, lowCovFormat)
        else:
            worksheetSNV.write_row(row, col, snv)
        row += 1
worksheetSNV.write(row+2,0,'Start position coverage <= '+str(minCov)+'x.',lowCovFormat)
########################################

######### Indel sheet ##################
# Add genes as info before the actual table. Just use bed table as input? Sort uniq

worksheetIndel.write('A1', 'Sample: '+str(sample), boldFormat)
with open(bedfile) as bed:
    genesDup = [line.split("\t")[3].strip() for line in bed]
    genes = set(genesDup)



for x in vcf_indel.header.records:
    if (x.key=='reference'):
        refI = x.value

worksheetIndel.write('A3', 'Reference used: '+str(refI))
worksheetIndel.write('A4', 'Genes looked at: '+str(genes))

tableheading = ['Chr','Start','End','SV length','AD','Ref','Alt']
worksheetIndel.write_row('A6',tableheading, boldFormat) #1 index
row = 6 #0 index
col=0

for indel in vcf_indel.fetch():
    if indel.filter.keys()==["PASS"]:
        # import pdb; pdb.set_trace()
        svlen = indel.info["SVLEN"]

        ads = indel.samples[sample]["AD"]
        ad = ', '.join([str(i) for i in ads])

        if len(indel.alts) == 1:
            alt=indel.alts[0]
        else:
            print(indel.alts)
            sys.exit() ##Fix!!

        indelRow = [indel.contig,indel.pos, indel.stop, svlen, ad, indel.ref, alt] ##Add gene, how?
        worksheetIndel.write_row(row,col,indelRow)
        row += 1

##########################################

######### Coverage sheet #################

## Heading in sheet
worksheetCov.write('A1','CarTool',boldFormat)
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetCov.write('A3', description)
covHeadings = ['Region Name','Chr','Start','Stop','Mean Coverage','Length of Region']
worksheetCov.write_row('A5',covHeadings,boldFormat) ## 1 index
row = 5 ## 0 index

with open(cartool) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    ## Skip header

    for line in readCSV:
        while not line[-1]: ## Remove empty fields at the end of line.
            line.pop()
        for i in range(1, int((len(line)-1)/5+1)): ##For each region/ with lower cov, create a new row.
            start = 1+5*(i-1)
            end = 1+5*i
            worksheetCov.write_row(row,col,[line[0]]+line[start:end])
            row += 1

##########################################

######### Prog Version sheet #################
worksheetVersions.write('A1','Version Log', boldFormat)
worksheetVersions.write('A3', 'Variant calling reference used: '+str(refV))
worksheetVersions.write('A4', 'Pindel reference used: '+str(refI))
worksheetVersions.write('A6', 'Containers used: ', boldFormat)


with open('containers.txt') as file:
    singularitys = [line.strip() for line in file]
    singularitys.pop() ##Last slurm is always the makeContainersList rule.
    row = 6
    col = 0
    for singularity in singularitys:
        worksheetVersions.write_row(row,col,[singularity])
        row += 1



workbook.close()
