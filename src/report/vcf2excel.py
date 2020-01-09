#!/bin/bash
import sys
import csv
from pysam import VariantFile
import xlsxwriter ## need to download? Can I use same as pysam singularity?
from datetime import date
import subprocess

## Define sys.argvs
vcf_snv = VariantFile(sys.argv[1])
vcf_indel = VariantFile(sys.argv[2])
cartool = sys.argv[3]
minCov = sys.argv[4]  ##grep thresholds ../somaticpipeline/qc/cartool/10855-17_Log.csv | cut -d ',' -f2
bedfile = sys.argv[5]
hotspotFile = sys.argv[6]
artefactFile = sys.argv[7]
germlineFile = sys.argv[8]
output = sys.argv[9]


 ## Create execl file and sheets.
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')
worksheetSNV = workbook.add_worksheet('SNVs')
worksheetIndel = workbook.add_worksheet('InDel') #.... sys.argv[2]
worksheetCov = workbook.add_worksheet('Coverage') #... sys.argv[3]
worksheetHotspot = workbook.add_worksheet('Hotspot')
worksheetIVA = workbook.add_worksheet('IVA')
worksheetVersions = workbook.add_worksheet('Version')
## Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

greenFormat = workbook.add_format({'bg_color': '#85e085', 'font_color': '#b3b3b3'})
orangeFormat = workbook.add_format({'bg_color': '#ffd280', 'font_color': '#b3b3b3'}) #font color gray 70%
green_italicFormat = workbook.add_format({'bg_color': '#85e085', 'font_color': '#b3b3b3', 'italic': 'True'})
orange_italicFormat = workbook.add_format({'bg_color': '#ffd280', 'font_color': '#b3b3b3', 'italic': 'True'}) #font color gray 70%

## Define sample based on annotated vcf
sample = list(vcf_snv.header.samples)[0]
today=date.today()
emptyList=['','','','','','']
######## Create low cov dict #########
#
# with open(cartool,'r') as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')
#     next(readCSV) ##Skip header
#     lowCovDict={}
#     for line in readCSV:
#         posList=list(range(int(line[2]),int(line[3])))
#         if line[1] in lowCovDict.keys():
#             lowCovDict[line[1]] = lowCovDict[line[1]] + posList
#         else:
#             lowCovDict[line[1]]=posList
######################################

########## Hotspot sheet (5)#############
worksheetHotspot.write('A1', 'Hotspot Coverage', headingFormat)  #headerFormat)
worksheetHotspot.write('A3', 'Sample: '+str(sample))
worksheetHotspot.set_column(1,2,10)
worksheetHotspot.write_row('A5',['Chr', 'Pos', 'Depth', 'Gene'], tableHeadFormat)

lowPos = 0
row=5
with open(cartool.replace("_MeanCoverageShortList.csv", "_coverageShort.tsv"),'r') as covFile:
    for dpLine in covFile:
        cov = dpLine.split("\t")
        chrCov = cov[0]
        posCov = cov[1]
        dp = cov[2].rstrip()
        with open(hotspotFile, "r") as hotspotBed:
            for hotspotLine in hotspotBed:
                hotspot = hotspotLine.split("\t")
                chrHS = hotspot[0]
                lowHS = hotspot[1]
                # highHS = hotspot[2]
                if chrCov == chrHS and lowHS == posCov: #Start with just one pos, other if highHS and lowHS not same do something else.  or cov[1] == highHS ):
                    hotspotTable = [chrCov,posCov, dp, hotspot[3].rstrip()]
                    worksheetHotspot.write_row(row,0,hotspotTable)
                    row += 1
                    if int(dp) <= 500 :
                        lowPos += 1

############################################


######### Overview sheet (1) #################

worksheetOver.write(0,0, sample, headingFormat)
worksheetOver.write(1,0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(2,0, emptyList,lineFormat)

worksheetOver.write(3,0, "Created by: ")
worksheetOver.write(3,4, "Valid from: ")
worksheetOver.write(4,0, "Signed by: ")
worksheetOver.write(4,4, "Document nr: ")
worksheetOver.write_row(5,0,emptyList,lineFormat)

worksheetOver.write(6,0,"Sheets:", tableHeadFormat)
worksheetOver.write_url(7,0,"internal:'SNVs'!A1", string='Variants analysis')
worksheetOver.write_url(8,0,"internal:'Indel'!A1", string = 'Indel variants')
worksheetOver.write_url(9,0,"internal:'Coverage'!A1", string = 'Positions with coverage lower than 100x')
worksheetOver.write_url(10,0,"internal:'Hotspot'!A1", string = 'Coverage of hotspot positions')
worksheetOver.write_url(11,0,"internal:'Version'!A1", string = 'Version Log')
worksheetOver.write_row(13,0,emptyList,lineFormat)

if lowPos == 0:
    worksheetOver.write(16,0,'Number of positions from the hotspot list not covered by at least 500x: ')
    worksheetOver.write(17,0, str(lowPos))
else:
    worksheetOver.write(16,0,'Number of positions from the hotspot list not covered by at least 500x: ')
    worksheetOver.write(17,0, str(lowPos), redFormat)
    worksheetOver.write_url(18,0,"internal:'Hotspot'!A1" ,string = 'For more detailed list see hotspotsheet ')

##Added after CARTools sheet done
# worksheetOver.write(19,0,'Number of regions not covered by at least 100x: ')
# worksheetOver.write(20,0, str(lowRegions))
worksheetOver.write(22,0,'Hotspotlist: '+hotspotFile)
worksheetOver.write(23,0,'Artefact file: '+artefactFile)
worksheetOver.write(24,0,'Germline file: '+germlineFile)

######################################

######### SNV sheet (2) ##################

worksheetSNV.set_column(1,3,10) #Width for position column
# worksheetSNV.set_column(1,11,9) #Width for MaxPopAf column.

for x in vcf_snv.header.records:
    if (x.key=='reference'):
        refV = x.value
    if (x.key=='VEP'):
        vepline = x.value

## Headers before variants
worksheetSNV.write('A1', 'Variants found', headingFormat)
worksheetSNV.write('A3', 'Sample: '+str(sample))
worksheetSNV.write('A4', 'Reference used: '+str(refV))
# worksheetSNV.write('A4', 'Callers used: VardictJava v.1,6 ? Dubbelkolla!, Pisces 5.2.11, Freebayes v1.1.0, LoFreq v.2.1.3.1, SNVer v.2.1.3.1')
worksheetSNV.write('A6', 'VEP: '+vepline ) #, textwrapFormat)
worksheetSNV.write('A8', 'The following filters were applied: ')
worksheetSNV.write('B9','Coverage >= 100x')
worksheetSNV.write('B10','Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetSNV.write('B11','Biotype is protein coding')
worksheetSNV.write('B12','Consequence not deemed relevant')

worksheetSNV.write('A14','Coverage below 500x', italicFormat)
worksheetSNV.write('A15','Variant in artefact list ', orangeFormat)
worksheetSNV.write('A16','Variant likely germline', greenFormat)

## Variant table
tableheading = ['Gene','Chr','Pos','Ref','Alt', 'AF', 'DP', 'Canonical Transcript','Mutation cds', 'Consequence','COSMIC ids on position','N COSMIC Hemato hits on position','Clinical significance', 'dbSNP','Max popAF','Max Pop']
worksheetSNV.write_row('A18', tableheading, tableHeadFormat) #1 index
row = 18 #0 index
col=0
white=[]
green=[]
orange=[]

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

  #       cosmicList = []
  #       geneList = []
  #       maxPopAfList = []
  #       maxPopList = []
  #       cliSigList = []
  #       rsList = []
  # ##Change all with the new pick!
  #       for csq in record.info["CSQ"]:
  #           geneList.append(csq.split("|")[3])
  #           existing = csq.split("|")[17] ##What about the rest?? What is tmp_esp??
  #           cliSigList.append(csq.split("|")[62] )
  #           cosmicList += [cosm for cosm in existing.split("&") if cosm.startswith('COSM')]
  #           rsList += [rs for rs in existing.split("&") if rs.startswith('rs')]
        #
        # geneList = list(dict.fromkeys(geneList))
        # gene = ', '.join(geneList)
        # cosmicList = list(dict.fromkeys(cosmicList))
        # cosmic = ', '.join(cosmicList)
        # cliSigList = list(dict.fromkeys(cliSigList))
        # clinical = ', '.join(cliSigList)
        # rsList = list(dict.fromkeys(rsList))
        # rs = ', '.join(rsList)

        csq = record.info["CSQ"][0]
        gene = csq.split("|")[3]
        clinical = csq.split("|")[62]
        existing = csq.split("|")[17].split("&")

        #rs IDs use more than just the first!
        rsList = [rs for rs in existing if rs.startswith('rs')]
        if len(rsList) == 0:
            rs=''
        else:
            rs = ', '.join(rsList)
            # rs = rsList[0]

        # Total number of cosmic hemato hits on the position. Vep reports all cosmicId for that position.
        cosmicVepList = [cosmic for cosmic in existing if cosmic.startswith('CO')]
        if len(cosmicVepList) == 0:
            cosmicVep=''
        else:
            cosmicVep=', '.join(cosmicVepList)

        if len(cosmicVepList) == 0:
            cosmicN = ''
        else:
            cosmicN = 0
            for cosmicId in cosmicVepList:
                cmdCosmic = 'grep -w '+cosmicId+' /gluster-storage-volume/data/ref_data/COSMIC/COSMIC_v90_hemato_counts.txt | cut -f 16 '
                cosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                if len(cosmicNew) == 0:
                    cosmicNew = 0
                cosmicN += int(cosmicNew)

        transcript = csq.split("|")[10].split(":")[0]
        codingName = csq.split("|")[10].split(":")[1]
        consequence = csq.split("|")[1]


        #Population allel freq
        maxPopAf = record.info["CSQ"][0].split("|")[60]
        if len(maxPopAf) > 1:
            maxPopAf = round(float(maxPopAf),4)
        maxPop = record.info["CSQ"][0].split("|")[61]

        # Should pos be shifted if del??
        snv = [gene, record.contig, record.pos, record.ref, alt, af, record.info["DP"], transcript, codingName, consequence, cosmicVep, cosmicN, clinical, rs, maxPopAf, maxPop]

        #Artefact_file
        cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
        artLine = subprocess.run(cmdArt, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip() ##What happens if two hits?
        if artLine and record.ref == artLine.split()[2] and alt == artLine.split()[3]: #if pos exists and match in artefact file.
            orange.append(snv)
        else:
            # Germline /gluster-storage-volume/projects/wp2/nobackup/Twist_Myeloid/Artefact_files/Low_VAF_SNVs.txt
            cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
            germLine = subprocess.run(cmdGerm, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()

            if germLine and record.ref == germLine.split()[2] and alt == germLine.split()[3]: #if exists in germline file
                green.append(snv)
            else:
                white.append(snv)

        # row += 1

### Actually writing to the excelsheet
for line in white:
    if line[6] < 500:
        worksheetSNV.write_row(row,col,line, italicFormat)
    else:
        worksheetSNV.write_row(row,col,line)
    row +=1
for line in green:
    if line[6] < 500:
        worksheetSNV.write_row(row,col,line, green_italicFormat)
    else:
        worksheetSNV.write_row(row,col,line, greenFormat)
    row +=1
for line in orange:
    if line[6] < 500:
        worksheetSNV.write_row(row,col,line, orange_italicFormat)
    else:
        worksheetSNV.write_row(row,col,line, orangeFormat)
    row +=1


########################################

######### Indel sheet (3)##################
# Add genes as info before the actual table. Just use bed table as input? Sort uniq
worksheetIndel.set_column(1,2,10)
worksheetIndel.set_column(1,3,10)
worksheetIndel.write('A1', 'Pindel results', headingFormat)
worksheetIndel.write_row(1,0,emptyList,lineFormat)
with open(bedfile) as bed:
    genesDup = [line.split("\t")[3].strip() for line in bed]
    genes = set(genesDup)

genesString = ['Genes looked at: ']+list(genes)

for x in vcf_indel.header.records:
    if (x.key=='reference'):
        refI = x.value
worksheetIndel.write('A3', 'Sample: '+str(sample))
worksheetIndel.write('A4', 'Reference used: '+str(refI))
worksheetIndel.write('A5','Genes included: ')
row=5
for gene in genes:
    worksheetIndel.write('B'+str(row),gene)
    row+=1

row+=1
tableheading = ['Gene','Chr','Start','End','SV length','Af','Ref','Alt']
worksheetIndel.write_row('A'+str(row),tableheading, tableHeadFormat) #1 index
# row = 7 #0 index
col=0

for indel in vcf_indel.fetch():
    if indel.filter.keys()==["PASS"]:
        svlen = indel.info["SVLEN"]

        ads = indel.samples[sample]["AD"]
        af = int(ads[1])/(int(ads[0]) + int(ads[1]))
        # ad = ', '.join([str(i) for i in ads])

        if len(indel.alts) == 1:
            alt=indel.alts[0]
        else:
            print(indel.alts)
            sys.exit() ##Fix!!

        csqIndel = indel.info["CSQ"][0] #VEP annotation

        indelGene = csqIndel.split("|")[3]

        indelRow = [indelGene,indel.contig,indel.pos, indel.stop, svlen, af, indel.ref, alt]
        worksheetIndel.write_row(row,col,indelRow)
        row += 1

##########################################

######### Coverage sheet (4)#################
worksheetCov.set_column(1,3,10)
worksheetCov.set_column(1,4,10)
## Heading in sheet
worksheetCov.write('A1', 'CARTools coverage analysis', headingFormat)
worksheetCov.write_row('A2',emptyList,lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetCov.write('A4', description)
covHeadings = ['Region Name','Chr','Start','Stop','Mean Coverage','Length of Region']
worksheetCov.write_row('A6',covHeadings,tableHeadFormat) ## 1 index
row = 6 ## 0 index

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
# Number of low cov regions for the Overview sheet.
lowRegions = row - 6

worksheetOver.write(19,0,'Number of regions not covered by at least 100x: ')
worksheetOver.write(20,0, str(lowRegions))

##############################################

######### IVA sheet (6)#################
worksheetIVA.write('A1', 'Results from Variant Analysis ', headingFormat)
worksheetIVA.write_row('A2',emptyList,lineFormat)

worksheetIVA.write('A5', "Analysen utfördes i enlighet med dokumentationen.")
worksheetIVA.write('A6', "Eventuella avikelser:")
iva = ['DNA nr', 'Chromosome', 'Position', 'Gene Region', 'Gene Symbol', 'Transcript ID', 'Transcript Variant', 'Protein Variant', 'Variant Findings', 'Sample Genotype Quality', 'Read Depth', 'Allele Fraction', 'Translation Impact', 'dbSNP ID','1000 Genomes Frequency', 'ExAC Frequency', 'HGMD', 'COSMIC ID', 'Artefacts_without_ASXL1','ASXL1_variant_filter']
worksheetIVA.write_row(9,0, iva, tableHeadFormat)


#################################################

######### Prog Version sheet (7)#################
worksheetVersions.write('A1', 'Version Log', headingFormat)
worksheetVersions.write_row(1,0,emptyList,lineFormat)
worksheetVersions.write('A3','Sample: '+str(sample))
worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))
worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)


with open('containers.txt') as file:
    singularitys = [line.strip() for line in file]
    # singularitys.pop() ##Last slurm is always the makeContainersList rule.
    row = 7
    col = 0
    for singularity in singularitys:
        worksheetVersions.write_row(row,col,[singularity])
        row += 1



workbook.close()
