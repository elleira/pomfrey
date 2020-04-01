trusight#!/bin/python3.6
import sys
import csv
from pysam import VariantFile
import xlsxwriter ## need to download? Can I use same as pysam singularity?
from datetime import date
import subprocess

## Define sys.argvs
vcf_snv = VariantFile(sys.argv[1])
vcf_indel = VariantFile(sys.argv[2])
seqID = sys.argv[3]
cartool = sys.argv[4]
minCov = int(sys.argv[5])  ##grep thresholds ../somaticpipeline/qc/cartool/10855-17_Log.csv | cut -d ',' -f2
medCov = int(sys.argv[6])
maxCov = int(sys.argv[7])
bedfile = sys.argv[8]
hotspotFile = sys.argv[9]
artefactFile = sys.argv[10]
germlineFile = sys.argv[11]
output = sys.argv[12]

 ## Create execl file and sheets.
workbook = xlsxwriter.Workbook(output)
worksheetOver = workbook.add_worksheet('Overview')
worksheetTruSight = workbook.add_worksheet('TruSight')
worksheetSNV = workbook.add_worksheet('SNVs')
worksheetIndel = workbook.add_worksheet('InDel') #.... sys.argv[2]
worksheetLowCov = workbook.add_worksheet('Low Coverage') #... sys.argv[3]
worksheetHotspot = workbook.add_worksheet('Hotspot')
worksheetCov = workbook.add_worksheet('Coverage')
worksheetIVA = workbook.add_worksheet('IVA')
worksheetVersions = workbook.add_worksheet('Version')
## Define formats to be used.
headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
redFormat = workbook.add_format({'font_color': 'red'})

greenFormat = workbook.add_format({'bg_color': '#85e085'}) #, 'font_color': '#b3b3b3'
orangeFormat = workbook.add_format({'bg_color': '#ffd280'}) #font color gray 70%
green_italicFormat = workbook.add_format({'bg_color': '#85e085', 'italic': 'True'}) #'font_color': '#b3b3b3',
orange_italicFormat = workbook.add_format({'bg_color': '#ffd280', 'italic': 'True'}) #'font_color': '#b3b3b3',

## Define sample based on annotated vcf
sample = list(vcf_snv.header.samples)[0]
today=date.today()
emptyList=['','','','','','']

trusightGenes=['ABL1','ASXL1','ATRX','BCOR','BCORL1','BRAF','CALR','CBL','CBLB','CBLC','CDKN2A','CEBPA','CSF3R','CUX1','DNMT3A','ETV6','TEL','EZH2','FBXW7','FLT3','GATA1','GATA2','GNAS','HRAS','IDH1','IDH2','IKZF1','JAK2','JAK3','KDM6A','KIT','KRAS','KMT2A','MPL','MYD88','NOTCH1','NPM1','NRAS','PDGFRA','PHF6','PTEN','PTPN11','RAD21','RUNX1','SETBP1','SF3B1','SMC1A','SMC3','SRSF2','STAG2','TET2','TP53','U2AF1','WT1','ZRSR2']


######### Prog Version sheet (7)###########
worksheetVersions.write('A1', 'Version Log', headingFormat)
worksheetVersions.write_row(1,0,emptyList,lineFormat)
worksheetVersions.write('A3','Sample: '+str(sample))
# worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
# worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))
worksheetVersions.write('A7', 'Containers used: ', tableHeadFormat)


with open('containers.txt') as file:
    singularitys = [line.strip() for line in file]
    # singularitys.pop() ##Last slurm is always the makeContainersList rule.
    row = 7
    col = 0
    for singularity in singularitys:
        worksheetVersions.write_row(row,col,[singularity])
        row += 1
########################################

######### IVA sheet (6)#################
worksheetIVA.write('A1', 'Results from Variant Analysis ', headingFormat)
worksheetIVA.write_row('A2',emptyList,lineFormat)

worksheetIVA.write('A5', "Analysen utfördes i enlighet med dokumentationen.")
worksheetIVA.write('A6', "Eventuella avikelser:")
iva = ['DNA nr', 'Chromosome', 'Position', 'Gene Region', 'Gene Symbol', 'Transcript ID', 'Transcript Variant', 'Protein Variant', 'Variant Findings', 'Sample Genotype Quality', 'Read Depth', 'Allele Fraction', 'Translation Impact', 'dbSNP ID','1000 Genomes Frequency', 'ExAC Frequency', 'HGMD', 'COSMIC ID', 'Artefacts_without_ASXL1','ASXL1_variant_filter']
worksheetIVA.write_row(9,0, iva, tableHeadFormat)


#########################################

########### Coverage ####################
##Lägga till hela MeanFullCoverage... som filtrerbar lista
#Number of lines in MeanFullCoverage
covFullFile = cartool.replace("_MeanCoverageShortList.csv", "_MeanCoverageFullList.csv")
def file_lengthy(fname):
        with open(fname) as f:
                for i, l in enumerate(f):
                        pass
        return i + 1

numRows = file_lengthy(covFullFile)
worksheetCov.write('A1','Average Coverage',headingFormat)
worksheetCov.write_row('A2',emptyList,lineFormat)
worksheetCov.write('A3', 'Sample: '+str(sample))
worksheetCov.write('A4', 'Averge coverage of each region in bedfile')

## Fixa data i table först

tableLines=[]
with open(covFullFile) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    ## Skip header
    for line in readCSV:
        while not line[-1]: ## Remove empty fields at the end of line.
            line.pop()
        for i in range(1, int((len(line)-1)/5+1)): ##For each region/ with lower cov, create a new row.
            start = 1+5*(i-1)
            end = 1+5*i
            covRow = [line[0]]+line[start:end]
            tableLines.append(covRow)
            # worksheetLowCov.write_row(row,col,[line[0]]+line[start:end])
# import pdb; pdb.set_trace()
tableArea='A6:F'+str(len(tableLines)+6) ##rows of full list
headerListDict = [{'header':'Region Name'}, {'header':'Chr'}, {'header':'Start'},{'header':'End'}, {'header':'Mean Coverage'}, {'header':'Length of Region'},]
worksheetCov.add_table(tableArea, {'data':tableLines,'columns':headerListDict,'style': 'Table Style Light 1'})
# worksheetCov.write_row('A6', )#, tableHeadFormat) #table header How sokbar??


#########################################

########## Hotspot sheet (5)#############
worksheetHotspot.write('A1', 'Hotspot Coverage', headingFormat)  #headerFormat)
worksheetHotspot.write('A3', 'Sample: '+str(sample))
worksheetHotspot.set_column(1,2,10)
worksheetHotspot.write_row('A5',['Chr', 'Pos', 'Depth', 'Gene'], tableHeadFormat)

lowPos = 0
row=5
hotspotTable=[]
with open(cartool.replace("_MeanCoverageShortList.csv", "_coverageShortHotspot.tsv"),'r') as hotFile: #Always the same as bedfile just without region
    for dpLine in hotFile:
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
                if chrCov == chrHS and lowHS == posCov: #Always the same as bedfile just without region, How does CARTool handle if hotspot is longer than 1bp?
                    hotspotTable.append([chrCov,posCov, dp, hotspot[3].rstrip()])

hotspotTable.sort(key=lambda x: float(x[2]))
for hotLine in hotspotTable:
    if int(hotLine[2]) <= medCov: ##How to get the number from configfile?!
        worksheetHotspot.write_row(row,0,hotLine,redFormat)
        row += 1
        lowPos += 1
    else:
        worksheetHotspot.write_row(row,0,hotLine)
        row += 1

############################################

######### Low Coverage sheet (4)################
worksheetLowCov.set_column(1,3,10)
worksheetLowCov.set_column(1,4,10)
## Heading in sheet
worksheetLowCov.write('A1', 'CARTools coverage analysis', headingFormat)
worksheetLowCov.write_row('A2',emptyList,lineFormat)
worksheetLowCov.write('A3', 'Sample: '+str(sample))
description = 'Gene Regions with coverage lower than '+str(minCov)+'x.'
worksheetLowCov.write('A4', description)
covHeadings = ['Region Name','Chr','Start','Stop','Mean Coverage','Length of Region']
worksheetLowCov.write_row('A6',covHeadings,tableHeadFormat) ## 1 index
row = 6 ## 0 index

covLines=[]
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
            covRow = [line[0]]+line[start:end]
            covLines.append(covRow)
            # worksheetLowCov.write_row(row,col,[line[0]]+line[start:end])
            # row += 1
    #sort based on Coverage
    covLines.sort(key=lambda x: x[4])
    for line in covLines:
        worksheetLowCov.write_row(row,col,line)
        row += 1
        # covLines
# Number of low cov regions for the Overview sheet.
lowRegions = row - 6

###########################################

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
tableheading = ['SeqID','DNAnr','Gene','Chr','Start','End','SV length','Af','Ref','Alt','Dp','Canonical','Mutation cds','ENSP','Max popAF','Max Pop']
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

        maxPopAfIndel = csqIndel.split("|")[56] #[57]
        if len(maxPopAfIndel) > 1:
            maxPopAfIndel = round(float(maxPopAfIndel),4)
        maxPopIndel = csqIndel.split("|")[57] #[61]

        indelTranscript = csqIndel.split("|")[10].split(":")[0]
        indelCodingName = csqIndel.split("|")[10].split(":")[1]
        indelEnsp = csqIndel.split("|")[11]
        # indelIGV =

        indelRow = [seqID,sample,indelGene,indel.contig,indel.pos, indel.stop, svlen, af, indel.ref, alt, indel.info["DP"], indelTranscript, indelCodingName,indelEnsp ,maxPopAfIndel, maxPopIndel]
        worksheetIndel.write_row(row,col,indelRow)
        row += 1

##########################################

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
tableheading = ['SeqID','DNAnr','Gene','Chr','Pos','Ref','Alt', 'AF', 'DP', 'Canonical Transcript','Mutation cds', 'ENSP' ,'Consequence','COSMIC ids on position','N COSMIC Hemato hits on position','Clinical significance', 'dbSNP','Max popAF','Max Pop','IGV image']
worksheetSNV.write_row('A18', tableheading, tableHeadFormat) #1 index
row = 18 #0 index
col=0
white=[]
green=[]
orange=[]
whiteIGV=[]
underFive=[] #put after green and orange but still white
underFiveIGV=[] #put after green and orange but still white
trusightSNV=[] #trusight genes only
trusightSNVigv=[] #trusight genes only

for record in vcf_snv.fetch():
    synoCosmicN = 0
    if   record.filter.keys()==["Syno"]: #Only if Syno not and popAF.   any(x in "Syno" for x in record.filter.keys()):
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
        if af >= 0.03:
            csq = record.info["CSQ"][0]
            gene = csq.split("|")[3]
            clinical = csq.split("|")[58] #[59]
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
                    cmdCosmic = 'grep -w '+cosmicId+' /apps/bio/singularities/gms_hematology/bedfiles/COSMIC_v90_hemato_counts.txt | cut -f 16 '
                    cosmicNew = subprocess.run(cmdCosmic, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                    if len(cosmicNew) == 0:
                        cosmicNew = 0
                    cosmicN += int(cosmicNew)

            transcript = csq.split("|")[10].split(":")[0]
            codingName = csq.split("|")[10].split(":")[1]
            consequence = csq.split("|")[1]
            ensp = csq.split("|")[11]

            #Population allel freq
            # maxPopAf = record.info["CSQ"][0].split("|")[57] #[60]
            # if len(maxPopAf) > 1:
            #     maxPopAf = round(float(maxPopAf),4)
            # maxPop = record.info["CSQ"][0].split("|")[58] #[61]

            popFreqsPop=['AF', 'AFR_AF','AMR_AF','EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF']
            popFreqAllRaw=record.info["CSQ"][0].split("|")[41:56] #[42:57]
            if any(popFreqAllRaw) and max([float(x) if x else 0 for x  in popFreqAllRaw]) != 0: #if all not empty
                popFreqAll=[float(x) if x else 0 for x  in popFreqAllRaw]
                maxPopAf=max(popFreqAll)
                maxIndex=[i for i, j in enumerate(popFreqAll) if j == maxPopAf]
                if len(maxIndex) == 1:
                    maxPop=popFreqsPop[maxIndex[0]]
                else:
                    popFreqPops=[popFreqsPop[x] for x in maxIndex]
                    maxPop="&".join(popFreqPops)
            else:
                maxPopAf=''
                maxPop=''
            ##IGV image path for each SNV
            igv="external:IGV/"+gene+"-"+record.contig+"_"+str(int(record.pos)-1)+"_"+str(int(record.pos)-1+len(alt))+".svg"

            snv = [seqID,sample,gene, record.contig, record.pos, record.ref, alt, af, record.info["DP"], transcript, codingName, ensp, consequence, cosmicVep, cosmicN, clinical, rs, maxPopAf, maxPop]
            #Append line with sample and rundate to rolling list of artefacts..
            with open("/apps/bio/singularities/gms_hematology/bedfiles/twistVariants.txt", "a") as appendfile:
                variants = snv+["\n"]
                appendfile.write('\t'.join(str(e) for e in variants))

            #Artefact_file
            cmdArt = 'grep -w '+str(record.pos)+' '+artefactFile
            artLines = subprocess.run(cmdArt, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip() ##What happens if two hits?!!!
            artefact_variant = 0
            # import pdb; pdb.set_trace()
            for artLine in artLines.split("\n"):
                if artLine and record.ref == artLine.split()[2] and alt == artLine.split()[3]: #if pos exists and match in artefact file.
                    orange.append(snv)
                    artefact_variant = 1
                    break
            if artefact_variant == 0:
                # Germline /gluster-storage-volume/projects/wp2/nobackup/Twist_Myeloid/Artefact_files/Low_VAF_SNVs.txt
                cmdGerm = 'grep -w '+str(record.pos)+' '+germlineFile
                germLines = subprocess.run(cmdGerm, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
                germline_variant = 0
                for germLine in germLines.split("\n"):
                    if germLine and record.ref == germLine.split()[2] and alt == germLine.split()[3]: #if exists in germline file
                        green.append(snv)
                        germline_variant = 1
                        break
                if germline_variant == 0:
                    if float(af) < 0.05:
                        underFive.append(snv)
                        underFiveIGV.append(igv)
                    else:
                        white.append(snv)
                        whiteIGV.append(igv)
                    if gene in trusightGenes:
                        trusightSNV.append(snv)
                        trusightSNVigv.append(igv)
### Actually writing to the excelsheet
i=0
for line in white:
    if line[8] < 500:
        worksheetSNV.write_row(row,col,line, italicFormat)
        worksheetSNV.write_url('S'+str(row+1), whiteIGV[i],string = "IGV image")
    else:
        worksheetSNV.write_row(row,col,line)
        # import pdb; pdb.set_trace()
        worksheetSNV.write_url('S'+str(row+1), whiteIGV[i],string = "IGV image")
    row +=1
    i+=1

for line in green:
    if line[8] < 500:
        worksheetSNV.write_row(row,col,line, green_italicFormat)
    else:
        worksheetSNV.write_row(row,col,line, greenFormat)
    row +=1

for line in orange:
    if line[8] < 500:
        worksheetSNV.write_row(row,col,line, orange_italicFormat)
    else:
        worksheetSNV.write_row(row,col,line, orangeFormat)
    row +=1

i=0
for line in underFive:
    if line[8] < 500:
        worksheetSNV.write_row(row,col,line, italicFormat)
        worksheetSNV.write_url('S'+str(row+1), underFiveIGV[i],string = "IGV image")
    else:
        worksheetSNV.write_row(row,col,line)
        worksheetSNV.write_url('S'+str(row+1), underFiveIGV[i],string = "IGV image")
    row +=1
    i+=1

#############################################

######### TruSight varianter ###############
# Variants or snv rows from SNV sheet.
worksheetTruSight.write('A1', 'TruSight variants found', headingFormat)
worksheetTruSight.write('A3', 'Sample: '+str(sample))
worksheetTruSight.write('A4', 'Reference used: '+str(refV))
# worksheetTruSight.write('A4', 'Callers used: VardictJava v.1,6 ? Dubbelkolla!, Pisces 5.2.11, Freebayes v1.1.0, LoFreq v.2.1.3.1, TruSighter v.2.1.3.1')
worksheetTruSight.write('A6', 'VEP: '+vepline ) #, textwrapFormat)
worksheetTruSight.write('A8', 'The following filters were applied: ')
worksheetTruSight.write('B9','Coverage >= 100x')
worksheetTruSight.write('B10','Population freq (KGP, gnomAD, NHLBI_ESP ) <= 2%')
worksheetTruSight.write('B11','Biotype is protein coding')
worksheetTruSight.write('B12','Consequence not deemed relevant')


worksheetTruSight.write('A14','Only variants in genes from TruSight panel:')
worksheetTruSight.write_row(14,0,trusightGenes)
worksheetTruSight.write('A16','For all variants see: ')
worksheetTruSight.write_url('B16', "internal:'SNVs'!A1",string='SNVs')

worksheetTruSight.write('A18','Coverage below 500x', italicFormat)
worksheetTruSight.write_row('A20', tableheading, tableHeadFormat) #1 index
row=20 # 0 index
i=0
for line in trusightSNV:
    if line[8] < 500:
        worksheetTruSight.write_row(row,col,line, italicFormat)
        worksheetTruSight.write_url('S'+str(row+1), trusightSNVigv[i],string = "IGV image")
    else:
        worksheetTruSight.write_row(row,col,line)
        # import pdb; pdb.set_trace()
        worksheetTruSight.write_url('S'+str(row+1), trusightSNVigv[i],string = "IGV image")
    row +=1
    i+=1




#############################################

######### Overview sheet (1) ################

worksheetOver.write(0,0, sample, headingFormat)
worksheetOver.write(1,0, "SeqID: "+seqID)
worksheetOver.write(2,0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(3,0, emptyList,lineFormat)

worksheetOver.write(4,0, "Created by: ")
worksheetOver.write(4,4, "Valid from: ")
worksheetOver.write(5,0, "Signed by: ")
worksheetOver.write(5,4, "Document nr: ")
worksheetOver.write_row(6,0,emptyList,lineFormat)

worksheetOver.write(7,0,"Sheets:", tableHeadFormat)
worksheetOver.write_url(8,0,"internal:'TruSight'!A1", string='TruSight Variants')
worksheetOver.write_url(9,0,"internal:'SNVs'!A1", string='Variants analysis')
worksheetOver.write_url(10,0,"internal:'Indel'!A1", string = 'Indel variants')
worksheetOver.write_url(11,0,"internal:'Coverage'!A1", string = 'Positions with coverage lower than 100x')
worksheetOver.write_url(12,0,"internal:'Hotspot'!A1", string = 'Coverage of hotspot positions')
worksheetOver.write_url(13,0,"internal:'Version'!A1", string = 'Version Log')
worksheetOver.write_row(15,0,emptyList,lineFormat)


##Add avg. cov and clonalisy
cartoolLog=cartool.replace("_MeanCoverageShortList.csv", "_Log.csv")
cmdAvgCov = 'grep Depth '+cartoolLog+' | cut -d"," -f2 | cut -f1 -d" "'
avgCov = subprocess.run(cmdAvgCov, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()

duplicationFile = cartool.replace("_MeanCoverageShortList.csv", "_DuplicationMetrics.txt")
cmdDupl = 'grep -A1 PERCENT '+duplicationFile+' | tail -1 | cut -f9'
duplicateLevel = subprocess.run(cmdDupl, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()

breadthCmd = 'grep "Mean Coverage Breadth:" '+cartoolLog + ' | cut -f2- -d"," '
breadth = subprocess.run(breadthCmd, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
#
# worksheetOver.write(17,0, "Avg. coverage [x]: ")
# worksheetOver.write(17,2, avgCov)
# worksheetOver.write(18,0, "Duplicationlevel [%]: ")
# worksheetOver.write(18,2, str(round(float(duplicateLevel)*100,2)))

worksheetOver.write_row(17,0,['SeqId', 'DNAnr', 'Avg. coverage [x]','Duplicationlevel [%]',str(minCov)+'x',str(medCov)+'x',str(maxCov)+'x'],tableHeadFormat)
worksheetOver.write_row(18,0,[seqID, sample, avgCov, str(round(float(duplicateLevel)*100,2))]+breadth.split(','))

if lowPos == 0: #From Hotspot sheet
    worksheetOver.write(21,0,'Number of positions from the hotspot list not covered by at least 500x: ')
    worksheetOver.write(22,0, str(lowPos))
else:
    worksheetOver.write(21,0,'Number of positions from the hotspot list not covered by at least 500x: ')
    worksheetOver.write(22,0, str(lowPos), redFormat)
    worksheetOver.write_url(23,0,"internal:'Hotspot'!A1" ,string = 'For more detailed list see hotspotsheet ')


worksheetOver.write(24,0,'Number of regions not covered by at least 100x: ') #From Cov sheet
worksheetOver.write(25,0, str(lowRegions)) #From Cov sheet
worksheetOver.write(28,0,'Hotspotlist: '+hotspotFile)
worksheetOver.write(29,0,'Artefact file: '+artefactFile)
worksheetOver.write(30,0,'Germline file: '+germlineFile)

######################################

####### Prog Version sheet (7), added last #######
worksheetVersions.write('A5', 'Variant calling reference used: '+str(refV))
worksheetVersions.write('A6', 'Pindel reference used: '+str(refI))

######################################
workbook.close()
