#!/bin/python3.6

import sys
import subprocess

rawDataFolder = sys.argv[1]
defaultConfig = sys.argv[2]
sequencerun = sys.argv[3]

samplesLines = []
startReading = 0
with open(rawDataFolder+'/SampleSheetUsed.csv', 'r') as file:
    lines = [line.strip() for line in file]
    for line in lines:
        if startReading == 1: ##Once reached [Data]
            samplesLines.append(line.split(','))
        if line == "[Data]":
            startReading = 1
    # samples.pop() #Remove any empty are there empty line at end?!
samplesLines = samplesLines[1:]

subprocess.run('cp '+defaultConfig+' '+sequencerun+'_config.yaml',stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
with open(sequencerun+'_config.yaml', 'a') as configfile:
    for line in samplesLines:
        # import pdb; pdb.set_trace()
        if line[7] == "TM": #and not line[1].startswith('HD829'): 
            r1cmd = 'ls '+rawDataFolder+'/'+line[1]+'*R1_001.fastq.gz'
            r1 = subprocess.run(r1cmd, stdout=subprocess.PIPE,shell = 'TRUE').stdout.decode('utf-8').strip()
            r1 = r1.replace('//','/')
            configfile.write('  \"'+line[1]+'": \"'+r1+'\"\n')
    configfile.write('seqID: \n')
    configfile.write('  sequencerun: \"'+sequencerun+'\"\n')
