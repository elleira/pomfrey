#!/bin/python3
import sys

SampleSheetUsed = sys.argv[1]
outFile = sys.argv[2]

samplesLines = []
startReading = 0
with open(SampleSheetUsed, 'r') as file:
    lines = [line.strip() for line in file]
    for line in lines:
        if line.startswith("Date"): #Containers
            date = line.split(',')[1]
        if startReading == 1: ##Once reached [Data]
            samplesLines.append(line.split(','))
        if line == "[Data]":
            startReading = 1
# samples.pop() #Remove any empty are there empty line at end?!
samplesLines = samplesLines[1:] #Remove header from SampleSheetUsed
# sampleSheetSamples = [string for string in samples if string != ""]
#How to remove trailing lines in samplesLines??

with open(outFile, 'w') as outfile:
    for line in samplesLines:
        experimentID = line[7]+line[2].split('_')[3]
        sample = line[1]
        #For each TM sample that is starts with D, FD or number only. Only add samples that match the clinics namning
        if line[7] == "TM" and (sample.startswith('D') or sample.startswith('FD') or sample[0].isdigit()):
            outfile.write('{ \"experiment.wp\": \"WP2\", \"experiment.prep\": \"Fresh\", \"@timestamp\": \"'+date+'T01:01:01.000Z\", ')
            outfile.write('\"experiment.method\": \"Twist_Myeloid\", \"experiment.rerun\": false, \"experiment.user\": \"Unknown\", ')
            outfile.write('\"experiment.tissue\": \"Hematology\", \"experiment.id\": \"'+experimentID+'\"experiment.sample\": \"'+sample+'\"}\n')

#
# {"experiment.wp": "WP2", "experiment.prep": "Fresh", "@timestamp": "2020-04-06T01:01:01.000Z",
# "experiment.method": "TruSight_Myeloid", "experiment.rerun": false, "experiment.user": "JHS",
# "experiment.tissue": "Hematology", "experiment.id": "TruSight_405_2", "experiment.sample": "D20-02190"}
