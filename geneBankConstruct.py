"""
File: geneBankConstruct.py
Generates four GeneBank objects representing genes in wt1, wt2, rny1, and rny2.
Save each object using cPickle."""

import csv
import math
import cPickle
from sequence import Sequence
from geneBank import GeneBank

fileList = ["wt1.csv", "wt2.csv", "rny1.csv", "rny2.csv"]
for fileName in fileList:
    geneFile = open('Spy49_allGenes.csv', 'rU')
    genes = csv.reader(geneFile)

    with open(fileName, 'rU') as readsFile:
        reads = csv.reader(readsFile, delimiter='\t')
        readsList = []
        for row in reads:
            readsList.append(float(row[2]))
            
    geneBank = GeneBank()
    for row in genes:
        name = row[0]
        start = int(row[1])
        end = int(row[2])
        orientation = row[3]
        seq = Sequence(readsList, (start, end), name, orientation) 
        geneBank.add(seq)
    geneBankName = "geneBank_" + fileName[0:-4] + ".txt"
    geneBank.save(geneBankName)
    geneFile.close()
