"""
File:
RNaseYProcessedOperonFinder.py

Identify operons that are specifically processed by RNase Y.
These operons should be present in the "processedOperons_wt1.csv" file
but absent in the "processedOperons_rny1.csv" file.

"""

import csv
import math
import numpy
from operon import Operon
from operonBank import OperonBank

operonBankWT1 = OperonBank("operonBankRD_wt1.txt")
operonBankRNY1 = OperonBank("operonBankRD_rny1.txt")

f1 = open("processedOperons_wt1.csv", 'rU')
processedGenes_wt1 = csv.reader(f1, delimiter=',')
processedGenesList_wt1 = []
for row in processedGenes_wt1:
     processedGenesList_wt1.append(row[0])
     
f2 = open("processedOperons_rny1.csv", 'rU')
processedGenes_rny1 = csv.reader(f2, delimiter=',')
processedGenesList_rny1 = []
for row in processedGenes_rny1:
     processedGenesList_rny1.append(row[0])

f3 = open("rnyProcessedOperons.csv", 'w')    
for i in xrange(len(operonBankWT1)):
    operonWT1 = operonBankWT1[i]
    operonRNY1 = operonBankRNY1[i]
    name = operonWT1.getName()
    if operonWT1.getOrientation() == "+":
        operonStart = operonWT1.getLeftBound()
        geneStart = operonWT1[0].getStart()
        operonEnd = operonWT1.getRightBound()
        geneEnd = operonWT1[-1].getEnd()
                
    else:
        operonStart = operonWT1.getRightBound()
        geneStart = operonWT1[-1].getEnd()
        operonEnd = operonWT1.getLeftBound()
        geneEnd = operonWT1[0].getStart()
        
    UTR5 = abs(operonStart - geneStart)
    UTR3 = abs(operonEnd - geneEnd)
         
    if (name in processedGenesList_wt1) and (name not in processedGenesList_rny1):
        string = str(operonWT1.getName() + "," + operonWT1[0].getName() + "," + \
                     operonWT1[-1].getName() + "," + str(UTR5) + "," + str(UTR3) + "\n")
        print string
        f3.write(string)      
f3.close()
##        print operonWT1.getName(), operonWT1[0].getName(), operonWT1[-1].getName(), str(UTR), \
##              round(operonWT1.getAverageRead(), 0), round(operonRNY1.getAverageRead(), 0)
