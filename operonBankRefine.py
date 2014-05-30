"""
File: operonBankRefine.py

Unify four operon banks.

1. Include all operons that are identified in the four datasets.
2. Re-define operons in the four operonBanks.

"""
import csv
import math
import cPickle
from sequence import Sequence
from geneBank import GeneBank
from operon   import Operon
from operonBank import OperonBank
from operonBankConstruct import operonJudge
from operonBankConstruct import boundFinder
from operonBankConstruct import operonFinder

"""Find all operon after-gaps that occur in the four operonBanks. """

geneBankFiles = ["geneBank_wt1.txt", "geneBank_wt2.txt", \
                 "geneBank_rny1.txt", "geneBank_rny2.txt"]
readsFiles = ["wt1.csv", "wt2.csv", "rny1.csv", "rny2.csv"]

operonBankFiles = [OperonBank("operonBank_wt1.txt"), \
                  OperonBank("operonBank_wt2.txt"), \
                  OperonBank("operonBank_rny1.txt"), \
                  OperonBank("operonBank_rny2.txt")]

gapList = []
n = 0
for i in xrange(4):
    operonBank = operonBankFiles[i]
    geneBank   = GeneBank(geneBankFiles[i])
    for operon in operonBank:
        gene = operon[-1]
        gap = geneBank.getIndex(gene.getName())
        if gap not in gapList:
            gapList.append(gap)
            n += 1
print "The number of after-gaps is: ", n
gapList.sort()
        

"""Group genes into operons based on the gap information """


for i in xrange(len(geneBankFiles)):
    geneBank = GeneBank(geneBankFiles[i])
    operonBankRD = OperonBank()
    with open(readsFiles[i], 'rU') as readsFile:
        reads = csv.reader(readsFile, delimiter='\t')
        readsList = []
        for row in reads:
            readsList.append(float(row[2]))
    
    
    n = 0
    operon = Operon("operon_" + str(n + 1))
    for gene in geneBank:        
        afterGap = gapList[n]
        if geneBank.getIndex(gene.getName()) < afterGap:
            operon.add(gene)
        elif geneBank.getIndex(gene.getName()) == afterGap:
            operon.add(gene)
            print gene.getName()
            n += 1
            operonBankRD.add(operon)
            operon.setLeftBound(operon[0].getStart())
            operon.setRightBound(operon[-1].getEnd())
            operon.setOrientation()
            operon.setReads(readsList)
            print operon.getName()
            print str(operon)
            operon = Operon("operon_" + str(n + 1))

    j = 0
    while j < len(operonBankRD) -1:
        operon1 = operonBankRD[j]
        operon2 = operonBankRD[j+1]        
        boundFinder(operon1, operon2, readsList)
        print operon1.getName(), (operon1.getLeftBound(), operon1.getRightBound()), \
              round(operon1.getAverageRead(), 2)
        
        
        j += 1
    operonBankRD[-1].setOrientation()
    
    
    fileName1 = 'operonBankRD_' + str(readsFiles[i])[0:-4] + '.txt'
    operonBankRD.save(fileName1)
    print fileName1 + " is generated."
    
    """Export operon boundaries to a new bed file """
    fileName2 = "operonBedRD_" + str(readsFiles[i])[0:-4] + ".bed"
    f2 = open(fileName2, 'w')
    for operon in operonBankRD:
        f2.write("Spy|"+ '\t' + str(operon.getLeftBound()) + '\t'
                 + str(operon.getRightBound()) + '\t'
                 + operon.getName() + '\t'
                 + '44' + '\t'
                 + operon.getOrientation() +'\n')
    f2.close()

    print fileName2 + " is generated."



