"""
File: operonBankConstruct2.py

Construct operonBank objects based on geneBank and readsList
Construct bed file based on operonBank.

"""
import csv
import math
import cPickle
from sequence import Sequence
from geneBank import GeneBank
from operon   import Operon
from operonBank import OperonBank


geneBankFiles = ["geneBank_wt1.txt", "geneBank_wt2.txt", \
                 "geneBank_rny1.txt", "geneBank_rny2.txt"]
readsFiles = ["wt1.csv", "wt2.csv", "rny1.csv", "rny2.csv"]

def operonJudge(seq1, seq2, readsList):
    """Define a function to determine whether two sequences belong to different operons.
    Criteria: 1. The average reads of the two sequences are more than 4-fold different.
              2. There is a 'dent' between the two sequences, the bottom of which
                 is lower than half of the average read of the less expressed gene.
              3. There two sequences have different orientations.
              4. The two sequences are over 100 bp apart.
    We consider that the two sequences lie in different operons if any of the criteria
    is satisfied.
    Precondition: seq1 and seq2 must be consecutive genes on the chromosome.
    Returns True if two sequences belong to two different operons, and False otherwise."""
        
    def expressJudge(seq1, seq2):
        """Returns True if seq1 and seq2 are differentially expressed,
        and False otherwise."""
        aveRead1 = seq1.getAverageRead()
        aveRead2 = seq2.getAverageRead()
        return aveRead1 >= aveRead2 * 4 or aveRead1 <= aveRead2 * 0.25             

    def dentJudge(seq1, seq2, readsList):
        """Returns True if there exists a dent between the two sequences,
        and False otherwise."""    
        IGRStart = seq1.getEnd() + 1
        IGREnd   = seq2.getStart() - 1
        if IGRStart >= IGREnd:
            return False
        else:
            IGR = Sequence(readsList, (IGRStart, IGREnd))
            minIGRRead = IGR.getMinRead()
##            print "minIGRRead = ", minIGRRead
            
            aveRead1 = seq1.getAverageRead()
            aveRead2 = seq2.getAverageRead()
            minAveRead = min(aveRead1, aveRead2)
##            print "minAveRead = ", minAveRead

            return minIGRRead <= minAveRead * 0.5
            
    def strandJudge(seq1, seq2):
        """Returns True if the two sequences are on the same strand,
        and False otherwise."""
        strandness1 = seq1.getOrientation()
        strandness2 = seq2.getOrientation()
        return strandness1 != strandness2
             
    def distanceJudge(seq1, seq2):
        """Returns True if the two sequences are over 100bp apart,
        and False otherwise."""
        distance = seq2.getStart() - seq1.getEnd()
        return distance > 100

    return expressJudge(seq1, seq2) or \
           dentJudge(seq1, seq2, readsList) or \
           strandJudge(seq1, seq2) or \
           distanceJudge(seq1, seq2)
        
def operonFinder(geneBank, readsList):
    """Group genes into operons."""
         
    operon1 = Operon("operon_1")
    gene1 = geneBank[0]
    operon1.add(gene1)
    operonBank = OperonBank()
    operonBank.add(operon1)
 
    j = 0    
    for i in xrange(1, len(geneBank)):
        currentOperon = operonBank[-1]
        currentGene  = geneBank[i]
        previousGene = geneBank[i-1]
        if operonJudge(previousGene, currentGene, readsList) == False:
            currentOperon.add(currentGene)
        else:
            j += 1
            name = "operon_" + str(j + 1)
            operon = Operon(name)
            operon.add(currentGene)
            operonBank.add(operon)

    for operon in operonBank:
        operon.setOrientation()
        operon.setLeftBound(operon[0].getStart())
        operon.setRightBound(operon[-1].getEnd())
        operon.setReads(readsList)
        
    return operonBank

def boundFinder(operon1, operon2, readsList):
    """Defines boundaries between two consecutive operons
    using sliding window method to find the region with maximum CV."""
    WINDOWWIDTH = 25

    upGene     = operon1[-1]
    downGene   = operon2[0]
    
    leftEdge1  = upGene.getEnd() - 200
    leftEdge2  = (upGene.getStart() + upGene.getEnd() )/ 2
    leftEdge   = max(leftEdge1, leftEdge2)
    rightEdge1 = downGene.getStart() + 200
    rightEdge2 = (downGene.getStart() + downGene.getEnd()) / 2
    rightEdge  = min(rightEdge1, rightEdge2)
    midPoint  = (upGene.getEnd() + downGene.getStart()) / 2
    workSpace = Sequence(readsList, (leftEdge, rightEdge))
    breakPoint = workSpace.getMinReadLocation()

    if len(workSpace) <= WINDOWWIDTH:        
        operon1.setRightBound(breakPoint)
        operon2.setLeftBound(breakPoint)
    else:
        refCV = min(upGene.getLogCV(), downGene.getLogCV())
        space1 = Sequence(readsList, (leftEdge, breakPoint))
        space2 = Sequence(readsList, (breakPoint, rightEdge))
        # Determine the turning point in the left space:     
        if len(space1) <= WINDOWWIDTH:
            maxWindow1 = space1
        else:
            start1 = space1.getStart()
            maxWindow1 = Sequence(readsList, (start1, start1 + WINDOWWIDTH))      
            while start1 + WINDOWWIDTH < space1.getEnd():
                window1 = Sequence(readsList, (start1, start1 + WINDOWWIDTH))
                if maxWindow1.getLogCV() < window1.getLogCV() and \
                   window1.getReads()[0] > window1.getReads()[-1]:
                    maxWindow1 = window1
                start1 += 1
        if maxWindow1.getLogCV() >= refCV * 2:
            turnPoint1 = maxWindow1.getMinReadLocation()
        else:
            turnPoint1 = -1
            
        # Determine the turning point in the right space:     
        if len(space2) <= WINDOWWIDTH:
            maxWindow2 = space2
        else:
            start2 = space2.getStart()
            maxWindow2 = Sequence(readsList, (start2, start2 + WINDOWWIDTH))
            while start2 + WINDOWWIDTH < space2.getEnd():
                window2 = Sequence(readsList, (start2, start2 + WINDOWWIDTH))
                if maxWindow2.getLogCV() < window2.getLogCV() and \
                   window2.getReads()[0] < window2.getReads()[-1]:
                    maxWindow2 = window2
                start2 += 1
        if maxWindow2.getLogCV() >= refCV * 2:
            turnPoint2 = maxWindow2.getRightMinReadLocation()
        else:
            turnPoint2 = -1
            
            
        # Assign turning points to each space.
        if turnPoint1 > 0 and turnPoint2 > 0:
            operon1.setRightBound(turnPoint1)
            operon1.setRightBoundPrecision(True)
            operon2.setLeftBound(turnPoint2)
            operon2.setLeftBoundPrecision(True)
            
        elif turnPoint1 > 0 and turnPoint2 < 0:
            operon1.setRightBound(turnPoint1)
            operon1.setRightBoundPrecision(True)
            operon2.setLeftBound(turnPoint1+1)
            operon2.setLeftBoundPrecision(False)
            
        elif turnPoint1 < 0 and turnPoint2 > 0:
            operon1.setRightBound(turnPoint2-1)
            operon1.setRightBoundPrecision(False)
            operon2.setLeftBound(turnPoint2)
            operon2.setLeftBoundPrecision(True)
            
        else:
            operon1.setRightBound(midPoint)
            operon1.setRightBoundPrecision(False)
            operon2.setLeftBound(midPoint+1)
            operon2.setLeftBoundPrecision(False)


def operonBankConstruct(geneBank, readsList, strainName):
    """Construct a operonBank object. Define the boundaries of each operon."""

    operonBank = operonFinder(geneBank, readsList)

    j = 0
    while j < len(operonBank) -1:
        operon1 = operonBank[j]
        operon2 = operonBank[j+1]
        operon1.setOrientation()
        boundFinder(operon1, operon2, readsList)
        print operon1.getName(), (operon1.getLeftBound(), operon1.getRightBound()), \
              operon1.getLeftBoundPrecision(), operon1.getRightBoundPrecision()
        
        j += 1

    bankName = "operonBank_" + strainName + ".txt"
    operonBank.save(bankName)
    print bankName + " is generated."
    
    """Export operon boundaries to a new bed file """
    bedName = "operonBank_" + strainName + ".bed"
    f2 = open(bedName, 'w')
    for operon in operonBank:
        f2.write("Spy|"+ '\t' + str(operon.getLeftBound()) + '\t'
                 + str(operon.getRightBound()) + '\t'
                 + operon.getName() + '\t'
                 + '44' + '\t'
                 + operon.getOrientation() +'\n')
    f2.close()

    print bedName + " is generated."


def main():
    for i in xrange(len(geneBankFiles)):
        geneBank = GeneBank(geneBankFiles[i])
        
        with open(readsFiles[i], 'rU') as readsFile:
            reads = csv.reader(readsFile, delimiter='\t')
            readsList = []
            for row in reads:
                readsList.append(float(row[2]))

        strainName = readsFiles[i][0:-4]
        operonBankConstruct(geneBank, readsList, strainName)

main()
