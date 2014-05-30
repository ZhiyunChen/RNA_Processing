"""
File: operon.py

Define an Operon class to describe aa operon.

"""

import math
import numpy
from sequence import Sequence

with open("SpyGenome.txt", 'r') as myfile:
    SpyGenome = myfile.read().replace('\n', '')

class Operon(object):
    """Defines an Operon class representing a cluster of genes controlled
    by the same promoter."""

    SEQUENCESTR = SpyGenome

    def __init__(self, name):
        """Constructs an Operon object """
        self._name         = name
        self._sequences    = []
        self._size         = 0
        self._orientation  = ''
        self._leftBound    = 0
        self._leftBoundPrecision = False
        self._rightBound   = -1
        self._rightBoundPrecision = False
        self._reads        = ()
        self._baseSequence = ''
        self._probeHalfLife = []
        

    def __str__(self):
        """Returns a string representation of gene names in an operon."""
        operonStr = ''
        operonStr += self._name + ": "
        for seq in self._sequences:
            operonStr += seq.getName() + ", "
        operonStr += self._orientation + ","
        operonStr += '(' + str(self._leftBound) + ", "
        operonStr += str(self._rightBound) + ')'
        return operonStr

    def __len__(self):
        """Returns the number of sequences in an operon."""
        return self._size

    def __getitem__(self, index):
        """Returns the sequence at indexed position. """
        return self._sequences[index]

    def getGenes(self):
        """Returns all genes in the operon. """
        return self._sequences

    def add(self, seq):
        """Append seq to the end of the sequence list. Returns None """
##        if self._orientation != '' and \
##           seq.getOrientation() != self._orientation:
##            raise ValueError, "Inconsistent sequence orientation."
        self._sequences.append(seq)
        self._size += 1    

    def setOrientation(self):
        """Preconditions: 1. The operon is not empty.
        2. All sequences in the operon have the same orientation.
        Defines the orientation of an operon on the chromosme. """
        if self._size == 0:
            raise ValueError, "The operon is empty."""
##        for i in xrange(self._size - 1):
##            seq1 = self._sequences[i]
##            seq2 = self._sequences[i+1]
##            if seq1.getOrientation() != seq2.getOrientation():
##                raise ValueError, "Inconsistent sequence orientation"
        self._orientation = self._sequences[0].getOrientation()
               
    def setLeftBound(self, leftBound):
        """Defines the left boundary of an operon. """
        self._leftBound = leftBound

    def setRightBound(self, rightBound):
        """Defines the right boundary of an operon. """
        self._rightBound = rightBound

    def setReads(self, readsList):
        """Generates a tuple showing RNA-seq reads of all bases in an operon.
        Returns None. """
        self._reads = tuple(readsList[self._leftBound:self._rightBound])

    def getReads(self):
        """Returns the RNA-seq reads in an operon. """
        return self._reads

    def _logTransform(self):
        """Transform reads to log-reads"""
        logReadsList = []
        for read in self._reads:
            if read == 0:
                read = 1
            logRead = math.log(read, 2)
            logReadsList.append(logRead)
        return logReadsList          
   
    def getAverageRead(self):
        """Returns the average read of the operon. The reads are log-transformed \
        before averaging and expo-transformed after averaging."""
        if self._reads == ():
            return 0
        else:
            averageLogRead = numpy.average(self._logTransform())
            averageRead = 2 ** averageLogRead
            return averageRead
        
    def setBaseSequence(self):
        """Defines the base sequence of an operon."""
        self._baseSequence = Operon.SEQUENCESTR[self._leftBound:\
                                                self._rightBound]

    def setProbeHalfLife(self, probeHalfLifeList):
        """Set probe half-lives to the operon """
        for i in xrange(len(probeHalfLifeList)):
            position = probeHalfLifeList[i][0]
            orientation = probeHalfLifeList[i][1]
            if position >= self._leftBound and \
               position <= self._rightBound and \
               orientation == self._orientation:
                self._probeHalfLife.append(probeHalfLifeList[i])

    def getProbeHalfLife(self):
        return self._probeHalfLife
    
    def _HLlogTransform(self):
        """Transform probel-level half-lives to log-half-lives.
        Delete the probe if the half-life is -1 (which indicates a poor probe). """
        HLLogList = []
        for HLTuple in self._probeHalfLife:
            name = HLTuple[0]
            orientation = HLTuple[1]
            halfLife = HLTuple[2]
            if halfLife > 0:
                HLLog = math.log(halfLife, 2)
                HLLogList.append(HLLog)
        return HLLogList      

    def getHLLogCV(self):
        """Returns the coefficient variation of probe-level
        half-lives in an operon. """
        if self._probeHalfLife == []:
            return 0
        else:
            averageHL = numpy.average(self._HLlogTransform())
            stdHL     = numpy.std(self._HLlogTransform())
            if average != 0:
                CV = stdHL / averageHL
            else:
                CV = 0
        return CV

    def getAverageHL(self):
        """Returns the average half-life of the operon """
        if self._probeHalfLife == []:
            return 0
        else:
            averageHL = numpy.average(self._HLlogTransform())
        return 2 ** averageHL

    def getName(self):
        """Returns an operon name """
        return self._name
    
    def getReads(self):
        """Returns a tuple of all base reads in an operon. """
        return self._reads

    def getOrientation(self):
        """Returns the orientation of an operon."""
        return self._orientation

    def getLeftBound(self):
        """Returns the left boundary of an operon."""
        return self._leftBound

    def getRightBound(self):
        """Returns the right boundary of an operon."""
        return self._rightBound

    def getBaseSequence(self):
        """Returns the base sequence of an operon."""
        return self._baseSequence
    
    def setLeftBoundPrecision(self, precision):
        """set the precsion of the left boundary to either
        True or False. """
        self._leftBoundPrecision = precision

    def getLeftBoundPrecision(self):
        """return the precsion of the left boundary"""
        return self._leftBoundPrecision

    def setRightBoundPrecision(self, precision):
        """set the precsion of the right boundary to either
        True or False. """
        self._rightBoundPrecision = precision

    def getRightBoundPrecision(self):
        """return the precsion of the right boundary"""
        return self._rightBoundPrecision

    
##def main():
##    readsList = [800, 800, 400, 800, 400, 100, 200, 50, 100, 100]
##    seq1 = Sequence(readsList, (0, 5), 'gene1', '-')
##    seq2 = Sequence(readsList, (6, 9), 'gene2', '-')
##    probeHalfLifeList = [(2, "-", 5), (4, "-", 10)]
##    operon1 = Operon('operon1')
##    operon1.add(seq1)
##    operon1.add(seq2)
##    leftBound = seq1.getStart()
##    rightBound = seq2.getEnd()
##    operon1.setLeftBound(leftBound)
##    operon1.setRightBound(rightBound)
##    operon1.setReads(readsList)
##    operon1.setOrientation()
##    operon1.setBaseSequence()
##    operon1.setProbeHalfLife(probeHalfLifeList)
##    operon1.setLeftBoundPrecision(True)
##    operon1.setRightBoundPrecision(False)
##    
##    print "The operon name is:            ", operon1.getName()
##    print "The operon reads are:          ", operon1.getReads()
##    print "The operon orientation is:     ", operon1.getOrientation()
##    print "The operon left boundary is:   ", operon1.getLeftBound()
##    print "The operon right boundary is:  ", operon1.getRightBound()    
##    print "The operon's base sequence is: ", operon1.getBaseSequence()
##    print "The operon contains:           ", str(operon1)
##    print "The operon's length is:        ", len(operon1)
##    print "The operon's probe half-life list is: ", operon1.getProbeHalfLife()
##    print "The CV of the operon's probe half-life is: ", operon1.getHLLogCV()
##    print "The average read of the operon is: ", operon1.getAverageRead()
##    print "The precision of the right boundary is: ", operon1.getRightBoundPrecision()
##
##main()
##        
