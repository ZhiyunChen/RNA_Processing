"""
File: sequence.py

Define a Sequence class to describe a DNA sequence.
This time the readsList is a must.
The sequence name is optional.

"""

##SpyGenome = 'ATGCCTGAAAGTCTAGCCTTAAG'

import math
import numpy

with open("SpyGenome.txt", 'r') as myfile:
    SpyGenome = myfile.read().replace('\n', '')

class Sequence(object):
    """ Defines a sequence class representing a string of bases"""

    SEQUENCESTR = SpyGenome

    def __init__(self, readsList, (start, end), name = '', orientation = '+'):
        """Constructor of a sequence object.
        readsList: a list containing RNA-seq reads in order.
        (start, end): a tuple containing two integers that indicate the
                start and end site of a sequence
        name: Optional, a string representing the name of the sequence
        orientation: Optional, a string ('+' or '-') representing the sequence
                orientation on the chromosome. The default is '+'. """

        self._name         = name
        self._start        = start
        self._end          = end
        self._position     = (start, end)
        self._orientation  = orientation
        self._reads        = readsList[start:end]
        self._baseSequence = Sequence.SEQUENCESTR[start:end]
        self._halfLife     = 0
        self._probeHalfLife = []

    def __str__(self):
        """Returns a string representation of the sequence."""
        seqStr = ''
        seqStr += "Sequence Name: " + self._name        + '\n'        
        seqStr += "Position:      " + str(self._start)  + ", " \
                  + str(self._end)  + '\n'        
        seqStr += "Orientation:   " + self._orientation + '\n'       
        seqStr += "Base sequence: " + self._baseSequence        
        return seqStr

    def getName(self):
        """Returns the name of the sequence"""
        return self._name

    def getPosition(self):
        """Returns a tuple showing the start and end sites of the sequence
        on the chromosome, regardless of sequence orientation."""
        return self._position
    
    def getStart(self):
        """Returns the start site of the sequence on the chromosome,
        regardless of sequence orientation."""
        return self._start
    
    def getEnd(self):
        """Returns the end site of the sequence on the chromosome,
        regardless of sequence orientation."""
        return self._end
    
    def getOrientation(self):
        """Returns '+' or '-' indicating that the sequence is oriented positive
        or negative on the chromosome."""
        return self._orientation

    def __len__(self):
        """Returns the number of bases in the sequence"""        
        return len(self._reads)

    def getBaseSequence(self):
        """Returns the base sequences of the sequence."""        
        return self._baseSequence
    
    def isEmpty(self):
        """Returns True if the sequence is empty, False otherwise."""
        return self._reads == ()

    def getReads(self):
        """Returns a tuple of all base reads in a sequence """
        return self._reads

    def getReadIndex(self, read):
        """Returns the location index of a read. """
        if read in self._reads:
            return self._reads.index(read)
        else:
            return -1

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
        """Returns the averaged base reads in the sequence. The reads are initially
        log-transformed, averaged, and then reverse-log-transformed."""
        if self._reads == ():
            return 0
        else:
            averageLogRead = numpy.average(self._logTransform())
            averageRead = 2 ** averageLogRead
            return averageRead

    def getLogCV(self):
        """Returns the coefficient variation of a sequence reads"""
        if self._reads == ():
            return 0
        else:
            average = numpy.average(self._logTransform())
            std     = numpy.std(self._logTransform())
            if average > 0:
                CV = std / average
            else:
                CV = 0
        return CV
        

    def getMinRead(self):
        """Returns the minimum read of the sequence."""
        if self._reads == () :
            raise IndexError, "The read list is empty."
        minRead = self._reads[0]
        for read in self._reads:
            if minRead > read:
                minRead = read
        return minRead

    def getMinReadLocation(self):
        """Returns the location of the minimum read on the chromosome. """
        i = 0
        j = 0
        while i < len(self):
            if self._reads[j] > self._reads[i]:
                j = i
            i += 1
        return self._start + j

    def getRightMinReadLocation(self):
        """Returns the right most location of the minimum read on the chromosome. """
        i = 0
        j = 0
        while i < len(self):
            if self._reads[j] >= self._reads[i]:
                j = i
            i += 1
        return self._start + j 
    

    def getMaxRead(self):
        """Returns the maximum read of the sequence."""
        if self._reads == ():
            raise IndexError, "The read list is empty."
        maxRead = self._reads[0]
        for read in self._reads:
            if maxRead < read:
                maxRead = read
        return maxRead

    def getHalfLife(self):
        return self._halfLife

    def setHalfLife(self, halfLifeDict):
        self._halfLife = halfLifeDict[self._name]
        

    def setProbeHalfLife(self, probeHalfLifeList):
        """Set probe half-lives to the sequence """
        for i in xrange(len(probeHalfLifeList)):
            position = probeHalfLifeList[i][0]
            orientation = probeHalfLifeList[i][1]
            if position >= self._start and \
               position < self._end and \
               orientation == self._orientation:
                self._probeHalfLife.append(probeHalfLifeList[i])
                
    def getProbeHalfLife(self):
        return self._probeHalfLife

    
##def main():
##    readsList = [800, 800, 400, 800, 400, 100, 200, 50, 100, 100]
##    sequence1 = Sequence(readsList, (0, 10), name = 'sequence1')
##    probeHalfLifeList = [(2, "+", 5), (4, "+", 10)]
##    halfLifeDict = {'sequence1': 2.5, 'sequence2': 5}
##    sequence1.setProbeHalfLife(probeHalfLifeList)
##    sequence1.setHalfLife(halfLifeDict)
##    print "The sequence is empty:           ", sequence1.isEmpty()
##    print "The sequence name is:            ", sequence1.getName()
##    print "The sequence position is:        ", sequence1.getPosition()
##    print "The sequence orientation is:     ", sequence1.getOrientation()
##    print "The sequence start site is:      ", sequence1.getStart()
##    print "The sequence end site is:        ", sequence1.getEnd()
##    print "The length of the sequence is:   ", len(sequence1)
##    print "The sequence reads are:          ", sequence1.getReads()
##    print "The read 200 locates at:         ", sequence1.getReadIndex(200)    
##    print "The sequence average read is:    ", sequence1.getAverageRead()
##    print "The sequence maximum read is:    ", sequence1.getMaxRead()
##    print "The sequence minimum read is:    ", sequence1.getMinRead()
##    print "The sequence minimum read locates at: ", sequence1.getMinReadLocation()
##    print "The sequence's base sequence is: ", sequence1.getBaseSequence()
##    print "The sequence is:                 ", str(sequence1)
##    print "The sequence is empty:           ", sequence1.isEmpty()
##    print "The sequence CV is:              ", sequence1.getLogCV()
##    print "The sequence half life is:       ", sequence1.getHalfLife()
##    print "The sequence probe half life is: ", sequence1.getProbeHalfLife()
##
##main()
##        

