"""
File: geneBank.py
Define GeneBank class that groups all genes on the chromosome

Use a regular list and Sequence to accomplish the task.
"""
import csv
import math
import cPickle
from sequence import Sequence
from exceptions import StopIteration

class GeneBank(object):
    """Groups all genes on the chromosome in a particular strain"""

    def __init__(self, fileName = None):
        """Creates a new dictionary to hold sequences. If
        a filename is provided, loads the accounts from a file
        of pickled sequences. """
        self._sequences = []
        self._fileName = fileName
        self._size = 0
        if fileName != None:
            fileObj = open(fileName, 'r')
            while True:
                try:
                    seq = cPickle.load(fileObj)
                    self.add(seq)
                except EOFError:
                    fileObj.close()
                    break

        
    def __str__(self):
        """Returns the string rep of all genes """
        seqStr = ''
        for sequence in self._sequences:
            seqStr += sequence.getName() + ', '
        return seqStr[0:-2]

    def isEmpty(self):
        """Returns True if the GeneBank is empty, False otherwise."""
        return self._size == 0

    def __len__(self):
        """Returns the length of the GeneBank """
        return self._size
    
    def __iter__(self):
        """Iterates through the GeneBank. """
        i = 0
        while True:
            if i == self._size:
                raise StopIteration
            yield self._sequences[i]
            i += 1
    
    def add(self, sequence):
        """Insert a sequence at the end of the list."""
        self._sequences.append(sequence)
        self._size += 1

    def __getitem__(self, index):
        """Takes the index. Returns the sequence at that index position."""
        if index >= 0 and index < self._size:
            return self._sequences[index]
        return None
    
    def get(self, name):
        """Returns a sequence with the given name."""
        for seq in self._sequences:
            if seq.getName() == name:
                return seq
        return None

    def getIndex(self, name):
        for i in xrange(self._size):
            seq = self._sequences[i]
            if seq.getName() == name:
                return i
        return -1                

    def getPrevious(self, sequence):
        """Precondition: The sequence exists.
        Returns its upstream sequence """
        if not sequence in self._sequences:
            return None
        i = 0
        while i < self._size:
            seq = self._sequences[i]
            if seq == sequence:
                break
            i += 1
        if i == 0:
            return None
        else:
            return self._sequences[i-1]
        
    def getNext(self, sequence):
        """Precondition: The sequence exists.
        Returns its downstream sequence """
        if not sequence in self._sequences:
            return None
        i = 0
        while i < self._size:
            seq = self._sequences[i]
            if seq == sequence:
                break
            i += 1
        if i == self._size-1:
            return None
        else:
            return self._sequences[i+1]



    def save(self, fileName = None):
        """Saves pickled sequences to a file. The parameter
        allows the user to change filenames."""
        if fileName != None:
            self._fileName = fileName
        elif self._fileName == None:
            return
        fileObj = open(self._fileName, 'w')
        for sequence in self._sequences:
            cPickle.dump(sequence, fileObj)
        fileObj.close()
        

##def main():
##    with open("wt1.csv", 'rU') as readsFile:
##        reads = csv.reader(readsFile, delimiter='\t')
##        readsList = []
##        for row in reads:
##            readsList.append(float(row[2]))
##            
##    seq1 = Sequence(readsList, (0, 10), "test1")
##    seq2 = Sequence(readsList, (20, 100), "test2")
##    seq3 = Sequence(readsList, (300, 1000), "test3")
##    seq4 = Sequence(readsList, (40, 100), "test4")
##
##    geneBank = GeneBank()
##    geneBank.add(seq1)
##    geneBank.add(seq2)
##    geneBank.add(seq3)
##    geneBank.add(seq4)
##
##    print "The geneBank is empty:         ", geneBank.isEmpty()
##    print "The length of the geneBank is: ", len(geneBank)
##    print "The geneBank contains:         "
##    for seq in geneBank:
##        print seq
##    print "The second seq in the geneBank is: ", geneBank[1]
##    print "test1 is at the position: ", geneBank.getIndex("test1")
##    print "The sequence downstream of seq2 is: ", geneBank.getNext(seq2)
##    print "The sequence upstream of seq2 is: ", geneBank.getPrevious(seq2)
##    
##    
##
##main()
##    

