"""
File: operonBank.py

Defines an operonBank that holds operons on a chromosome.

A regular list is used to hold operons.

"""
import math
import cPickle
import csv
from sequence    import Sequence
from geneBank   import GeneBank
from operon      import Operon
##from operonJudge import operonJudge


##geneBank = GeneBank("geneBank_wt1.txt")
##
##with open("wt1.csv", 'rU') as readsFile:
##    reads = csv.reader(readsFile, delimiter='\t')
##    readsList = []
##    for row in reads:
##        readsList.append(float(row[2]))
        

class OperonBank(object):


    def __init__(self, fileName = None):
        """Construct an operonBank object.
        Import the object if a fileName is give."""
        self._operons = []
        self._fileName = fileName
        self._size = 0
        if fileName != None:
            fileObj = open(fileName, 'r')
            while True:
                try:
                    operon = cPickle.load(fileObj)
                    self.add(operon)
                except EOFError:
                    fileObj.close()
                    break

    def __str__(self):
        """Returns the string rep of all operons """
        operonStr = ''
        for operon in self._operons:
            operonStr += operon.getName() + ', '
        return operonStr[0:-2]

    def __len__(self):
        """Returns the length of the operonBank """
        return self._size

    def __iter__(self):
        """Iterates through the operonBank. """
        i = 0
        while True:
            if i == self._size:
                raise StopIteration
            yield self._operons[i]
            i += 1

    def __getitem__(self, index):
        """Takes the index. Returns the sequence at that index position."""
        if index < self._size:
            return self._operons[index]
        return None
    
    def getName(self):
        """Returns the name of the OperonBank."""
        return self._fileName

    def add(self, operon):
        """Insert an operon to the end of operonBank."""
        self._operons.append(operon)
        self._size += 1
        
    def save(self, fileName = None):
        """Saves pickled sequences to a file. The parameter
        allows the user to change filenames."""
        if fileName != None:
            self._fileName = fileName
        elif self._fileName == None:
            return
        fileObj = open(self._fileName, 'w')
        for operon in self._operons:
            cPickle.dump(operon, fileObj)
        fileObj.close()



##operonBank = operonFinder(geneBank, readsList)
##for operon in operonBank:
##    print str(operon)
