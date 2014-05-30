"""
File:
HLShiftFinder.py

Detects probe-level half-life shifts within an operon.

"""

import csv
import math
import numpy
from operon import Operon
from operonBank import OperonBank


operonBankList = [OperonBank("operonBankRD_wt1.txt"), \
                  OperonBank("operonBankRD_wt2.txt"), \
                  OperonBank("operonBankRD_rny1.txt"), \
                  OperonBank("operonBankRD_rny2.txt")]
nameList = ["wt1", "wt2", "rny1", "rny2"]
probeHLFileList = ["probehalfLifeWT.csv", "probehalfLifeWT.csv", \
                   "probehalfLifeRNY.csv", "probehalfLifeRNY.csv"]

##Construct a probe-level half-life list
for i in xrange(0, 4):
    fileName = "processedGenes_" + nameList[i] + ".csv"
    processedGenes = open(fileName, 'w')
    probeHLFile = open(probeHLFileList[i], 'rU')
    probeHLTable = csv.reader(probeHLFile)
    probeHLList = []
    for row in probeHLTable:
        position = int(row[0])
        orientation = row[1]
        probeHL = float(row[2])
        probeHLList.append((position, orientation, probeHL))

    """Look for regions with differential RNA stabilities,
    using a sliding-window algorithm for averageRead change."""    
    count = 0
    operonBank = operonBankList [i]
    for operon in operonBank:
        operon.setProbeHalfLife(probeHLList)

        HLLogList = []
        for HLTuple in operon.getProbeHalfLife():
            halfLife = HLTuple[2]
            if halfLife > 0:
                HLLog = math.log(halfLife, 2)
                HLLogList.append(HLLog)
        if HLLogList != []:
            refAverageLogHL = numpy.average(HLLogList)
            workWidth = 3
            threshold = 1.0
            
            if len(HLLogList) > workWidth:
                start = 0
                end   = start + workWidth
                maxAverageLogHL = refAverageLogHL
                minAverageLogHL = refAverageLogHL
                while end < len(HLLogList):            
                    workSpace = HLLogList[start:end]
                    workAverageLogHL = numpy.average(workSpace)
                    if operon.getName() == "operon_118":
                        print "The work AverageLogHL is: ", round(workAverageLogHL, 2)
                    if workAverageLogHL > maxAverageLogHL:
                        maxAverageLogHL = workAverageLogHL
                    if workAverageLogHL < minAverageLogHL:
                        minAverageLogHL = workAverageLogHL 
                    start += 1
                    end    = start + workWidth
                
                if maxAverageLogHL >= refAverageLogHL + threshold \
                   or minAverageLogHL <= refAverageLogHL - threshold \
                   and operon.getAverageRead() > 500:
                    processedGenes.write(operon.getName() + \
                                         "," + operon[0].getName() + \
                                         "," + operon[-1].getName() + \
                                         "," + str(operon.getAverageRead()) + \
                                         "," + str(refAverageLogHL) + \
                                         "," + str(maxAverageLogHL) + \
                                         "," + str(minAverageLogHL) + \
                                         '\n')
                    count += 1
                    print str(operon)
                    print "The average RNA-seq read is: ", round(operon.getAverageRead())
    ##              print HLLogList
                    print "The ref AverageLogHL is:     ", round(refAverageLogHL, 2)
                    print "The max AverageLogHL is:     ", round(maxAverageLogHL, 2)
                    print "The min AverageLogHL is:     ", round(minAverageLogHL, 2), "\n"
    print count
    print fileName, " is generated" 
    
    operonBank.save()
        
            
        
