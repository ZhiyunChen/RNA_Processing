SampleCode
==========
Documentation of SampleCodes

1. Brief introduction

| File name	| Function |
|---------------|  --------- |
| geneBank.py	| Defines the **GeneBank** class |
| geneBankConstruct.py|	Constructs a **GeneBank** object|
| HLShiftFinder.py|	Identifies processed genes |
| operon.py|	Defines the **Operon** class|
| operonBank.py|	Defines the **OperonBank** class|
| operonBankConstruct.py|	Constructs an **OperonBank** object|
| operonBankRefine.py|	Unifies different OperonBanks|
| RNaseYProcessedGeneFinder.py|	Identifies RNase Y processed genes|
| sequence.py|	Defines the **Sequence** class|

2. Procedure

2.1 Operon definition

2.1.1 GeneBank Construction

Running script:  **geneBankConstruct.py**

	Input data:
		wt1.csv
		wt2.csv
		rny1.csv
		rny2.csv
		Spy49_allGenes.csv
		SpyGenome.txt

	Output data:
		geneBank_wt1.txt
		geneBank_wt2.txt
		geneBank_rny1.txt
		geneBank_rny2.txt

2.1.2 OperonBank Construction

Running script: **operonBankConstruct.py**

	Input data:
		wt1.csv
		wt2.csv
		rny1.csv
		rny2.csv
		geneBank_wt1.txt
		geneBank_wt2.txt
		geneBank_rny1.txt
		geneBank_rny2.txt

	Output data:
		operonBank_wt1.txt
		operonBank _wt2.txt
		operonBank _rny1.txt
		operonBank _rny2.txt
		operonBank_wt1.bed
		operonBank _wt2.bed
		operonBank _rny1.bed
		operonBank _rny2.bed

2.1.3 Unification of four OperonBanks

Running script: **operonBankRefine.py**

	Input data:
		wt1.csv
		wt2.csv
		rny1.csv
		rny2.csv
		geneBank_wt1.txt
		geneBank_wt2.txt
		geneBank_rny1.txt
		geneBank_rny2.txt
		operonBank_wt1.txt
		operonBank _wt2.txt
		operonBank _rny1.txt
		operonBank _rny2.txt

	Output data:
		operonBankRD_wt1.txt
		operonBankRD _wt2.txt
		operonBankRD _rny1.txt
		operonBankRD _rny2.txt
		operonBankRD_wt1.bed
		operonBankRD _wt2.bed
		operonBankRD _rny1.bed
		operonBankRD _rny2.bed

2.2 Identification of RNase Y processed genes

2.2.1 Identification of processed genes in the WT and rny mutant

Running script: **HLShiftFinder.py**

	Input data:
		operonBankRD_wt1.txt
		operonBankRD _wt2.txt
		operonBankRD _rny1.txt
		operonBankRD _rny2.txt
		probehalfLifeWT.csv
		probehalfLifeRNY.csv

	Output data:
		processedGenes_wt1.csv
		processedGenes_wt2.csv
		processedGenes_rny1.csv
		processedGenes_rny2.csv

2.2.2 Identification of processed genes that are present in the WT but absent in the rny mutant

Running script: **RNaseYProcessedGeneFinder.py**

	Input data:
		processedGenes_wt1.csv
		processedGenes_rny1.csv
		operonBankRD_wt1.txt
		operonBankRD_rny1.txt

	Output data:
		rnyProcessedGenes.csv
