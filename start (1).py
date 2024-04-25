import random, glob, os, argparse
import pandas as pd
import gzip as gz
import cudf
import CalculateContigSizes

buscofileloc = 'anof_funestus_new_buscooutputs.tsv.gz'
hapalignmentfile = 'anof_funestus_new_self_aln.hap.gz'
myasmFileName = 'primary_new.fasta.gz'

n_best_sol = 5
maxzeros = 10
niterations = 5
weight_missing = 1
weight_duplicate = 1
weight_single = 1
weight_fragmented = 1

contigsDictSet = set()
contigsDictionary = dict()
buscosDictionary = dict() #key: busco, value = [count of complete, count of fragmented], each busco can only be added to one contig

mydf = cudf.read_csv(hapalignmentfile, sep='\t', header=None, names=['qName', 'tName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct'], dtype={'qName': object, 'tName': object})
myContigSizeDict = CalculateContigSizes.CalculateContigSizes(myasmFileName)
myAllContigsSet = set(myContigSizeDict)

for elem in myAllContigsSet:
    contigsDictionary[elem] = set()
BUSCOS2CTGSDICT = dict() 
for line in gz.open(buscofileloc): 
    line = line.decode()
    line = line.strip().split()

    if line[1][0] != 'M':
        if line[2] not in contigsDictionary.keys(): 
            contigsDictionary[line[2]] = set() 
        contigsDictionary[line[2]].add(line[0])
    if len(line) >= 1:
        myBUSCO = line[0]
        myBUSCOtype = line[1][0]
        if len(line) > 2:
            myCtg = line[2]
        else:
            myCtg = ''
        if myBUSCO not in BUSCOS2CTGSDICT.keys():
            BUSCOS2CTGSDICT[myBUSCO] = []
        if myBUSCOtype != 'M':
            BUSCOS2CTGSDICT[myBUSCO].append([myCtg, myBUSCOtype])

contigsDictSet = set(contigsDictionary.keys())



