#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import argparse, glob, gzip, os, datetime, sys
from math import exp, log, ceil
from random import seed, randint, uniform
import pandas as pd
import multiprocessing as mp

# usage haplotigreduction.py mypslfile.psl myfastafile.fasta buscoresults.tsv
parser = argparse.ArgumentParser(description='Process alignments and BUSCO"s for selecting reduced assembly candidates', epilog='-p/--psl and -a/--paf are mutually exclusive')
parser.add_argument('-i', '--input', help='Input Fasta file', type=str, required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
group.add_argument('-a', '--paf', help='Minimap2 PAF alignnment file. Note. paf file functionality is currently experimental', type=str)
parser.add_argument('-b', '--buscos', help='Location BUSCO output directories. i.e. buscoN/', type=str, required=True)
parser.add_argument('-m', '--maxzeros', help='Max # of times cost function delta can consecutively be 0. Default = 10', type=str, required=False)
parser.add_argument('-t', '--threads', help='# of threads. Multiplies iterations by threads. Default = 1', type=int, required=False)
parser.add_argument('-n', '--niterations', help='# of total iterations to run per gradient descent. Default = 1000', type=int, required=False)
parser.add_argument('-B', '--Bestn', help='# of best candidate assemblies to return using gradient descent. Default = 1', type=int, required=False)
parser.add_argument('-S', '--thetaS', help='Weight for single BUSCOs in linear fxn. Default = 1.0', type=int, required=False)
parser.add_argument('-D', '--thetaD', help='Weight for duplicate BUSCOs in linear fxn. Default = 1.0', type=int, required=False)
parser.add_argument('-F', '--thetaF', help='Weight for fragmented BUSCOs in linear fxn. Default = 0.0', type=int, required=False)
parser.add_argument('-M', '--thetaM', help='Weight for missing BUSCOs in linear fxn. Default = 1.0', type=int, required=False)
# parser.add_argument('-T', '--thetaS', help='Weight for single BUSCOs in linear fxn. Default = 1.0', type=int, required=False)

args = parser.parse_args()

useprimaryformula = True

myasmFileName = args.input
pslalignmentfile = args.psl
pafalignmentfile = args.paf
buscofileloc = args.buscos
maxzeros = args.maxzeros
threads = args.threads
iterations = args.niterations
bestnscores = args.Bestn
thetaS = args.thetaS
thetaD = args.thetaD
thetaM = args.thetaM
thetaF = args.thetaF
# thetaT = args.thetaT

#if alignmentfile == None:
#    print('Please assign an alignment file. Either by using Minimap2 or Blat/PBlat')
#    print('This can be done by doing either --paf or --psl. Please only submit one file.')
#    quit(1)
if maxzeros == None:
    maxzeros = 10
if threads == None:
    threads = 1
if iterations == None:
    iterations = 1000
if bestnscores == None:
    bestnscores = 1
if thetaS == None:
    thetaS = 1.0
if thetaD == None:
    thetaD = 1.0
if thetaM == None:
    thetaM = 1.0
if thetaF == None:
    thetaF = 0.0
# if thetaT == None:
    # thetaT = 0.0

# alignmentfile = 'chardonnay_quiv2x_qm3x_fu_qm2xfu_pilon.self_blat.psl'
# myasmFileName = 'chardonnay_quiv2x_qm3x_fu_qm2xfu_pilon.fasta'

# maxASMSize = 600 * 1000000

dumpscores = True
stepsize = 0.0001
buscotypes = ['C', 'S', 'D', 'F', 'M']
resolution = 0.0001
mypddf = pd.DataFrame()
missingrefcontigset = set()
qrycontigset = set()
myMinContigSize = 1000
busco2contigdict = dict()
contigs2buscodict = dict()
pythonversion = sys.version_info[0]

if pythonversion != 2:
    print("Please run the correct version of Python. You are currently running Python " + str(pythonversion))
    print("HapSolo is compatible with Python 2")
    quit(1)

######################################
def CalculateContigSizes(asmFileName):
    # contigsDict[contigname] = [contiglen,headerpos,startseqpos,endseqpos]
    fin = open(asmFileName)
    lastPos = headerPos = fin.tell()
    totalLines = sum(1 for line in fin)
    fin.seek(lastPos)
    seqLen = 0
    seqName = ''
    lastPos = 0
    count = 0
    myContigSizeDict = dict()
    # print('begin for loop')
    while count < totalLines:
        # print('for loop executed')
        lastPos = headerPos = fin.tell()
        line = fin.readline().replace('\n', '')
        count = count + 1
        if line[0:1] == '>':
            # print('found seq_name ' + line)
            seqName = line.split(" ")[0].replace('>', '').replace('/', '_')
            # seqName = line.split("_")[0]
            lastPos = startPos = fin.tell()
            line = fin.readline().replace('\n', '')
            count = count + 1
            # print('begin while loop on seq ' + line)
            while line[0:1] != '>' and line[0:1] != '':
                seqLen = seqLen + len(line)
                endPos = lastPos
                lastPos = fin.tell()
                line = fin.readline().replace('\n', '')
                count = count + 1
            if line[0:1] == '>' or line[0:1] == '':
                myContigSizeDict[seqName] = [seqLen, headerPos, startPos, endPos]
                # print(len(seq_read.replace("\n", "")))
                seqName = ''
                seqLen = 0
                count = count - 1
                fin.seek(lastPos)
    fin.close()
    return myContigSizeDict


def calculateasmstats(bestcontigset):
    mycontiglist = list()
    for contig in bestcontigset:
        if contig in myContigsDict.keys():
            mycontiglist.append(myContigsDict[contig][0])
    mycontiglist.sort(reverse=True)
    largestcontig = mycontiglist[0]
    asmsize = sum(mycontiglist)
    topn50contigs = 0
    n50 = 0
    for i in range(len(mycontiglist)):
        n50 = mycontiglist[i]
        topn50contigs = topn50contigs + mycontiglist[i]
        if topn50contigs > asmsize / 2.0:
            break
    l50 = mycontiglist.index(n50)
    return asmsize, n50, l50, largestcontig


def importBuscos(buscofileloc):
    contignames = set()
    buscoids = set()
    mybuscofiles = glob.glob(buscofileloc + '/busco*/*/full_table_*.tsv')
    global busco2contigdict
    global contigs2buscodict
    # propogate busco ids into a set
    for line in open(mybuscofiles[0]):
        if line[0] != '#':
            buscoids.add(line.strip().split()[0])
    # propogate contig names into a set
    for i in range(0, len(mybuscofiles)):
        contignames.add(mybuscofiles[i].split('/')[-1].replace('full_table_', '').split('_new')[0])
    if len(contignames) != len(set(contignames)):
        print('duplicate contig names exist. Please fix contig names so that no duplicates exist and rerun HapSolo')
        quit(1)
    # propogate dictionaries
    for buscoid in buscoids:
        busco2contigdict[buscoid] = dict()
        for buscotype in buscotypes:
            busco2contigdict[buscoid][buscotype] = list()
    for contigname in contignames:
        contigs2buscodict[contigname] = dict()
        for buscotype in buscotypes:
            contigs2buscodict[contigname][buscotype] = list()
    # create a data structure for duplicate, single and fragmented busco id lookups
    # should start with contigs? or buscoids? maybe both
    for file in mybuscofiles:
        mylines = list()
        for line in open(file):
            if line[0] != '#':
                mylines.append(line.strip().split())
        for i in range(0, len(mylines)):
            buscoid = mylines[i][0]
            contigname = file.split('/')[-1].replace('full_table_', '').split('_new')[0]
            buscotype = mylines[i][1][0]
            if buscotype != 'M':
                busco2contigdict[buscoid][buscotype].append(contigname)
                contigs2buscodict[contigname][buscotype].append(buscoid)
    return busco2contigdict, contigs2buscodict


def calculateBuscos(mycontigslist, busco2contigdict, contigs2buscodict):
    # how should we deal with fragmented busco exists but exists as complete elsewhere?
    duplicatebuscos = 0
    singlebuscos = 0
    fragmentedbuscos = 0
    buscotypecounts = dict()
    buscoids = busco2contigdict.keys()
    # count the number of complete buscos
    completebuscoidcounts = dict()
    fragmentedbuscoidcounts = dict()
    for buscoid in buscoids:
        completebuscoidcounts[buscoid] = 0
        fragmentedbuscoidcounts[buscoid] = 0
    for buscotype in buscotypes:
        buscotypecounts[buscotype] = 0
    mycontigset = set(mycontigslist).union(missingrefcontigset)
    for contig in mycontigset:
        if contig in contigs2buscodict.keys():
            for buscotype in contigs2buscodict[contig]:
                buscosize = len(contigs2buscodict[contig][buscotype])
                if buscotype == 'C' and buscosize > 0:
                    for buscoid in contigs2buscodict[contig][buscotype]:
                        # print('Found complete busco: ' + buscoid)
                        completebuscoidcounts[buscoid] += 1
        # else:
            # print('contig: ' + contig + ' not found in contigs2busco dictionary.')
    for contig in mycontigset:
        if contig in contigs2buscodict.keys():
            for buscotype in contigs2buscodict[contig]:
                buscosize = len(contigs2buscodict[contig][buscotype])
                if buscosize > 0 and buscotype == 'F':
                    for buscoid in contigs2buscodict[contig][buscotype]:
                        if completebuscoidcounts[buscoid] == 0 and fragmentedbuscoidcounts[buscoid] == 0:
                            fragmentedbuscos += 1
    for buscoid in completebuscoidcounts:
        mybuscocount = completebuscoidcounts[buscoid]
        if mybuscocount == 1:
            singlebuscos += 1
        elif mybuscocount > 1:
            duplicatebuscos += 1
    buscotypecounts['D'] = duplicatebuscos
    buscotypecounts['S'] = singlebuscos
    buscotypecounts['C'] = singlebuscos + duplicatebuscos
    buscotypecounts['F'] = fragmentedbuscos
    buscotypecounts['M'] = len(completebuscoidcounts) - buscotypecounts['D'] - buscotypecounts['S'] - buscotypecounts['F']
    return buscotypecounts


def ReduceASM(myPID, myQPctMin, myQRPctMin):
    # mypddf.columns = Index([u'matches', u'misMatches', u'repMatches', u'nCount', u'qNumInsert',
    # u'qBaseInsert', u'tNumInsert', u'tBaseInsert', u'strand', u'qName',
    # u'qSize', u'qStart', u'qEnd', u'tName', u'tSize', u'tStart', u'tEnd',
    # u'blockCount', u'blockSizes', u'qStarts', u'tStarts', u'qMin', u'qMax',
    # u'tMin', u'tMax', u'qAlignLen', u'rAlignLen', u'QRAlignLenPct', u'QPct',
    # u'QRPct'], dtype='object')
    myQRPctMax = CalculateInverseProportion(myQRPctMin)
    temppd0 = mypddf[mypddf['PID'] >= myPID]
    temppd1 = temppd0[temppd0['QPct'] >= myQPctMin]
    temppd0 = temppd1[temppd1['qSize'] >= myMinContigSize]
    temppd1 = temppd0[temppd0['QRAlignLenPct'] >= myQRPctMin]
    temppd0 = temppd1[temppd1['QRAlignLenPct'] <= myQRPctMax]
    goodcontigset = set(mypddf['qName']) - set(temppd0['qName'])
    return goodcontigset


def hillclimbing(job_args):
    mythread = job_args[0]
    numofiterations = job_args[1]
    # todo: remove resolution as it is already a global var. adjust job args
    resolution = job_args[2]
    myPID = job_args[3]
    myQPctMin = job_args[4]
    myQRPctMin = job_args[5]
    # parameters
    bestnscoreslist = list()
    mysteps = list()
    mysteps.append(0)
    mysteps.append(0)
    mysteps.append(0)
    # begin to populate data structure for best scores
    print('Starting threshholds. Thread: ' + str(mythread) + ' PID: ' + str(myPID) + ' QPctMin: ' + str(myQPctMin) + ' QRPctMin: ' + str(myQRPctMin))
    costfxn = [0.0 for ii in range(numofiterations)]
    costfxndelta = [0.0 for ii in range(numofiterations)]
    # 0. Calculate scores for original assembly
    allmycontigs = qrycontigset.union(missingrefcontigset)
    allcontigsbuscoscore = calculateBuscos(allmycontigs, busco2contigdict, contigs2buscodict)
    allsinglebuscos = allcontigsbuscoscore['S']
    allmissingbuscos = allcontigsbuscoscore['M']
    alldupebuscos = allcontigsbuscoscore['D']
    allfragbuscos = allcontigsbuscoscore['F']
    totalbuscos = allcontigsbuscoscore['C'] + allcontigsbuscoscore['M'] + allcontigsbuscoscore['F']
    # allasmsize = 0
    if allsinglebuscos == 0:
        oldasmscorefxn = 5000.0
    else:
        oldasmscorefxn = myLinearFxn(allmissingbuscos, allsinglebuscos, alldupebuscos, allfragbuscos, totalbuscos)
    bestcontigset = allmycontigs.copy()
    bestpurgedset = set()
    # bestscore = oldasmscorefxn
    bestbuscos = allcontigsbuscoscore.copy()
    # use fxn uniquepriorityqueue(pqlist, myvalues) for returning a sorted unique priority list
    myvalues = [oldasmscorefxn, bestcontigset, bestpurgedset, bestbuscos, [0.0, 0.0, 0.0]]
    # uniquepriorityqueue(bestnscoreslist[bestnscoreidx], [score, setofgoodcontigs, setofmissingcontigs, listofparameters, buscos])
    bestnscoreslist.append(myvalues)
    # process:
    # 1. Make a step
    # 2. Calculate new assembly
    mygoodcontigs = ReduceASM(myPID, myQPctMin, myQRPctMin)
    mygoodcontigs = mygoodcontigs.union(missingrefcontigset)
    purgedcontigs = qrycontigset - mygoodcontigs
    numofcontigs = len(mygoodcontigs)
    # 3. Calculate new busco scores
    mygoodcontigsbuscoscore = calculateBuscos(mygoodcontigs, busco2contigdict, contigs2buscodict)
    newsinglebuscos = mygoodcontigsbuscoscore['S']
    newmissingbuscos = mygoodcontigsbuscoscore['M']
    newdupebuscos = mygoodcontigsbuscoscore['D']
    newfragbuscos = mygoodcontigsbuscoscore['F']
    newasmsize = 0
    if newsinglebuscos == 0:
        newasmscorefxn = 5000.0
    else:
        newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
    costfxn[0] = newasmscorefxn
    costfxndelta[0] = 1.0
    myvalues = [newasmscorefxn, mygoodcontigs, purgedcontigs, mygoodcontigsbuscoscore, [myPID, myQPctMin, myQRPctMin]]
    bestnscoreslist = uniquepriorityqueue(bestnscoreslist, myvalues)
    #if useprimaryformula:
    #    newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
    #else:
    #    newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
    # 4. use cost function of previous busco scores with new busco scores
    # 5. make new step trying to minimize cost function of busco scores.
    for i in range(1, numofiterations):
        if costfxndelta[i] <= resolution and costfxndelta[i] > 0.0:
            break
        # forward stepping of GD
        if (myPID > 1.0 and myQPctMin > 1.0 and CalculateInverseProportion(myQRPctMin) > 1.0) or (i >= maxzeros and sum(costfxndelta[i-maxzeros:i+1]) == 0):
            # reassign myQ's
            myPID = uniform(0.6, 1.0)
            myQPctMin = uniform(0.6, 1)
            myQRPctMin = uniform(0.6, 1)
        elif myQPctMin > 1.0 and CalculateInverseProportion(myQRPctMin) > 1.0:
            myQPctMin = uniform(0.6, 1)
            myQRPctMin = uniform(0.6, 1)
        elif myPID > 1.0 and myQPctMin > 1.0:
            myPID = uniform(0.6, 1.0)
            myQPctMin = uniform(0.6, 1)
        elif myPID > 1.0 and myQRPctMin > 1.0:
            myPID = uniform(0.6, 1.0)
            myQRPctMin = uniform(0.6, 1)
        elif myPID > 1.0:
            # reassign %ID
            myPID = uniform(0.6, 1.0)
        elif myQPctMin > 1.0:
            # reassign myQPctMin
            myQPctMin = uniform(0.6, 1)
        elif myQRPctMin > 1.0:
            # reassign myQRPctMin
            myQRPctMin = uniform(0.6, 1)
        else:
            mystepindex = randint(0, len(mysteps)-1)
            mysteps[mystepindex] = stepsize
            while True:
                # print(str(mysteps[0]) + ':' + str(mysteps[1]))
                for j in range(0, len(mysteps)):
                    myrand = randint(0, len(mysteps)-1)
                    mysteps[j] = stepsize * myrand
                if mysteps[0] != 0.0 or mysteps[1] != 0.0 or mysteps[2] != 0.0:
                    break
            myPID = myPID + mysteps[0]
            myQPctMin = myQPctMin + mysteps[1]
            myQRPctMin = myQRPctMin + mysteps[2]
        mygoodcontigs = ReduceASM(myPID, myQPctMin, myQRPctMin)
        # include missing contigs > 10Mb
        mygoodcontigs = mygoodcontigs.union(missingrefcontigset)
        purgedcontigs = qrycontigset - mygoodcontigs
        numofcontigs = len(mygoodcontigs)
        mygoodcontigsbuscoscore = calculateBuscos(mygoodcontigs, busco2contigdict, contigs2buscodict)
        newsinglebuscos = mygoodcontigsbuscoscore['S']
        newmissingbuscos = mygoodcontigsbuscoscore['M']
        newdupebuscos = mygoodcontigsbuscoscore['D']
        newfragbuscos = mygoodcontigsbuscoscore['F']
        newasmsize = 0
        # cost function here
        if newsinglebuscos == 0:
            newasmscorefxn = 5000.0
        else:
            newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
        costfxn[i] = newasmscorefxn
        costfxndelta[i] = costfxn[i-1] - costfxn[i]
        print('Thread: ' + str(mythread) + ' Iteration: ' + str(i) + ' PID: ' + str(myPID) + ' QPctMin: ' + str(myQPctMin) + ' QRPctMin: ' + str(myQRPctMin) + ' CostFxnDelta: ' + str(costfxndelta[i]) + ' ASMScoreFxn: ' + str(newasmscorefxn))
        if costfxndelta[i] < 0.0 and costfxndelta[i-1] < 0.0:
            while mysteps[0] != 0 or mysteps[1] != 0:
                for j in range(0, len(mysteps)):
                    myrand = randint(0, len(mysteps)-1)
                    mysteps[j] = stepsize * myrand
        # store new data into data structure unique priority queue
        myvalues = [newasmscorefxn, mygoodcontigs, purgedcontigs, mygoodcontigsbuscoscore,
                    [myPID, myQPctMin, myQRPctMin]]
        if len(bestnscoreslist) < bestnscores:
            bestnscoreslist = uniquepriorityqueue(bestnscoreslist, myvalues)
        else:
            if bestnscoreslist[bestnscores-1][0] >= newasmscorefxn:
                bestnscoreslist = uniquepriorityqueue(bestnscoreslist, myvalues)
    return [bestnscoreslist, costfxn, costfxndelta]


def CalculatePctAlign(myAlignLen, myTotalLen):
    if myTotalLen == 0:
        myFloat = 0.0000
    else:
        myFloat = float(myAlignLen) / float(myTotalLen)
    return myFloat


def CalculateInverseProportion(myPct):
    inversePct = exp(-1.0 * log(myPct, 2))
    return inversePct


def uniquepriorityqueue(pqlist, myvalue):
    pqlist = pqlist[:]
    pqlist.append(myvalue)
    pqlist = sorted(pqlist, key=lambda x: x[0])
    while True:
        #print(pqlist)
        mydupes = set()
        realdupes = list()
        mysamesizesets = list()
        mymasterdupelist = list()
        myfinaldupelist = list()
        for i in range(0, len(pqlist)):
            mysamesizesets.append([len(pqlist[i][1]), i])
        mysamesizesets = sorted(mysamesizesets)
        for i in range(1, len(mysamesizesets)):
            if mysamesizesets[i-1][0] == mysamesizesets[i][0]:
                mydupes.add(mysamesizesets[i-1][1])
                mydupes.add(mysamesizesets[i][1])
            else:
                if len(mydupes) > 0:
                    mymasterdupelist.append(mydupes)
                    mydupes = set()
        if len(mydupes) > 0:
            mymasterdupelist.append(mydupes)
        for i in range(0, len(mymasterdupelist)):
            templist = list(mymasterdupelist[i])
            for j in range(0, len(templist)):
                for k in range(j, len(templist)):
                    if j != k:
                        myfinaldupelist.append((templist[j], templist[k]))
        # here we remove the dupes. Keep lowest score
        for i in range(0, len(myfinaldupelist)):
            if pqlist[myfinaldupelist[i][0]][1] == pqlist[myfinaldupelist[i][1]][1]:
                realdupes.append(i)
        if len(realdupes) == 0:
            return pqlist[0:bestnscores]
        else:
            for i in range(0, len(realdupes)):
                # compare scores
                if pqlist[myfinaldupelist[i][0]][0] < pqlist[myfinaldupelist[i][1]][0]:
                    #print('remove ' + str(pqlist[myfinaldupelist[i][1]]))
                    pqlist = pqlist[0:myfinaldupelist[i][1]] + pqlist[myfinaldupelist[i][1] + 1:]
                elif pqlist[myfinaldupelist[i][0]][0] >= pqlist[myfinaldupelist[i][1]][0]:
                    #print('remove ' + str(pqlist[myfinaldupelist[i][0]]))
                    pqlist = pqlist[0:myfinaldupelist[i][0]] + pqlist[myfinaldupelist[i][0] + 1:]


def CreateMM2AlignmentDataStructure(alignmentfile):
    global mypddf
    fileext = alignmentfile.split('.')[-1]
    if fileext == 'gz':
        fin = gzip.open(alignmentfile, 'r')
    elif fileext == 'paf':
        fin = open(alignmentfile, 'r')
    myLines = fin.readlines()
    fin.close()
    for lineNum in range(0, len(myLines)):
        myLines[lineNum] = myLines[lineNum].strip().split('\t')
    if len(myLines) == 0:
        print('Empty PAF file. Please fix and rerun')
        quit(1)
    for lineNum in range(0, len(myLines)):
        if len(myLines[lineNum]) > 18:
            print('Invalid PAF format. Line Number: ' + str(lineNum + 1) + ' contains more than 18 fields. ' + str(
                len(myLines[lineNum])) + ' fields to be exact! Please correct.')
    # pandas time!
    mypddf = pd.DataFrame(myLines[:],
                          columns=['qName', 'qSize', 'qStart', 'qEnd', 'strand', 'tName', 'tSize', 'tStart', 'tEnd',
                                   'matches', 'gaps+matches', 'mappingqv', 'alignmenttype', 'numofminschain',
                                   'chainingscore', 'secondchainingscore', 'approxdivergence', 'lqrhrepseeds'])
    myLines = list() #clear this var to release RAM
    mypddf['qStart'] = pd.to_numeric(mypddf['qStart'])
    mypddf['qEnd'] = pd.to_numeric(mypddf['qEnd'])
    mypddf['tStart'] = pd.to_numeric(mypddf['tStart'])
    mypddf['tEnd'] = pd.to_numeric(mypddf['tEnd'])
    mypddf['qSize'] = pd.to_numeric(mypddf['qSize'])
    mypddf['tSize'] = pd.to_numeric(mypddf['tSize'])
    mypddf['matches'] = pd.to_numeric(mypddf['matches'])
    # mypddf['misMatches'] = pd.to_numeric(mypddf['misMatches'])
    # mypddf['repMatches'] = pd.to_numeric(mypddf['repMatches'])
    # mypddf['nCount'] = pd.to_numeric(mypddf['nCount'])
    # mypddf['qBaseInsert'] = pd.to_numeric((mypddf['qBaseInsert']))
    # mypddf['qName'] = mypddf['qName'].str.replace('|','_').str[0:13]
    # mypddf['tName'] = mypddf['tName'].str.replace('|','_').str[0:13]
    # mypddf['qName'] = mypddf['qName'].str.split('|').str[0]
    # mypddf['tName'] = mypddf['tName'].str.split('|').str[0]
    mypddf['qMin'] = mypddf[['qStart', 'qEnd']].min(axis=1)
    mypddf['qMax'] = mypddf[['qStart', 'qEnd']].max(axis=1)
    mypddf['tMin'] = mypddf[['tStart', 'tEnd']].min(axis=1)
    mypddf['tMax'] = mypddf[['tStart', 'tEnd']].max(axis=1)
    mypddf['qAlignLen'] = mypddf['qMax'] - mypddf['qMin']
    mypddf['rAlignLen'] = mypddf['tMax'] - mypddf['tMin']
    mypddf['QRAlignLenPct'] = mypddf[['qAlignLen', 'rAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['QPct'] = mypddf[['qAlignLen', 'qSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['PID'] = mypddf[['matches', 'qAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['RPct'] = mypddf[['rAlignLen', 'tSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    # Here we create a mask where qryname != refname
    lenbeforemask = len(mypddf)
    mypddf = mypddf[mypddf['qName'] != mypddf['tName']]
    lenaftermask = len(mypddf)
    print(str(lenbeforemask - lenaftermask) + ' alignments Purged where query = reference')
    return mypddf[['qName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct']]


# Create a dictionary based on the alignment file
def CreateBlatAlignmentDataStruture(alignmentfile):
    global mypddf
    counter = 0
    fileext = alignmentfile.split('.')[-1]
    if fileext == 'gz':
        fin = gzip.open(alignmentfile, 'r')
    elif fileext == 'psl':
        fin = open(alignmentfile, 'r')
    myLines = fin.readlines()
    fin.close()
    mystop = 0
    for lineNum in range(0, len(myLines)):
        if myLines[lineNum][0] == '-':
            mystop = lineNum + 1
            break
    if (mystop - 1) == len(myLines) or mystop == 0:
        print('Invalid PSL format. Missing ----\n. Please fix and rerun')
        quit(1)
    # pandas datastructre here
    #  columns = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts']
    for lineNum in range(mystop, len(myLines)):
        myLines[lineNum] = myLines[lineNum].strip().split('\t')
    for lineNum in range(mystop, len(myLines)):
        if len(myLines[lineNum]) > 21 and len(myLines[lineNum]) == 41:
            print('Invalid PSL format. Line Number: ' + str(lineNum+1) + ' contains more than 21 fields. ' + str(len(myLines[lineNum])) + ' fields to be exact! Please correct.')
            newlist = list()
            for extras in range(0, int(ceil(len(myLines[lineNum])/21.0))):
                if extras == 0:
                    templist = myLines[lineNum][0:17]
                    if len(templist) >= 17:
                        #print(len(templist))
                        newlist.append(templist)
                else:
                    templist = [myLines[lineNum][(21*extras)-1].split(',')[-1]] + myLines[lineNum][21*extras:21*extras+16]
                    if len(templist) >= 17:
                        #print(len(templist))
                        newlist.append(templist)
            myLines[lineNum] = newlist[0]
            for ex in range(1, extras+1):
                myLines.append(newlist[ex])
                counter += 1
        elif len(myLines[lineNum]) > 21 and myLines[lineNum] != 41:
            myLines[lineNum] = myLines[lineNum][0:17]
            print('Invalid PSL format. Line Number: ' + str(lineNum+1) + ' contains more than 21 fields. ' + str(len(myLines[lineNum])) + ' fields to be exact! Please correct.')
            counter += 1
        else:
            myLines[lineNum] = myLines[lineNum][0:17]
    print(str(counter) + ' number of extra lines found in ' + alignmentfile + '.')
    mypddf = pd.DataFrame(myLines[mystop-1:],
                          columns=['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
                                   'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName',
                                   'tSize', 'tStart', 'tEnd'])
    myLines = list() #clear this var to release RAM
    mypddf['qStart'] = pd.to_numeric(mypddf['qStart'])
    mypddf['qEnd'] = pd.to_numeric(mypddf['qEnd'])
    mypddf['tStart'] = pd.to_numeric(mypddf['tStart'])
    mypddf['tEnd'] = pd.to_numeric(mypddf['tEnd'])
    mypddf['qSize'] = pd.to_numeric(mypddf['qSize'])
    mypddf['tSize'] = pd.to_numeric(mypddf['tSize'])
    mypddf['matches'] = pd.to_numeric(mypddf['matches'])
    # mypddf['misMatches'] = pd.to_numeric(mypddf['misMatches'])
    # mypddf['repMatches'] = pd.to_numeric(mypddf['repMatches'])
    # mypddf['nCount'] = pd.to_numeric(mypddf['nCount'])
    # mypddf['qBaseInsert'] = pd.to_numeric((mypddf['qBaseInsert']))
    # mypddf['qName'] = mypddf['qName'].str.replace('|','_').str[0:13]
    # mypddf['tName'] = mypddf['tName'].str.replace('|','_').str[0:13]
    # mypddf['qName'] = mypddf['qName'].str.split('|').str[0]
    # mypddf['tName'] = mypddf['tName'].str.split('|').str[0]
    mypddf['qMin'] = mypddf[['qStart', 'qEnd']].min(axis=1)
    mypddf['qMax'] = mypddf[['qStart', 'qEnd']].max(axis=1)
    mypddf['tMin'] = mypddf[['tStart', 'tEnd']].min(axis=1)
    mypddf['tMax'] = mypddf[['tStart', 'tEnd']].max(axis=1)
    mypddf['qAlignLen'] = mypddf['qMax'] - mypddf['qMin']
    mypddf['rAlignLen'] = mypddf['tMax'] - mypddf['tMin']
    mypddf['QRAlignLenPct'] = mypddf[['qAlignLen', 'rAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['QPct'] = mypddf[['qAlignLen', 'qSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['PID'] = mypddf[['matches', 'qAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    mypddf['RPct'] = mypddf[['rAlignLen', 'tSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    # Here we create a mask where qryname != refname
    lenbeforemask = len(mypddf)
    mypddf = mypddf[mypddf['qName'] != mypddf['tName']]
    lenaftermask = len(mypddf)
    print(str(lenbeforemask - lenaftermask) + ' alignments Purged where query = reference')
    return mypddf[['qName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct']]


def myLinearFxn(mbusco, sbusco, dbusco, fbusco, cbusco):
    myValue = float(thetaF * fbusco + thetaD * dbusco + thetaM * mbusco) / float(thetaS * sbusco)
    # todo: add custom linear function option
    return myValue


def WriteNewAssembly(myasmFileName, newASMFileName, myGoodContigsSet):
    fin = open(myasmFileName, 'r')
    mydirectory = 'asms'
    outfile = mydirectory + '/' + newASMFileName
    if not os.path.exists(mydirectory):
        os.makedirs(mydirectory)
    fout = open(outfile, 'w')
    # contigsDict[key] = [contiglen,headerpos,startseqpos,endseqpos]
    for contig in myGoodContigsSet:
        myContigPositionsList = myContigsDict[contig]
        fin.seek(myContigPositionsList[1])  # extract headerpos
        fout.write(fin.readline())
        newPos = fin.tell()
        mySeq = fin.readline().replace('\n', '')
        while newPos != myContigPositionsList[3]:
            newPos = fin.tell()
            mySeq = mySeq + fin.readline().replace('\n', '')
        fout.write(mySeq + '\n')
    fout.close()


if __name__ == '__main__':
    seed(1)
    try:
        myContigsDict = CalculateContigSizes(myasmFileName)
    except:
        print('Invalid assembly file. Please check and try again')
        quit(1)
    if pslalignmentfile == None:
        mypddf = CreateMM2AlignmentDataStructure(pafalignmentfile)
    elif pafalignmentfile == None:
        mypddf = CreateBlatAlignmentDataStruture(pslalignmentfile)
    qrycontigset = set(mypddf['qName'])
    missingrefcontigset = set(myContigsDict.keys()) - qrycontigset
    busco2contigdict, contigs2buscodict = importBuscos(buscofileloc)
    # execute Hill Climbing here.
    job_args = list()
    if threads == 1:
        job_args = [0, iterations, resolution, uniform(0.0, 1), uniform(0.0, 1), uniform(0.0, 1)]
        mylist = hillclimbing(job_args)
        mybestnscoreslist = mylist[0]
        for i in range(0, bestnscores):
            for j in range(0, len(mybestnscoreslist[0][4])):
                mybestnscoreslist[i][4][j] = '%.4f' % mybestnscoreslist[i][4][j]
            newasmfilename = myasmFileName.replace('.fasta', '') + '_' + str(myMinContigSize) + '_' + str(
                mybestnscoreslist[i][4][0]) + '_' + str(mybestnscoreslist[i][4][2]) + 'to' + str(
                '%.4f' % CalculateInverseProportion(float(mybestnscoreslist[i][4][2]))) + '_' + str(
                mybestnscoreslist[i][4][1]) + '_primary.fasta'
            print('Writing ' + newasmfilename + ' with score: ' + str(mybestnscoreslist[i][0]))
            WriteNewAssembly(myasmFileName, newasmfilename, mybestnscoreslist[i][1])
            WriteNewAssembly(myasmFileName, newasmfilename.replace('_primary.fasta', '_secondary.fasta'), mybestnscoreslist[i][2])
        if dumpscores:
            fout = open(myasmFileName.replace('.fasta', '_' + str(datetime.datetime.today()).replace(' ', '_').replace('-', '_').replace(':', '_').split('.')[0] + '.scores'), 'w')
            fout.write(str(mylist[1][0]))
            for i in range(1, iterations):
                fout.write(',' + str(mylist[1][i]))
            fout.close()
            fout = open(myasmFileName.replace('.fasta', '_' + str(datetime.datetime.today()).replace(' ', '_').replace('-', '_').replace(':', '_').split('.')[0] + '.deltascores'), 'w')
            fout.write(str(mylist[2][0]))
            for i in range(1, iterations):
                fout.write(',' + str(mylist[2][i]))
            fout.close()
    elif threads > 1:
        for i in range(threads):
            job_args.append([i, iterations, resolution, uniform(0.0, 1), uniform(0.0, 1), uniform(0.0, 1)])
        pool = mp.Pool(processes=threads)
        mylist = pool.map(hillclimbing, job_args)
        mybestnscoreslist = list()
        mybestnscoreslist.append(mylist[0][0][0])
        for i in range(0, threads):
            for j in range(0, bestnscores):
                mybestnscoreslist = uniquepriorityqueue(mybestnscoreslist, mylist[i][0][j])
        for i in range(0, bestnscores):
            for j in range(0, len(mybestnscoreslist[0][4])):
                mybestnscoreslist[i][4][j] = '%.4f' % mybestnscoreslist[i][4][j]
            newasmfilename = myasmFileName.replace('.fasta', '') + '_' + str(myMinContigSize) + '_' + str(
                mybestnscoreslist[i][4][0]) + '_' + str(mybestnscoreslist[i][4][2]) + 'to' + str(
                '%.4f' % CalculateInverseProportion(float(mybestnscoreslist[i][4][2]))) + '_' + str(
                mybestnscoreslist[i][4][1]) + '_primary.fasta'
            print('Writing ' + newasmfilename + ' with score: ' + str(mybestnscoreslist[i][0]))
            WriteNewAssembly(myasmFileName, newasmfilename, mybestnscoreslist[i][1])
            WriteNewAssembly(myasmFileName, newasmfilename.replace('_primary.fasta', '_secondary.fasta'), mybestnscoreslist[i][2])
        if dumpscores:
            fout = open(myasmFileName.replace('.fasta', '_' + str(datetime.datetime.today()).replace(' ', '_').replace('-', '_').replace(':', '_').split('.')[0] + '.scores'), 'w')
            for i in range(0, threads):
                fout.write(str(mylist[i][1][0]))
                for j in range(1, iterations):
                    fout.write(',' + str(mylist[i][1][j]))
                fout.write('\n')
            fout.close()
            fout = open(myasmFileName.replace('.fasta', '_' + str(datetime.datetime.today()).replace(' ', '_').replace('-', '_').replace(':', '_').split('.')[0] + '.deltascores'), 'w')
            for i in range(0, threads):
                fout.write(str(mylist[i][2][0]))
                for j in range(1, iterations):
                    fout.write(',' + str(mylist[i][2][j]))
                fout.write('\n')
            fout.close()
    else:
        print('Invalid # of threads set. Please use a positive integer for threads')
        quit(1)
