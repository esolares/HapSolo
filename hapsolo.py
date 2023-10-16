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

#mode = parser.add_mutually_exclusive_group(required=True)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
parser.add_argument('--mode', help='HapSolo run mode. 0 = Random walking, 1 = No optimization with defaults, 2 = Optimized walking, Default = 0', type=int, required=False)

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
parser.add_argument('-P', '--minPID', help='Restrict values of PID to be >= the value set here. Default = 0.2', type=int, required=False)
parser.add_argument('-Q', '--minQ', help='Restrict values of Q to be >= the value set here. Default = 0.2', type=int, required=False)
parser.add_argument('-R', '--minQR', help='Restrict values of QR to be >= the value set here. Cannot be 0. Default = 0.2', type=int, required=False)
parser.add_argument('--min', help='Minimum size of contigs for Primary assembly. Default = 1000', type=int, required=False)


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
mode = args.mode
myMinPID = args.minPID
myMinQPctMin = args.minQ
myMinQRPctMin = args.minQR
myMinContigSize = args.min
aimode = 0

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
if mode == None:
    mode = 0
elif mode == 1:
    bestnscores = 1
    customMinPID = 0.7
    customMinQPctMin = 0.7
    customMinQRPctMin = 0.7
if myMinPID == None:
    myMinPID = 0.2
if myMinQPctMin == None:
    myMinQPctMin = 0.2
if myMinQRPctMin == None:
    myMinQRPctMin = 0.2
elif myMinQRPctMin < 0.02:
    myMinQRPctMin = 0.02
    print('-R/--minQR set to a value less than 0.02. using 0.02 instead.')

if myMinContigSize == None or myMinContigSize < 0:
    myMinContigSize = 1000

# maxASMSize = 600 * 1000000

dumpscores = True
stepsize = 0.0001
buscotypes = ['C', 'S', 'D', 'F', 'M']
resolution = 0.0001
mypddf = pd.DataFrame()
missingrefcontigset = set()
qrycontigset = set()
allcontigsset = set()
smallcontigset = set()
busco2contigdict = dict()
contigs2buscodict = dict()
pythonversion = sys.version_info[0]
# special_chars are !@#$%^&*()-=+,./\[{}]|;:'><?
special_chars = '!@#$%^&*()-=+,./\\[{]}|;:"\'><?'

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
            header = line[1:]
            special_char = False
            for char in header:
                if char in special_chars:
                    special_char = True
                    break
            # print('found seq_name ' + line)
            if len(header.split(" ")) > 1 or special_char:
                print('Spaces found in contig headers. Please remove spaces from contig names before proceeding with any analysis. Spaces, -"s, //"s and other special characters are not allowed in contig names.')
                print('Special characters except _ cause isues in aligners and BUSCO analysis. This will cause HapSolo to fail, as .')
                quit(1)
            seqName = header.split(" ")[0].replace('/', '_')
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
        mylinecounter = 0
        for line in open(mybuscofiles[i]):
            mylinecounter+=1
            if line[0] == '#' and mylinecounter < 4:
                if mylinecounter == 3:
                    contignames.add(line.split()[8].split('/')[-1].replace('.fasta',''))
            elif mylinecounter > 3:
                break
        #contignames.add(mybuscofiles[i].split('/')[-1].replace('full_table_', '').split('_new')[0])
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
        mylinecounter = 0
        for line in open(file):
            if line[0] != '#':
                mylines.append(line.strip().split())
            elif line[0] == '#' and mylinecounter < 4:
                mylinecounter+=1
                if mylinecounter == 3:
                    contigname = line.split()[8].split('/')[-1].replace('.fasta','')
        for i in range(0, len(mylines)):
            buscoid = mylines[i][0]
            #contigname = file.split('/')[-1].replace('full_table_', '').split('_new')[0]
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
    mycontigset = set(mycontigslist).union(missingrefcontigset) - smallcontigset
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
    temppd0 = temppd1[temppd1['QRAlignLenPct'] >= myQRPctMin]
    temppd1 = temppd0[temppd0['QRAlignLenPct'] <= myQRPctMax]
    goodcontigset = allcontigsset - set(temppd1['qName'])
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
    allmycontigs = qrycontigset.union(missingrefcontigset) - smallcontigset - {''}
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
    bestpurgedset = allcontigsset - bestcontigset  - {''}
    # bestscore = oldasmscorefxn
    bestbuscos = allcontigsbuscoscore.copy()
    # use fxn uniquepriorityqueue(pqlist, myvalues) for returning a sorted unique priority list
    myvalues = [oldasmscorefxn, bestcontigset, bestpurgedset, bestbuscos, [0.0, 0.0, 0.0]]
    # uniquepriorityqueue(bestnscoreslist[bestnscoreidx], [score, setofgoodcontigs, setofmissingcontigs, listofparameters, buscos])
    if mode != 1:
        bestnscoreslist.append(myvalues)
    # process:
    # 1. Make a step
    # 2. Calculate new assembly
    mygoodcontigs = ReduceASM(myPID, myQPctMin, myQRPctMin)
    mygoodcontigs = mygoodcontigs.union(missingrefcontigset) - smallcontigset - {''}
    purgedcontigs = allcontigsset - mygoodcontigs - {''}
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
    costfxndelta[0] = newasmscorefxn
    myvalues = [newasmscorefxn, mygoodcontigs, purgedcontigs, mygoodcontigsbuscoscore, [myPID, myQPctMin, myQRPctMin]]
    bestnscoreslist = uniquepriorityqueue(bestnscoreslist, myvalues)
    #if useprimaryformula:
    #    newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
    #else:
    #    newasmscorefxn = myLinearFxn(newmissingbuscos, newsinglebuscos, newdupebuscos, newfragbuscos, totalbuscos)
    # 4. use cost function of previous busco scores with new busco scores
    # 5. make new step trying to minimize cost function of busco scores.
    if mode == 1:
        return [bestnscoreslist, costfxn, costfxndelta]
    # Implement random forward walking, optimized decision based walking, other AI's here
    #if aimode == 0:
        # call rfw function
    #elif aimode == 1:
        # call odbw function
    #elif aimode == 2:
        # call another function
    #else:
        # throw an error
    for i in range(1, numofiterations):
        if costfxndelta[i] <= resolution and costfxndelta[i] > 0.0:
            break
        # forward stepping of GD
        if (myPID > 1.0 and myQPctMin > 1.0 and myQRPctMin > 1.0) or (i >= maxzeros and sum(costfxndelta[i-maxzeros:i+1]) == 0):
            # reassign myQ's
            myPID = uniform(myMinPID, 1.0)
            myQPctMin = uniform(myMinQPctMin, 1)
            myQRPctMin = uniform(myMinQRPctMin, 1)
        elif myQPctMin > 1.0 and myQRPctMin > 1.0:
            myQPctMin = uniform(myMinQPctMin, 1)
            myQRPctMin = uniform(myMinQRPctMin, 1)
        elif myPID > 1.0 and myQPctMin > 1.0:
            myPID = uniform(myMinPID, 1.0)
            myQPctMin = uniform(myMinQPctMin, 1)
        elif myPID > 1.0 and myQRPctMin > 1.0:
            myPID = uniform(myMinPID, 1.0)
            myQRPctMin = uniform(myMinQRPctMin, 1)
        elif myPID > 1.0:
            # reassign %ID
            myPID = uniform(myMinPID, 1.0)
        elif myQPctMin > 1.0:
            # reassign myQPctMin
            myQPctMin = uniform(myMinQPctMin, 1)
        elif myQRPctMin > 1.0:
            # reassign myQRPctMin
            myQRPctMin = uniform(myMinQRPctMin, 1)
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
        mygoodcontigs = mygoodcontigs.union(missingrefcontigset) - smallcontigset - {''}
        purgedcontigs = allcontigsset - mygoodcontigs - {''}
        numofcontigs = len(mygoodcontigs)
        mygoodcontigsbuscoscore = calculateBuscos(mygoodcontigs, busco2contigdict, contigs2buscodict)
        newsinglebuscos = mygoodcontigsbuscoscore['S']
        newmissingbuscos = mygoodcontigsbuscoscore['M']
        newdupebuscos = mygoodcontigsbuscoscore['D']
        newfragbuscos = mygoodcontigsbuscoscore['F']
        newasmsize = 0
        # cost function here
        if newsinglebuscos == 0:
            newasmscorefxn = 50000000.0
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
    if myPct < 0.02:
        inversePct = myPct
    else:
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
        #fin = gzip.open(alignmentfile, 'r')
        newalignfile = alignmentfile.replace('.paf.gz','.hap')
        if not os.path.exists(newalignfile):
            fout = open(newalignfile, 'w')
            fcounter = 0
            mcounter = 0
            for line in gzip.open(alignmentfile):
                line = line.strip().split()
                if len(line) < 11:
                    print('Error in reading PAF file')
                    quit(1)
                else:
                    fcounter+=1
                    #fout.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] + '\t' + line[9] + '\t' + line[10] + '\n')
                    fqAlignLen = max(int(line[2]), int(line[3])) - min(int(line[2]), int(line[3]))
                    frAlignLen = max(int(line[7]), int(line[8])) - min(int(line[7]), int(line[8]))
                    fQRAlignLenPct = CalculatePctAlign(fqAlignLen, frAlignLen)
                    fQPct = CalculatePctAlign(fqAlignLen, int(line[1]))
                    fPID = CalculatePctAlign(int(line[9]), fqAlignLen)
                    #fRPct = CalculatePctAlign(rAlignLen, int(line[6]))
                    if( (line[0] != line[5]) and (int(line[1] >= myMinContigSize)) and (fPID >= myMinPID) and (fQPct >= myMinQPctMin) and (fQRAlignLenPct >= myMinQRPctMin) ):
                        mcounter+=1
                        fout.write('"' + line[0] + '"' + '\t' + '"' + line[5] + '"' + '\t' + line[1] + '\t' + str(fQPct) + '\t' + str(fPID) + '\t' + str(fQRAlignLenPct) + '\n')
            print(str(fcounter - mcounter) + ' alignments Purged due to Search Space constraints')
            fout.close()
    elif fileext == 'paf':
        #fin = open(alignmentfile, 'r')
        newalignfile = alignmentfile.replace('.paf','.hap')
        if not os.path.exists(newalignfile):
            fcounter = 0
            mcounter = 0
            fout = open(newalignfile.replace('.paf','.hap'), 'w')
            for line in open(alignmentfile):
                line = line.strip().split()
                if len(line) < 11:
                    print('Error in reading PAF file')
                    quit(1)
                else:
                    fcounter+=1
                    #fout.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] + '\t' + line[9] + '\t' + line[10] + '\n')
                    fqAlignLen = max(int(line[2]), int(line[3])) - min(int(line[2]), int(line[3]))
                    frAlignLen = max(int(line[7]), int(line[8])) - min(int(line[7]), int(line[8]))
                    fQRAlignLenPct = CalculatePctAlign(fqAlignLen, frAlignLen)
                    fQPct = CalculatePctAlign(fqAlignLen, int(line[1]))
                    fPID = CalculatePctAlign(int(line[9]), fqAlignLen)
                    #fRPct = CalculatePctAlign(rAlignLen, int(line[6]))
                    if( (line[0] != line[5]) and (int(line[1] >= myMinContigSize)) and (fPID >= myMinPID) and (fQPct >= myMinQPctMin) and (fQRAlignLenPct >= myMinQRPctMin) ):
                        mcounter+=1
                        fout.write('"' + line[0] + '"' + '\t' + '"' + line[5] + '"' + '\t' + line[1] + '\t' + str(fQPct) + '\t' + str(fPID) + '\t' + str(fQRAlignLenPct) + '\n')
            print(str(fcounter - mcounter) + ' alignments Purged due to Search Space constraints')
            fout.close()
    #myLines = fin.readlines()
    #fin.close()
    #for lineNum in range(0, len(myLines)):
    fin = open(newalignfile)
    myline = fin.readline()
    if len(myline) <= 3:
        print('Empty HAP file. Please fix and rerun')
        quit(1)
    #myline = myline.strip().split('\t')
    #if len(myline) > 10:
    #     print('Invalid PAF format. Contains more than 18 fields. ' + str(
    #         len(myline)) + ' fields to be exact! Please correct.')
    #elif len(myline) != 17 and len(myline) != 18:
    #    print('Error in PAF file. expected 17 or 18 columns but received ' + str(len(myline)) + ' columns.')
    fin.close()
    # pandas time!
    #mypddf = pd.read_csv(newalignfile, sep='\t', header=None, names=['qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'matches', 'gaps+matches'])
    mypddf = pd.read_csv(newalignfile, sep='\t', header=None, names=['qName', 'tName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct'], dtype={'qName': object, 'tName': object})
    #if len(myline == 18:
        #mypddf = pd.DataFrame(myLines[:],
        #                  columns=['qName', 'qSize', 'qStart', 'qEnd', 'strand', 'tName', 'tSize', 'tStart', 'tEnd',
        #                           'matches', 'gaps+matches', 'mappingqv', 'alignmenttype', 'numofminschain',
        #                           'chainingscore', 'secondchainingscore', 'approxdivergence', 'lqrhrepseeds'])
    #elif len(myline) == 17:
        #mypddf = pd.DataFrame(myLines[:],
        #                  columns=['qName', 'qSize', 'qStart', 'qEnd', 'strand', 'tName', 'tSize', 'tStart', 'tEnd',
        #                           'matches', 'gaps+matches', 'mappingqv', 'alignmenttype', 'numofminschain',
        #                           'chainingscore', 'secondchainingscore', 'approxdivergence'])
    #else:
        #print('Error in PAF file. expected 17 or 18 columns but received ' + str(len(myLines[lineNum])) + ' columns.')
    #myLines = list() #clear this var to release RAM
    #mypddf['qStart'] = pd.to_numeric(mypddf['qStart'])
    #mypddf['qEnd'] = pd.to_numeric(mypddf['qEnd'])
    #mypddf['tStart'] = pd.to_numeric(mypddf['tStart'])
    #mypddf['tEnd'] = pd.to_numeric(mypddf['tEnd'])
    mypddf['qSize'] = pd.to_numeric(mypddf['qSize'])
    mypddf['QPct'] = pd.to_numeric(mypddf['QPct'])
    mypddf['PID'] = pd.to_numeric(mypddf['PID'])
    mypddf['QRAlignLenPct'] = pd.to_numeric(mypddf['QRAlignLenPct'])
    #mypddf['tSize'] = pd.to_numeric(mypddf['tSize'])
    #mypddf['matches'] = pd.to_numeric(mypddf['matches'])
    #mypddf['gaps+matches'] = pd.to_numeric(mypddf['gaps+matches'])
    # mypddf['qName'] = mypddf['qName'].str.replace('|','_').str[0:13]
    # mypddf['tName'] = mypddf['tName'].str.replace('|','_').str[0:13]
    # mypddf['qName'] = mypddf['qName'].str.split('|').str[0]
    # mypddf['tName'] = mypddf['tName'].str.split('|').str[0]
    #mypddf['qMin'] = mypddf[['qStart', 'qEnd']].min(axis=1)
    #mypddf['qMax'] = mypddf[['qStart', 'qEnd']].max(axis=1)
    #mypddf['tMin'] = mypddf[['tStart', 'tEnd']].min(axis=1)
    #mypddf['tMax'] = mypddf[['tStart', 'tEnd']].max(axis=1)
    #mypddf['qAlignLen'] = mypddf['qMax'] - mypddf['qMin']
    #mypddf['rAlignLen'] = mypddf['tMax'] - mypddf['tMin']
    #mypddf['QRAlignLenPct'] = mypddf[['qAlignLen', 'rAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['QPct'] = mypddf[['qAlignLen', 'qSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['PID'] = mypddf[['matches', 'qAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['RPct'] = mypddf[['rAlignLen', 'tSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    # Here we create a mask where qryname != refname
    #lenbeforemask = len(mypddf)
    #mypddf = mypddf[mypddf['qName'] != mypddf['tName']]
    #lenaftermask = len(mypddf)
    #print(str(lenbeforemask - lenaftermask) + ' alignments Purged where query = reference')
    #lenbeforemask = len(mypddf)
    #mypddf = mypddf[['qName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct']]
    #mypddf = mypddf[mypddf['qSize'] >= myMinContigSize]
    #mypddf = mypddf[mypddf['PID'] >= myMinPID]
    #mypddf = mypddf[mypddf['QPct'] >= myMinQPctMin]
    #mypddf = mypddf[mypddf['QRAlignLenPct'] >= myMinQRPctMin]
    #mypddf = mypddf[mypddf['QRAlignLenPct'] <= CalculateInverseProportion(myMinQRPctMin)]
    #lenaftermask = len(mypddf)
    #print(str(lenbeforemask - lenaftermask) + ' alignments Purged due to Search Space constraints')
    return mypddf


# Create a dictionary based on the alignment file
def CreateBlatAlignmentDataStruture(alignmentfile):
    global mypddf
    fileext = alignmentfile.split('.')[-1]
    if fileext == 'gz':
        newalignfile = alignmentfile.replace('.psl.gz','.hap')
        if not os.path.exists(newalignfile):
            fout = open(newalignfile, 'w')
            fcounter = 0
            mcounter = 0
            mylinenum = 0
            for line in gzip.open(alignmentfile):
                mylinenum += 1
                line = line.strip().split()
                if len(line) < 21 and mylinenum > 5:
                    print('Error in reading PSL file. Length of line < 21')
                    print(line)
                    quit(1)
                elif mylinenum > 5:
                    fcounter+=1
                    fqAlignLen = max(int(line[11]), int(line[12])) - min(int(line[11]), int(line[12]))
                    frAlignLen = max(int(line[15]), int(line[16])) - min(int(line[15]), int(line[16]))
                    fQRAlignLenPct = CalculatePctAlign(fqAlignLen, frAlignLen)
                    fQPct = CalculatePctAlign(fqAlignLen, int(line[10]))
                    fPID = CalculatePctAlign(int(line[0]), fqAlignLen)
                    #fRPct = CalculatePctAlign(rAlignLen, int(line[6]))
                    if( (line[9] != line[13]) and (int(line[10] >= myMinContigSize)) and (fPID >= myMinPID) and (fQPct >= myMinQPctMin) and (fQRAlignLenPct >= myMinQRPctMin) ):
                        mcounter+=1
                        fout.write('"' + line[9] + '"' + '\t' + '"' + line[13] + '"' + '\t' + line[1] + '\t' + str(fQPct) + '\t' + str(fPID) + '\t' + str(fQRAlignLenPct) + '\n')
            print(str(fcounter - mcounter) + ' alignments Purged due to Search Space constraints')
            fout.close()
    elif fileext == 'psl':
        #fin = open(alignmentfile, 'r')
        newalignfile = alignmentfile.replace('.psl','.hap')
        if not os.path.exists(newalignfile):
            fcounter = 0
            mcounter = 0
            mylinenum = 0
            fout = open(newalignfile.replace('.psl','.hap'), 'w')
            for line in open(alignmentfile):
                mylinenum += 1
                line = line.strip().split()
                if len(line) < 21 and mylinenum > 5:
                    print('Error in reading PSL file. Length of line < 21.')
                    print(line)
                    quit(1)
                elif mylinenum > 5:
                    fcounter+=1
                    fqAlignLen = max(int(line[11]), int(line[12])) - min(int(line[11]), int(line[12]))
                    frAlignLen = max(int(line[15]), int(line[16])) - min(int(line[15]), int(line[16]))
                    fQRAlignLenPct = CalculatePctAlign(fqAlignLen, frAlignLen)
                    fQPct = CalculatePctAlign(fqAlignLen, int(line[10]))
                    fPID = CalculatePctAlign(int(line[0]), fqAlignLen)
                    #fRPct = CalculatePctAlign(rAlignLen, int(line[6]))
                    if( (line[9] != line[13]) and (int(line[10] >= myMinContigSize)) and (fPID >= myMinPID) and (fQPct >= myMinQPctMin) and (fQRAlignLenPct >= myMinQRPctMin) ):
                        mcounter+=1
                        fout.write('"' + line[0] + '"' + '\t' + '"' + line[5] + '"' + '\t' + line[1] + '\t' + str(fQPct) + '\t' + str(fPID) + '\t' + str(fQRAlignLenPct) + '\n')
            print(str(fcounter - mcounter) + ' alignments Purged due to Search Space constraints')
            fout.close()
    #myLines = fin.readlines()
    #fin.close()
    #for lineNum in range(0, len(myLines)):
    fin = open(newalignfile)
    myline = fin.readline()
    if len(myline) <= 3:
        print('Empty HAP file. Please fix and rerun')
        quit(1)
    fin.close()
    # pandas time!
    #mypddf = pd.read_csv(newalignfile, sep='\t', header=None, names=['qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'matches', 'gaps$
    mypddf = pd.read_csv(newalignfile, sep='\t', header=None, names=['qName', 'tName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct'], dtype={'qName': object, 'tName': object})
    #for lineNum in range(0, len(myLines)):
    #    if myLines[lineNum][0] == '-':
    #        mystop = lineNum + 1
    #        break
    #if (mystop - 1) == len(myLines) or mystop == 0:
    #    print('Invalid PSL format. Missing ----\n. Please fix and rerun')
    #    quit(1)
    #mypddf = pd.DataFrame(myLines[mystop-1:],
    #                      columns=['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert',
    #                               'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName',
    #                               'tSize', 'tStart', 'tEnd'])
    #myLines = list() #clear this var to release RAM
    mypddf['qSize'] = pd.to_numeric(mypddf['qSize'])
    mypddf['QPct'] = pd.to_numeric(mypddf['QPct'])
    mypddf['PID'] = pd.to_numeric(mypddf['PID'])
    mypddf['QRAlignLenPct'] = pd.to_numeric(mypddf['QRAlignLenPct'])
    # mypddf['misMatches'] = pd.to_numeric(mypddf['misMatches'])
    # mypddf['repMatches'] = pd.to_numeric(mypddf['repMatches'])
    # mypddf['nCount'] = pd.to_numeric(mypddf['nCount'])
    # mypddf['qBaseInsert'] = pd.to_numeric((mypddf['qBaseInsert']))
    # mypddf['qName'] = mypddf['qName'].str.replace('|','_').str[0:13]
    # mypddf['tName'] = mypddf['tName'].str.replace('|','_').str[0:13]
    # mypddf['qName'] = mypddf['qName'].str.split('|').str[0]
    # mypddf['tName'] = mypddf['tName'].str.split('|').str[0]
    #mypddf['qMin'] = mypddf[['qStart', 'qEnd']].min(axis=1)
    #mypddf['qMax'] = mypddf[['qStart', 'qEnd']].max(axis=1)
    #mypddf['tMin'] = mypddf[['tStart', 'tEnd']].min(axis=1)
    #mypddf['tMax'] = mypddf[['tStart', 'tEnd']].max(axis=1)
    #mypddf['qAlignLen'] = mypddf['qMax'] - mypddf['qMin']
    #mypddf['rAlignLen'] = mypddf['tMax'] - mypddf['tMin']
    #mypddf['QRAlignLenPct'] = mypddf[['qAlignLen', 'rAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['QPct'] = mypddf[['qAlignLen', 'qSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['PID'] = mypddf[['matches', 'qAlignLen']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    #mypddf['RPct'] = mypddf[['rAlignLen', 'tSize']].apply(lambda x: CalculatePctAlign(*x), axis=1)
    # Here we create a mask where qryname != refname
    #lenbeforemask = len(mypddf)
    #mypddf = mypddf[mypddf['qName'] != mypddf['tName']]
    #lenaftermask = len(mypddf)
    #print(str(lenbeforemask - lenaftermask) + ' alignments Purged where query = reference')
    #lenbeforemask = len(mypddf)
    #mypddf = mypddf[['qName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct']]
    #mypddf = mypddf[mypddf['qSize'] >= myMinContigSize]
    #mypddf = mypddf[mypddf['PID'] >= myMinPID]
    #mypddf = mypddf[mypddf['QPct'] >= myMinQPctMin]
    #mypddf = mypddf[mypddf['QRAlignLenPct'] >= myMinQRPctMin]
    #mypddf = mypddf[mypddf['QRAlignLenPct'] <= CalculateInverseProportion(myMinQRPctMin)]
    #lenaftermask = len(mypddf)
    #print(str(lenbeforemask - lenaftermask) + ' alignments Purged due to Search Space constraints')
    return mypddf


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
    myGoodContigsSet = myGoodContigsSet - {''}
    # contigsDict[key] = [contiglen,headerpos,startseqpos,endseqpos]
    if len(myGoodContigsSet - set(myContigsDict.keys())) != 0:
        print('Error: HapSolo has two seperate set of contigs! Please submit bug report and sent bugreport.log file.')
        foutlogfile = open('bugreport.log','w')
        foutlogfile.write('Begin ContigsDict kyes:\n')
        for key in myContigsDict.keys():
            foutlogfile.write(str(key)) #  + ',' + str(myContigsDict[key][0]) + ',' + str(myContigsDict[key][1]) + ',' + str(myContigsDict[key][2]) + ',' + str(myContigsDict[key][3]) + '\n')
        foutlogfile.write('\nEnd ContigsDict keys\n\n')
        foutlogfile.write('Begin good contig set:\n')
        for contig in myGoodContigsSet:
            foutlogfile.write(str(contig) + ',')
        foutlogfile.write('\nEnd good contig set\n\n')
        foutlogfile.write('Begin non-matching contig set:\n')
        for contig in myGoodContigsSet - set(myContigsDict.keys()):
            foutlogfile.write(str(contig) + ',')
        foutlogfile.write('\nEnd non-matching contig set\n\n')
        foutlogfile.close()
        quit(1)
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
    for key in myContigsDict.keys():
        if myContigsDict[key][0] < myMinContigSize:
            smallcontigset.add(key)
    if pslalignmentfile == None:
        mypddf = CreateMM2AlignmentDataStructure(pafalignmentfile)
    elif pafalignmentfile == None:
        mypddf = CreateBlatAlignmentDataStruture(pslalignmentfile)
    qrycontigset = set(mypddf['qName'])
    missingrefcontigset = set(myContigsDict.keys()) - qrycontigset
    allcontigsset = set(myContigsDict.keys())
    busco2contigdict, contigs2buscodict = importBuscos(buscofileloc)
    # execute Hill Climbing here.
    job_args = list()
    if mode != 1:
        if threads == 1:
            job_args = [0, iterations, resolution, uniform(myMinPID, 1), uniform(myMinQPctMin, 1), uniform(myMinQRPctMin, 1)]
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
                job_args.append([i, iterations, resolution, uniform(myMinPID, 1), uniform(myMinQPctMin, 1), uniform(myMinQRPctMin, 1)])
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
    elif mode == 1:
        job_args = [0, 1, resolution, customMinPID, customMinQPctMin, customMinQRPctMin]
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
