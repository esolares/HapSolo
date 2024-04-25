import random, glob, os, argparse
import pandas as pd
import gzip as gz

def CalculateContigSizes(asmFileName):
    # contigsDict[contigname] = [contiglen,headerpos,startseqpos,endseqpos]
    fin = gz.open(asmFileName)
    lastPos = headerPos = fin.tell()
    totalLines = sum(1 for line in fin)
    fin.seek(lastPos)
    seqLen = 0
    seqName = ''
    lastPos = 0
    count = 0
    myContigSizeDict = dict()
    while count < totalLines:
        lastPos = headerPos = fin.tell()
        line = fin.readline().decode().replace('\n', '')
        count = count + 1
        if line[0:1] == '>':
            seqName = line.split(" ")[0].replace('>', '').replace('/', '_')
            lastPos = startPos = fin.tell()
            line = fin.readline().decode().replace('\n', '')
            count = count + 1
            while line[0:1] != '>' and line[0:1] != '':
                seqLen = seqLen + len(line)
                endPos = lastPos
                lastPos = fin.tell()
                line = fin.readline().decode().replace('\n', '')
                count = count + 1
            if line[0:1] == '>' or line[0:1] == '':
                myContigSizeDict[seqName] = [seqLen, headerPos, startPos, endPos]
                seqName = ''
                seqLen = 0
                count = count - 1
                fin.seek(lastPos)
    fin.close()
    return myContigSizeDict
