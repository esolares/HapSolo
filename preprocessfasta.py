#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import re, os, argparse

# usage preprocessfasta.py -i myfastafile.fasta
parser = argparse.ArgumentParser(description='Preprocess FASTA file and outputs a clean FASTA and seperates contigs based on unique headers. Removes special chars')
parser.add_argument('-i', '--input', help='Input FASTA file', type=str, required=True)
parser.add_argument('-m', '--maxcontig', help='Max size of contig in Mb to output in contigs folder.', type=str, required=False)

args = parser.parse_args()
asmfilename = args.input
maxcontigsize = args.maxcontig

if maxcontigsize == None:
    maxcontigsize = 10*1000000
else:
    maxcontigsize = 1*1000000

fin = open(asmfilename)
mylines = fin.readlines()
myheaderidxs = list()
for i in range(len(mylines)):
    mylines[i] = mylines[i].strip()
    if mylines[i][0] == '>':
        myheaderidxs.append(i)


maxlen = -1
for i in range(len(myheaderidxs)):
    mylines[myheaderidxs[i]] = re.sub('[^a-zA-Z0-9\n\.]', '_', mylines[myheaderidxs[i]][1:])
    mylen = len(mylines[myheaderidxs[i]])
    if mylen > maxlen:
        maxlen = mylen

uniqueheadersize = -1
for i in range(1, maxlen+1):
    mytempset = set()
    for j in range(len(myheaderidxs)):
        mytempset.add(mylines[myheaderidxs[j]][0:i])
    if len(mytempset) == len(myheaderidxs):
        uniqueheadersize = i
        break

if uniqueheadersize == -1:
    print('This FASTA file does not contain unique headers. Please fix and rerun again.')
    quit(1)

for i in range(len(myheaderidxs)):
    mylines[myheaderidxs[i]] = mylines[myheaderidxs[i]][0:uniqueheadersize]

mytempset = set()
for i in range(len(myheaderidxs)):
    mytempset.add(mylines[myheaderidxs[i]].split('_')[0])

if mytempset == len(myheaderidxs):
    mylines[myheaderidxs[i]] = mylines[myheaderidxs[i]].split('_')[0]

fileext = asmfilename.split('.')[-1]
fout = open(asmfilename.replace('.' + fileext, '') + '_new.' + fileext, 'w')
outdir = "contigs"
try:
    os.stat(outdir)
except:
    os.mkdir(outdir)

for j in range(len(myheaderidxs)):
    fout.write('>' + mylines[myheaderidxs[j]] + '\n')
    mystr = ''
    if j < len(myheaderidxs)-1:
        for i in range(myheaderidxs[j]+1, myheaderidxs[j+1]):
            mystr = mystr + mylines[i]
    else:
        for i in range(myheaderidxs[j]+1, len(mylines)):
            mystr = mystr + mylines[i]
    fout.write(mystr + '\n')
    if len(mystr) < maxcontigsize:
        contigout = open(outdir + '/' + mylines[myheaderidxs[j]] + '.fasta', 'w')
        contigout.write('>' + mylines[myheaderidxs[j]] + '\n')
        contigout.write(mystr + '\n')
        contigout.close()

fout.close()
