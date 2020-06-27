#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import os, sys, glob, csv, argparse
#import pandas as pd

parser = argparse.ArgumentParser(description='Process alignments and BUSCO"s for selecting reduced assembly candidates')
parser.add_argument('-i', '--input', help='Input Fasta file', type=str, required=True)
parser.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str, required=True)
parser.add_argument('-b', '--buscos', help='The directory of all contig BUSCO full output tsv files', type=str, required=True)
args = parser.parse_args()

fullasmbuscofile = args.input
buscotsvdir = args.buscos

if buscotsvdir == None:
    buscotsvdir = ''

contigbuscofiles = glob.glob(buscotsvdir + "/full*.tsv")
mysetcids = set()
mydictcontigbuscos = dict()
mydictfullbuscos = dict()
#fin = open(contigbuscofiles[0])
#read until no # in first char
for file in contigbuscofiles:
    for line in open(file):
        if line[0] != '#':
            line = line.strip().split()
            if len(line) > 1:
                mysetcids.add(line[0])

for i in mysetcids:
    mydictcontigbuscos[i] = list()

mysetfids = set()
for file in fullbuscofiles:
    for line in open(file):
        if line[0] != '#':
            line = line.strip().split()
            if len(line) > 1:
                mysetfids.add(line[0])

for file in contigbuscofiles:
    for line in open(file):
        if line[0] != '#':
            line = line.strip().split()
            if len(line) > 2:
                mydictcontigbuscos[line[0]].append([line[2].replace('_','')] + [line[1]] + line[3:])

for i in mysetfids:
    mydictfullbuscos[i] = list()

for file in fullbuscofiles:
    for line in open(file):
        if line[0] != '#':
            line = line.strip().split()
            if len(line) > 2:
                mydictfullbuscos[line[0]].append([line[2].replace('_','')] + [line[1]] + line[3:])

if len(mysetcids - mysetfids) > 0:
    print("Error: More BUSCO's may exist in contigs than in full fasta file")
    quit()
elif len(mysetfids - mysetcids) > 0:
    print("Warning: More BUSCO's exist in full fasta file")
    print("Testing each individual BUSCO for extra information:")
    counter = 0
    for i in (mysetfids - mysetcids):
        if i in mydictfullbuscos.keys():
            if len(mydictfullbuscos[i]) > 0:
                counter+=1
                print("Please check the following BUSCO ID: " + str(i) + "on contig: " + mydictfullbuscos[i][0][0])
    if counter > 0:
        print("Test failed: Number of BUSCO's between contigs and full fasta are not equal. Please check your input files again and try again.")
        print("RAD-C will proceed.")
    else:
        print("Test Passed! Number of BUSCO's between contigs and full fasta are equal.")
#test
mydictfullbuscos['EOG090W0JFZ']
mydictfullbuscos['EOG090W0JMX']
