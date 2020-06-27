#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import glob

mytsvlist = glob.glob("busco/busco*/run*/full*.tsv")

buscoidset = set()
buscoiddict = dict()
myflag = True

for tsvfile in mytsvlist:
    counter = 1
    totalbuscos = 0
    if myflag == True:
        for line in open(tsvfile):
            if line[0] != '#':
                buscoidset.add(line.split()[0])
            else:
                if counter == 2:
                    totalbuscos = int(line.split(' ')[-1].replace(')', ''))
                counter+=1
    if len(buscoidset) != 0:
        if len(buscoidset) == totalbuscos:
            myflag = False
            for buscoid in buscoidset:
                buscoiddict[buscoid] = list()
    for line in open(tsvfile):
        if line[0] != '#':
            # Buscoid, Status, Contig, Start, End, Score, Length
            buscoid = ''
            buscostatus = ''
            contigname = ''
            startpos = 0
            endpos = 0
            alnscore = 0.00
            bidlen = 0
            mylist = line.strip().split()
            buscolen = len(mylist)
            if buscolen == 2:
                buscoid = mylist[0]
                buscostatus = mylist[1]
            elif buscolen > 2:
                buscoid = mylist[0]
                buscostatus = mylist[1]
                contigname = mylist[2]
                startpos = int(mylist[3])
                endpos = int(mylist[4])
                alnscore = float(mylist[5])
                bidlen = int(mylist[6])
                buscoiddict[buscoid].append([contigname, startpos, endpos, alnscore, buscostatus, bidlen])

singlebuscolist = list()
nonsinglebuscolist = list()
missingbuscolist = list()
for buscoid in buscoidset:
    if len(buscoiddict[buscoid]) == 1:
        singlebuscolist.append(buscoid)
    elif len(buscoiddict[buscoid]) > 1:
        nonsinglebuscolist.append(buscoid)
    else:
        missingbuscolist.append(buscoid)

for buscoid in buscoidset:
    if len(buscoiddict[buscoid]) == 1:
        singlebuscolist.append(buscoid)

fragmentedbuscos = set()
completebuscos = set()
for buscoid in singlebuscolist:
    if buscoiddict[buscoid][0][4] != 'Complete':
        fragmentedbuscos.add(buscoid)
    else:
        completebuscos.add(buscoid)
