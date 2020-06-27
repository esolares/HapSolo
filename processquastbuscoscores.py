#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import os, sys, glob, csv
import pandas as pd

def writeCSVFile(myList):
    with open('haplotigreduction.csv', 'wb') as fout:
        writer = csv.writer(fout, lineterminator='\n')
        writer.writerows(myList)

offset = 1
fastaprefix = 'GCF_000238955.4_M_zebra_UMD2a_genomic'
prefix = fastaprefix
fastaprefixoffset = len(fastaprefix.split('_'))
asmloc = 10 + offset
n50loc = 20 + offset
buscoscoreloc = 32 + offset

myoutputfiles = glob.glob(prefix + '*.txt')
myscores = dict()
for file in myoutputfiles:
    counter = 1
    for line in open(file):
        if counter == asmloc:
            asmsize = int(line.strip().split()[-1])
        elif counter == n50loc:
            n50size = int(line.strip().split()[-1])
        elif counter == buscoscoreloc:
            buscoscores = line.strip()
            print('processing busco scores: ' + str(buscoscores) + 'in file: ' + str(file))
            cbusco = round(float(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[0][2:])/100, 3)
            sbusco = round(float(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[1][2:])/100, 3)
            fbusco = round(float(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[3][2:])/100, 3)
            dbusco = round(float(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[2][2:])/100, 3)
            mbusco = round(float(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[4][2:])/100, 3)
            nbuscos = int(buscoscores.replace('[', '').replace(']', '').replace(',', '').split('%')[5][2:])
        counter += 1
    myscores[file] = [asmsize, n50size, cbusco, sbusco, dbusco, fbusco, mbusco, nbuscos]

mypctset = set()
mymincontiglen = set()
myovlrangepct = set()
for file in myscores.keys():
    filelist = file.split('_')
    mymincontiglen.add(int(filelist[1 + fastaprefixoffset-1]))
    mypctset.add(float(filelist[3 + fastaprefixoffset-1]))
    myovlrangepct.add(filelist[2 + fastaprefixoffset-1])

mypctlist = sorted(list(mypctset))
mymincontiglenlist = sorted(list(mymincontiglen))
myovlrangepctlist = list(myovlrangepct)

asmsize = -1
n50size = -1
buscoscores = ''
#myfile = fastaprefix + str(mymincontiglenlist[0]) + '_' + myovlrangepctlist[0] + '_' + str(mypctlist[0]) + '_new_scoresreport.txt'
mymasterlist = list()
counter = 0
lcounter = 0
for mincontig in mymincontiglenlist:
    for ovlrange in myovlrangepctlist:
        for pct in mypctlist:
            myfile = fastaprefix + '_' + str(mincontig) + '_' + ovlrange + '_' + str(pct) + '_new_scoresreport.txt'
            counter += 1
            if myfile in myscores.keys():
                lcounter += 1
                #[asmsize, n50size, cbusco, sbusco, dbusco, fbusco, mbusco, nbuscos, mincontig, ovlrange, pct]
                mylist = myscores[myfile]
                mylist.append(mincontig)
                mylist.append(ovlrange)
                mylist.append(pct)
                mymasterlist.append(mylist)

writeCSVFile(mymasterlist)
mydf = pd.DataFrame(mymasterlist)
mydf.columns = ['asmsize', 'n50size', 'cbusco', 'sbusco', 'dbusco', 'fbusco', 'mbusco', 'nbuscos', 'mincontig', 'ovlrange', 'pct']
#sort by asmsize,n50,singlebusco,completebusco,nbuscos
#mydf.sort_values(['cbusco','sbusco','asmsize'],ascending=[False,False,False])
#mydf.sort_values(['mbusco','n50size','cbusco'],ascending=[True,False,False])
#mydf.sort_values(['mbusco','dbusco','cbusco','n50size','mincontig','asmsize'],ascending=[True,True,False,False,True,True])
#mydf.sort_values(['mbusco','dbusco','sbusco','asmsize'],ascending=[True,True,False,False]).head(6)

myfinalasmzerominctg = mydf[mydf['mincontig']==0].sort_values(['mbusco','dbusco','sbusco','asmsize'],ascending=[True,True,False,False]).head(1)
myfinalasm = mydf.sort_values(['mbusco','dbusco','sbusco','asmsize'],ascending=[True,True,False,False]).head(1)
#myprefix = 'hch.pilon2'
myprefix = 'GCF_000238955.4_M_zebra_UMD2a_genomic'
#print str(myfinalasm)
myfinalasmfasta = myprefix + '_' + str(myfinalasm['mincontig'].iat[0]) + '_' + str(myfinalasm['ovlrange'].iat[0]) + '_' + str(myfinalasm['pct'].iat[0]) + '_new.fasta'
print(myfinalasmfasta)
#myprefix = 'GCF_000238955.4_M_zebra_UMD2a_genomic'
'''
#extra stuff here
mypd = pd.read_csv('haplotigreduction.csv', header=None)
mypd.columns = ['asmsize', 'n50size', 'cbusco', 'sbusco', 'dbusco', 'fbusco', 'mbusco', 'nbuscos', 'mincontig', 'ovlrange', 'pct']
mypd.insert(11,'LinearFxn',((mypd['dbusco']+mypd['mbusco'])/mypd['sbusco']),True)
rank1 = mypd.sort_values(by=['LinearFxn'],ascending=True).head(20)
rank2 = mypd.sort_values(['mbusco','dbusco','sbusco','asmsize'],ascending=[True,True,False,False]).head(20)

myprefix = 'hch.pilon2'
myfinalasmfastasr1 = list()
for i in range(1,21):
    myfinalasm = rank1.head(i).tail(1)
    myfinalasmfastasr1.append(myprefix + '_' + str(myfinalasm['mincontig'].iat[0]) + '_' + str(myfinalasm['ovlrange'].iat[0]) + '_' + str(myfinalasm['pct'].iat[0]) + '_new.fasta')

myfinalasmfastasr2 = list()
for i in range(1,21):
    myfinalasm = rank2.head(i).tail(1)
    myfinalasmfastasr2.append(myprefix + '_' + str(myfinalasm['mincontig'].iat[0]) + '_' + str(myfinalasm['ovlrange'].iat[0]) + '_' + str(myfinalasm['pct'].iat[0]) + '_new.fasta')
'''
