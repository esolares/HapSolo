#!/usr/bin/env python3
import random, glob, os, argparse
import pandas as pd
import gzip as gz

# usage haplotigreduction.py mypslfile.psl myfastafile.fasta buscoresults.tsv
parser = argparse.ArgumentParser(description='Process alignments and BUSCO"s for selecting reduced assembly candidates', epilog='-p/--psl, -a/--paf and -H/--hap are mutually exclusive')
parser.add_argument('-i', '--input', help='Input Fasta file', type=str, required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
group.add_argument('-a', '--paf', help='Minimap2 PAF alignnment file. Note. paf file functionality is currently experimental', type=str)
group.add_argument('-H', '--hap', help='Import .hap summary alignment file instead', type=str)

#mode = parser.add_mutually_exclusive_group(required=True)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
#mode.add_argument('-p', '--psl', help='BLAT PSL alignnment file', type=str)
parser.add_argument('--mode', help='HapSolo run mode. 0 = Random walking, 1 = No optimization with defaults, 2 = Optimized walking, Default = 0', type=int, required=False)

parser.add_argument('-b', '--buscos', help='Location BUSCO output directories. i.e. buscoN/', type=str, required=True)
parser.add_argument('-m', '--maxzeros', help='Max # of times cost function delta can consecutively be 0. Default = 10', type=str, required=False)
parser.add_argument('-t', '--threads', help='# of threads. Multiplies iterations by threads. Default = 1', type=int, required=False)
#parser.add_argument('-t', '--threads', help='# of threads. Multiplies iterations by threads. Default = 1', type=int, required=False)
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
hapalignmentfile = args.hap
buscofileloc = args.buscos
maxzeros = args.maxzeros
threads = args.threads
iterations = args.niterations
n_best_sol = args.Bestn
weight_single = args.thetaS
weight_duplicate = args.thetaD
weight_missing = args.thetaM
weight_fragmented = args.thetaF
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
if n_best_sol == None:
    n_best_sol = 1
if weight_single == None:
    weight_single = 1.0
if weight_duplicate == None:
    weight_duplicate = 1.0
if weight_missing == None:
    weight_missing = 1.0
if weight_fragmented == None:
    weight_fragmented = 0.0
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
    print('-R/--minQR set to a value less than 0.02. using 0.02 instead')

if myMinContigSize == None or myMinContigSize < 0:
    myMinContigSize = 1000

#mydir = './mosquito/'
#l = glob.glob(mydir + '*.hap.gz')
#myfn = glob.glob(mydir + '*.hap.gz')[0]
#os.chdir('/Users/mansiagrawal/Downloads/mosquito')
#glob.glob('*')
#myfn_1 = 'anofun_buscooutput.tsv.gz'


#var intialization
try:
   mydf = pd.read_csv(hapalignmentfile, sep='\t', header=None, names=['qName', 'tName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct'], dtype={'qName': object, 'tName': object})
except:
   print('Error reading', myfn + '.', 'Please check file and location and try again')
   quit()

#weight_missing = 1.0
#weight_duplicate = 1.0
#weight_fragmented = 0.0
#weight_single = 1.0
#n_best_sol = 10
myqrycontigs = set(mydf['qName'])

contigsDictionary = dict()
buscosDictionary = dict() #key: busco, value = [count of complete, count of fragmented], each busco can only be added to one contig


class PriorityQueue(object):
    def __init__(self): #set of cost, contigs, pid, qpct
        self.pq = []
        self.max_size = n_best_sol

    def length(self):
        return len(self.pq)

    def delete_node(self):
        max_value = 0
        for i in range(len(self.pq)):
            if self.pq[i][0] > self.pq[max_value][0]:
                max_value = i
        item = self.pq[max_value]
        del self.pq[max_value]
        return item

    def print_queue(self):
        for item in self.pq:
            print(item)

    def add_set(self):
        myGoodSet = set()
        for item in self.pq:
            myGoodSet.update(set(item[1]))
        return myGoodSet


    def head(self):
        max_value = 0
        for i in range(len(self.pq)):
            if self.pq[i][0] > self.pq[max_value][0]:
                max_value = i
        item = self.pq[max_value]
        return item[0]


    def insert_node(self, data):
        if (self.length() < self.max_size):
            self.pq.append(data)
        else:
            max_value = 0
            for i in range(len(self.pq)):
                if self.pq[i][0] > self.pq[max_value][0]:
                    max_value = i
            if data[0] < self.pq[max_value][0]:
                self.delete_node()
                self.pq.append(data)

#disct with buscos as keys, dict with contigs as keys.
#1) static contigs to buscos, static buscos to contigs, global dict with contigs

for elem in myqrycontigs:
    contigsDictionary[elem] = set()
#busco, busco type ,contig
temp = set()
BUSCOS2CTGSDICT = dict() #populate contigs dict from allignment file in mydf
for line in gz.open(buscofileloc): #myfn_1):
    line = line.decode()
    line = line.strip().split()

    if line[1][0] != 'M':
        if line[2] not in contigsDictionary.keys(): #look up runtime of .keys()
            contigsDictionary[line[2]] = set() #use disctionary O(1)
            #print(line[2])
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


def find_count(mydf, step1, step2,myqthreshold,mypidthreshold):
    missing = 0
    duplicate = 0
    complete = 0
    fragmented = 0

    mytemp1 = mydf[mydf['QPct'] >= myqthreshold + step1]
    myfiltereddf = mytemp1[mytemp1['PID'] >= mypidthreshold + step2]
    myfiltqrycontigs = set(myfiltereddf['qName'])
    mynewset = myqrycontigs - myfiltqrycontigs
    for contig in mynewset: #lookup constant time (repopulate the contigsDictionary and BUSCOS2CTGSDICT) (try except for error checking if needed)
        for busco in contigsDictionary[contig]: #each busco needs a counter in it, need to make copy of busco dict
            for contig_and_type in BUSCOS2CTGSDICT[busco]: #need busco count dict, populate intially with empty lists (0 vals) and copy when we make step
                if contig_and_type[1] == 'C':
                    complete += 1
                elif contig_and_type[1] == 'D':
                    duplicate += 1
                elif contig_and_type[1] == 'M':
                    missing+=1
                else:
                    fragmented+=1
    return {"fragmented": fragmented, "missing": missing, "complete": complete, "duplicate": duplicate}


def cost(find_count_dict):
    cost = (weight_missing*find_count_dict["missing"])+(weight_duplicate*find_count_dict["duplicate"])+(find_count_dict["fragmented"]*weight_fragmented)
    if (weight_single*find_count_dict["complete"]) != 0:
        cost = cost/(weight_single*find_count_dict["complete"])
    return cost

def generated_df(mydf, step1, step2,myqthreshold,mypidthreshold):
    mytemp1 = mydf[mydf['QPct'] >= myqthreshold + step1]
    myfiltereddf = mytemp1[mytemp1['PID'] >= mypidthreshold + step2]
    return (myfiltereddf, step1, step2)

def find_neighbors_cost(step,myqthreshold,mypidthreshold):
    cost1 = cost(find_count(mydf, step, step,myqthreshold,mypidthreshold)) # 1st quadrant
    min_cost = (cost1,generated_df(mydf, step, step,myqthreshold,mypidthreshold)[1],generated_df(mydf, step, step,myqthreshold,mypidthreshold)[2])
    cost2= cost(find_count(mydf, -step, step,myqthreshold,mypidthreshold))  # 2nd quadrant
    if (cost2<min_cost[0]):
        min_cost = (cost2, generated_df(mydf, -step, step,myqthreshold,mypidthreshold)[1],generated_df(mydf, -step, step,myqthreshold,mypidthreshold)[2])
    cost3 = cost(find_count(mydf, -step, -step,myqthreshold,mypidthreshold))  # 3rd quadrant
    if (cost3<min_cost[0]):
        min_cost = (cost3, generated_df(mydf, -step, -step,myqthreshold,mypidthreshold)[1],generated_df(mydf, -step, -step,myqthreshold,mypidthreshold)[2])
    cost4 = cost(find_count(mydf, step, -step,myqthreshold,mypidthreshold))  # 4th quadrant
    if (cost4<min_cost[0]):
        min_cost = (cost4, generated_df(mydf, step, -step,myqthreshold,mypidthreshold)[1],generated_df(mydf, step, -step,myqthreshold,mypidthreshold)[2])
    return min_cost # (minimum cost, step 1, step 2)


def random_restart(myqthreshold,mypidthreshold):
    step1 = random.uniform(-myqthreshold, 1.0 -myqthreshold)
    step2 = random.uniform(-mypidthreshold, 1.0 - mypidthreshold)
    return (step1+myqthreshold, step2+mypidthreshold)



def hill_climbing2(mydf,niterations):
    priority_queue = PriorityQueue()
    myqthreshold = 0.7
    mypidthreshold = 0.7
    step = 0.005
    searching_plateau = False
    current_cost = cost(find_count(mydf,step,step,myqthreshold,mypidthreshold))
    my_q_step = 0
    my_pid_step = 0
    steps_searched_plateau = 0
    max_steps_in_plateau =10
    for i in range(1,niterations):
        step += 0.05
        current_cost, my_q_step,my_pid_step = find_neighbors_cost(step,myqthreshold,mypidthreshold)
        myqthreshold +=my_q_step
        mypidthreshold +=my_pid_step
        if mydf.empty:
            myqthreshold, mypidthreshold  = random_restart(myqthreshold,mypidthreshold)
            current_cost = cost(find_count(mydf, step,step,myqthreshold,mypidthreshold))
            continue
        if myqthreshold < 0 or mypidthreshold<0:
            myqthreshold, mypidthreshold  = random_restart(myqthreshold, mypidthreshold)
            current_cost = cost(find_count(mydf, step, step, myqthreshold, mypidthreshold))
            continue

        if priority_queue.length() != 0 and (current_cost == priority_queue.head()):
            searching_plateau = True

        if priority_queue.length() == 0:
            filtered_df = mydf[mydf['QPct'] >= myqthreshold]
            filtered_df = filtered_df[filtered_df['PID'] >= mypidthreshold]

            priority_queue.insert_node([current_cost, filtered_df['qName'],myqthreshold,mypidthreshold])
        if current_cost >  priority_queue.head():
            if searching_plateau == True:
                searching_plateau = False
                steps_searched_plateau = 0
            myqthreshold, mypidthreshold  = random_restart(myqthreshold,mypidthreshold)
            current_cost = cost(find_count(mydf, step,step,myqthreshold,mypidthreshold))
            continue
        if current_cost < priority_queue.head():
            if searching_plateau == True:
                searching_plateau = False
                steps_searched_plateau = 0
            filtered_df = mydf[mydf['QPct'] >= myqthreshold]
            filtered_df = filtered_df[filtered_df['PID'] >= mypidthreshold]
            priority_queue.insert_node([current_cost,filtered_df['qName'],myqthreshold,mypidthreshold])

        if (searching_plateau == True):
            steps_searched_plateau+=1
        if (searching_plateau == True and steps_searched_plateau == max_steps_in_plateau):
            searching_plateau = False
            steps_searched_plateau = 0
            myqthreshold, mypidthreshold = random_restart(myqthreshold, mypidthreshold)
            current_cost = cost(find_count(mydf, step, step, myqthreshold, mypidthreshold))
            continue
    myGoodContigsSet = priority_queue.add_set()
    WriteNewAssembly(myasmFileName,"./primary.hap.fasta",myGoodContigsSet)

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
    while count < totalLines:
        lastPos = headerPos = fin.tell()
        line = fin.readline().replace('\n', '')
        count = count + 1
        if line[0:1] == '>':
            seqName = line.split(" ")[0].replace('>', '').replace('/', '_')
            lastPos = startPos = fin.tell()
            line = fin.readline().replace('\n', '')
            count = count + 1
            while line[0:1] != '>' and line[0:1] != '':
                seqLen = seqLen + len(line)
                endPos = lastPos
                lastPos = fin.tell()
                line = fin.readline().replace('\n', '')
                count = count + 1
            if line[0:1] == '>' or line[0:1] == '':
                myContigSizeDict[seqName] = [seqLen, headerPos, startPos, endPos]
                seqName = ''
                seqLen = 0
                count = count - 1
                fin.seek(lastPos)
    fin.close()
    return myContigSizeDict

def WriteNewAssembly(myasmFileName, newASMFileName, myGoodContigsSet):
    fin = open(myasmFileName, 'r')
    fout = open(newASMFileName, 'w')
    myContigSizeDict = CalculateContigSizes(myasmFileName)
    for contig in myGoodContigsSet:
        myContigPositionsList = myContigSizeDict[contig]
        fin.seek(myContigPositionsList[1])  # extract headerpos
        fout.write(fin.readline())
        newPos = fin.tell()
        mySeq = fin.readline().replace('\n', '')
        while newPos != myContigPositionsList[3]:
            newPos = fin.tell()
            mySeq = mySeq + fin.readline().replace('\n', '')
        fout.write(mySeq + '\n')
    fout.close()

hill_climbing2(mydf, iterations)
