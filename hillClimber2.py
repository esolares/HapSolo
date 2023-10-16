import random
import glob, os
import pandas as pd
import gzip as gz

mydir = './mosquito/'
l = glob.glob(mydir + '*.hap.gz')
myfn = glob.glob(mydir + '*.hap.gz')[0]
mydf = pd.read_csv(myfn, sep='\t', header=None, names=['qName', 'tName', 'qSize', 'QPct', 'PID', 'QRAlignLenPct'], dtype={'qName': object, 'tName': object})
os.chdir('/Users/mansiagrawal/Downloads/mosquito')
glob.glob('*')
myfn_1 = 'anofun_buscooutput.tsv.gz'

weight_missing = 1.0
weight_duplicate = 1.0
weight_fragmented = 0.0
weight_complete = 1.0
k_best_sol = 10
myqrycontigs = set(mydf['qName'])

contigsDictionary = dict()
buscosDictionary = dict() #key: busco, value = [count of complete, count of fragmented], each busco can only be added to one contig


class PriorityQueue(object):
    def __init__(self): #set of cost, contigs, pid, qpct
        self.pq = []
        self.max_size = k_best_sol

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
for line in gz.open(myfn_1):
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
    if (weight_complete*find_count_dict["complete"]) != 0:
        cost = cost/(weight_complete*find_count_dict["complete"])
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
    WriteNewAssembly("/Users/mansiagrawal/PycharmProjects/hillClimbing/primary_new.fasta","/Users/mansiagrawal/PycharmProjects/hillClimbing/primary.fasta",myGoodContigsSet)

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

hill_climbing2(mydf,100)
