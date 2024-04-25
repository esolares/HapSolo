import start
class PriorityQueue(object):
    def __init__(self): #set of cost, contigs, pid, qpct
        self.pq = []
        self.max_size = start.n_best_sol

    def length(self):
        return len(self.pq)

    def delete_node(self):
        max_value = 0
        for i in range(len(self.pq)):
            if self.pq[i][0] > self.pq[max_value][0]:
                max_value = i
        item = self.pq[max_value]
        self.pq.pop(max_value)
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
        min_value = 0
        for i in range(len(self.pq)):
            if self.pq[i][0] < self.pq[min_value][0]:
                min_value = i
        item = self.pq[min_value]
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

    def return_pq(self):
        return self.pq
