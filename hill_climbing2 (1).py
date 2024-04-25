import PriorityQueue
import find_count
import cost
import generated_df
import find_neighbors_cost
import random_restart
import start


def hill_climbing_2(mytuple):
    myqthreshold = mytuple[0]
    mypidthreshold = mytuple[1]
    qr_align = mytuple[2]
    priority_queue = PriorityQueue.PriorityQueue()
    step = 0.005
    searching_plateau = False
    current_cost = cost.cost(find_count.find_count(start.mydf,step,step,step,myqthreshold,mypidthreshold, qr_align))
    my_q_step = 0
    my_pid_step = 0
    qr_align_step = 0
    steps_searched_plateau = 0
    max_steps_in_plateau =10
    
    for i in range(1,start.niterations):
        step += 0.0005
        current_cost, my_q_step,my_pid_step, qr_align_step = find_neighbors_cost.find_neighbors_cost(step,myqthreshold,mypidthreshold, qr_align)
        myqthreshold +=my_q_step
        mypidthreshold +=my_pid_step
        qr_align+=qr_align_step
        if start.mydf.empty:
            myqthreshold, mypidthreshold,qr_align  = random_restart.random_restart(myqthreshold,mypidthreshold, qr_align)
            current_cost = cost.cost(find_count.find_count(start.mydf, step,step,step,myqthreshold,mypidthreshold, qr_align))
            continue
        if myqthreshold < 0 or mypidthreshold<0:
            myqthreshold, mypidthreshold,qr_align  = random_restart.random_restart(myqthreshold, mypidthreshold, qr_align)
            current_cost = cost.cost(find_count.find_count(start.mydf, step, step, step, myqthreshold, mypidthreshold, qr_align))
            continue

        if priority_queue.length() != 0 and (current_cost == priority_queue.head()):
            searching_plateau = True

        if priority_queue.length() == 0:
            filtered_df = start.mydf[start.mydf['QPct'] >= myqthreshold]
            filtered_df = filtered_df[filtered_df['PID'] >= mypidthreshold]
            filtered_df = filtered_df[filtered_df['QRAlignLenPct'] >= qr_align]
            priority_queue.insert_node([current_cost, filtered_df['qName'],myqthreshold,mypidthreshold, qr_align])
        if current_cost >  priority_queue.head():
            if searching_plateau == True:
                searching_plateau = False
                steps_searched_plateau = 0
            myqthreshold, mypidthreshold, qr_align  = random_restart.random_restart(myqthreshold,mypidthreshold, qr_align)
            current_cost = cost.cost(find_count.find_count(start.mydf, step,step,step,myqthreshold,mypidthreshold, qr_align))
            continue
        if current_cost < priority_queue.head():
            if searching_plateau == True:
                searching_plateau = False
                steps_searched_plateau = 0
            filtered_df = start.mydf[start.mydf['QPct'] >= myqthreshold]
            filtered_df = filtered_df[filtered_df['PID'] >= mypidthreshold]
            filtered_df = filtered_df[filtered_df['QRAlignLenPct'] >= qr_align]
            priority_queue.insert_node([current_cost,filtered_df['qName'],myqthreshold,mypidthreshold, qr_align])
        if (searching_plateau == True):
            steps_searched_plateau+=1
        if (searching_plateau == True and steps_searched_plateau == max_steps_in_plateau):
            searching_plateau = False
            steps_searched_plateau = 0
            myqthreshold, mypidthreshold, qr_align = random_restart.random_restart(myqthreshold, mypidthreshold, qr_align)
            current_cost = cost.cost(find_count.find_count(start.mydf, step, step, step,myqthreshold, mypidthreshold, qr_align))
            continue

    return priority_queue
