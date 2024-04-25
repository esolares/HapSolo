import cost
import find_count
import generated_df
import start
def find_neighbors_cost(step,myqthreshold,mypidthreshold, qr_align):
    cost1 = cost.cost(find_count.find_count(start.mydf, step, step, -step,  myqthreshold, mypidthreshold, qr_align))
    min_cost = (cost1,generated_df.generated_df(start.mydf, step, step,-step, myqthreshold,mypidthreshold, qr_align)[1],generated_df.generated_df(start.mydf, step, step, -step,myqthreshold,mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, step, step,-step, myqthreshold,mypidthreshold, qr_align)[3])

    cost2 = cost.cost(find_count.find_count(start.mydf, step, step, step, myqthreshold, mypidthreshold, qr_align))  # 1st quadrant
    min_cost = (cost1, generated_df.generated_df(start.mydf, step, step, step, myqthreshold, mypidthreshold, qr_align)[1], generated_df.generated_df(start.mydf, step, step, step, myqthreshold, mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, step, step,step, myqthreshold,mypidthreshold, qr_align)[3])


    cost3= cost.cost(find_count.find_count(start.mydf, -step, step, step, myqthreshold,mypidthreshold, qr_align))
    if (cost3<min_cost[0]):
        min_cost = (cost3, generated_df.generated_df(start.mydf, -step, step,step, myqthreshold,mypidthreshold, qr_align)[1],generated_df.generated_df(start.mydf, -step, step,step, myqthreshold,mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, -step, step,step, myqthreshold,mypidthreshold, qr_align)[3])

    cost4= cost.cost(find_count.find_count(start.mydf, -step, step, -step, myqthreshold,mypidthreshold, qr_align))
    if (cost4<min_cost[0]):
        min_cost = (cost4, generated_df.generated_df(start.mydf, -step, step,-step, myqthreshold,mypidthreshold, qr_align)[1],generated_df.generated_df(start.mydf, -step, step,-step, myqthreshold,mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, -step, step,-step, myqthreshold,mypidthreshold, qr_align)[3])

    cost5 = cost.cost(find_count.find_count(start.mydf, -step, -step,step,myqthreshold,mypidthreshold, qr_align))
    if (cost5<min_cost[0]):
        min_cost = (cost5, generated_df.generated_df(start.mydf, -step, -step,step,myqthreshold,mypidthreshold,qr_align)[1],generated_df.generated_df(start.mydf, -step, -step,step,myqthreshold,mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, -step, -step,step, myqthreshold,mypidthreshold, qr_align)[3])

    cost6 = cost.cost(find_count.find_count(start.mydf, -step, -step, -step, myqthreshold, mypidthreshold, qr_align))
    if (cost6<min_cost[0]):
        min_cost = (cost6, generated_df.generated_df(start.mydf, -step, -step,-step,myqthreshold,mypidthreshold,qr_align)[1],generated_df.generated_df(start.mydf, -step, -step,-step,myqthreshold,mypidthreshold, qr_align)[2], generated_df.generated_df(start.mydf, -step, -step,-step, myqthreshold,mypidthreshold, qr_align)[3])

    cost7 = cost.cost(find_count.find_count(start.mydf, step, -step,step,myqthreshold,mypidthreshold,qr_align))
    if (cost7<min_cost[0]):
        min_cost = (cost7, generated_df.generated_df(start.mydf, step, -step,step,myqthreshold,mypidthreshold,qr_align)[1],generated_df.generated_df(start.mydf, step, -step,step,myqthreshold,mypidthreshold,qr_align)[2], generated_df.generated_df(start.mydf, step, -step,step, myqthreshold,mypidthreshold, qr_align)[3])

    cost8 = cost.cost(find_count.find_count(start.mydf, step, -step,-step,myqthreshold,mypidthreshold,qr_align))
    if (cost8<min_cost[0]):
        min_cost = (cost8, generated_df.generated_df(start.mydf, step, -step,-step,myqthreshold,mypidthreshold,qr_align)[1],generated_df.generated_df(start.mydf, step, -step,-step,myqthreshold,mypidthreshold,qr_align)[2], generated_df.generated_df(start.mydf, step, -step,-step, myqthreshold,mypidthreshold, qr_align)[3])
    return min_cost