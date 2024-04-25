import start

def generated_df(mydf, step1, step2,step3,myqthreshold,mypidthreshold,qr_align):
    myfiltereddf = mydf[mydf['QPct'] >= myqthreshold + step1]
    mynewfiltereddf = myfiltereddf[myfiltereddf['PID'] >= mypidthreshold + step2]
    my_final_filtered_df = mynewfiltereddf[mynewfiltereddf['QRAlignLenPct'] >= qr_align + step3]
    return (my_final_filtered_df, step1, step2,step3)