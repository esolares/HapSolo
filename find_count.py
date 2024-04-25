import start
def find_count(mydf, step1, step2,step3, myqthreshold,mypidthreshold, qr_align):
    missing = 0
    duplicate = 0
    complete = 0
    fragmented = 0

    myfiltereddf = mydf[mydf['QPct'] >= myqthreshold + step1]
    mynewfiltereddf = myfiltereddf[myfiltereddf['PID'] >= mypidthreshold + step2]
    qr_align_filtered = mynewfiltereddf[mynewfiltereddf['QRAlignLenPct'] >= qr_align + step3]
    myfiltqrycontigs = set(qr_align_filtered['qName'].to_pandas())
    mynewset = start.myAllContigsSet - myfiltqrycontigs
    for contig in mynewset: 
        for busco in start.contigsDictionary[contig]: 
            for contig_and_type in start.BUSCOS2CTGSDICT[busco]: 
                if contig_and_type[1] == 'C':
                    complete += 1
                elif contig_and_type[1] == 'D':
                    duplicate += 1
                elif contig_and_type[1] == 'M':
                    missing+=1
                else:
                    fragmented+=1
    return {"fragmented": fragmented, "missing": missing, "complete": complete, "duplicate": duplicate}
