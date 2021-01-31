#!/bin/bash
# Please also cite GNU Parallel if you use this script

module load blat

### Set number of cores available here
CORES=48

REF=myfastafile.fasta
### create this file in preprocessing step: faToTwoBit ${REF} $(basename ${REF} .fasta).2bit
REF=$(basename ${REF} .fasta).2bit

#JOBFILE contains a list of contig fasta files
JOBFILE=jobfile.txt

ls *.fa | sed 's/.fa//' > ${JOBFILE}

cat ${JOBFILE} | parallel -j $CORES echo {}.fa {}.2bit

cat ${JOBFILE} | parallel -j $CORES echo $REF {}.2bit {}.psl
