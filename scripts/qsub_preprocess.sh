#!/bin/bash
#$ -N mypreprocessrun
#$ -t 1-1
#$ -pe openmp 2
#$ -R Y
#$ -ckpt blcr
#$ -q myqueue

module load blat

REF=myassemblyfile.fasta

#create this file: faToTwoBit ${REF} $(basename ${REF} .fasta).2bit
faToTwoBit ${REF} $(basename ${REF} .fasta).2bit

python preprocessfasta.py -i ${REF}

JOBFILE=jobfile.txt
ls contigs/*.fasta > ${JOBFILE}
wait
