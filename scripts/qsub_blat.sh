#!/bin/bash
#$ -N myblatrun
#$ -t 1-1000
#$ -pe openmp 2
#$ -R Y
#$ -ckpt blcr
#$ -q bigmemory

module load blat
#PATH=/networkshare/bin/blat:$PATH

REF=primary_new.fasta

#create this file: faToTwoBit ${REF} $(basename ${REF} .fasta).2bit
#REF=$(basename ${REF} .fasta).2bit

#JOBFILE contains list of contig fasta files
JOBFILE=jobfile.txt
#ls contigs/*.fa > ${JOBFILE}
SEED=$(head -n ${SGE_TASK_ID} ${JOBFILE} | tail -n 1)
OUTPUT=job${SGE_TASK_ID}_$(basename ${REF} .fasta)_self_aln.psl
faToTwoBit ${SEED} $(basename ${SEED} .fasta).2bit
QRY=$(basename ${SEED} .fasta).2bit
echo ${SEED}
echo "begin blat ${REF} ${QRY} ${OUTPUT}"
blat ${REF} ${QRY} ${OUTPUT}
echo "job ${SGE_TASK_ID} on ${QRY} alignment finished"
