#!/bin/bash

REF=YOURASM.FASTA

# make sure to start from your root working directory that contains your assmebly
python preprocessfasta.py -i ${REF} -m 0

cd contigs
# make sure to run in ${WORKDIR}/contigs/ directory
ls *.fasta > buscojobfile.txt
