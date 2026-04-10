#!/bin/bash

### Set number of cores available here
CORES=48

JOBFILE="buscojobfile.txt"

mkdir -p busco

# make sure to make changes to bash_quastbusco.sh file for your binaries, and busco odb9 libraries
parallel -j ${CORES} -k bash bash_quastbusco.sh {} {#} $CORES ::: $(cat {${JOBFILE})
