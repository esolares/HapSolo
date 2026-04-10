#!/bin/bash

#run example: bash bash_quastbusco.sh MYFASTA.FASTA JOBID

#please load your Quast, BUSCO, Augustus, BRAKER binaries here
source /networkshare/.mybashrc
export AUGUSTUS_CONFIG_PATH="/networkshare/bin/augustus-3.2.2/config"

#please copy to contigs directory and run in that directory
#please run bash_buscopreprocess.sh prior to running this

INPUTTYPE="geno"
#please enter the directory for your ODB9 libraries here
MYLIBDIR="/networkshare/bin/busco/lineages/"

#drosophila
MYLIB="diptera_odb9"
SPTAG="fly"

#vitis vinifera & vignis radiata
#MYLIB="embryophyta_odb9"
#SPTAG="arabidopsis"

#zea mayz
#MYLIB="embryophyta_odb9"
#SPTAG="mayz"

#heliconius
#MYLIB="insecta_odb9"
#SPTAG="heliconius_melpomene1"

#bony fish
#MYLIB="actinopterygii_odb9"
#SPTAG="zebrafish"

OPTIONS="-l ${MYLIBDIR}${MYLIB} -sp ${SPTAG}"
JOBFILE="buscojobfile.txt"

SEED=$1
JOBID=$2
CORES=$3

TMPDIR="./busco/busco${JOBID}"
CWDIR=$(pwd)

[ -d busco/busco${JOBID} ] && rm -rf busco/busco${JOBID}
mkdir -p busco/busco${JOBID}
TMPDIR="./busco/busco${JOBID}"
CWDIR=$(pwd)

SEED=$(head -n ${JOBID} ${JOBFILE} | tail -n 1)
cd ${TMPDIR}

echo "Begin analysis on ${SEED}"
#removes escape chars and spaces. bug fix for mummer. mummer will not take escape characters and spaces in fasta headers
#echo "Begin removing invalid characters in header on ${SEED}"
ln -sf ../../${SEED}
QRY=${SEED}

echo "Begin quast analysis on ${QRY}"
quastrun="quast.py -t ${CORES} ${QRY} -o quast_$(basename ${QRY} .fasta)"
echo $quastrun
$quastrun
echo "End quast analysis, cat results and begin busco run"
cat quast_$(basename ${QRY} .fasta)/report.txt > ${CWDIR}/$(basename ${QRY} .fasta)_scoresreport.txt
buscorun="BUSCO.py -c ${CORES} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} .fasta)_${MYLIB}_${SPTAG} ${OPTIONS} -t ./run_$(basename ${QRY} .fasta)_${MYLIB}_${SPTAG}/tmp"
echo $buscorun
$buscorun
echo "End busco run and cat results"
cat run_$(basename ${QRY} .fasta)_${MYLIB}_${SPTAG}/short*.txt >> ${CWDIR}/$(basename ${QRY} .fasta)_scoresreport.txt
cd ..
#tar czf busco${JOBID}.tar.gz busco${JOBID}
#rm -rf busco${JOBID}
echo "Finished on ${QRY}"
