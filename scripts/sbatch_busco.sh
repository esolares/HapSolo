#!/bin/bash
#SBATCH -J buscocha		# jobname
#SBATCH -o buscocha.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e buscocha.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-3			# start and stop of the array start-end
###SBATCH -n 1			# -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 16			# -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1priority,gcpriority			# queue (partition) -- normal, development, largemem, etc.
#SBATCH --mem-per-cpu=1750
###SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes

#please load your BUSCO, Augustus, BRAKER binaries here
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

mkdir -p busco
[ -d busco/busco${SLURM_ARRAY_TASK_ID} ] && rm -rf busco/busco${SLURM_ARRAY_TASK_ID}
mkdir -p busco/busco${SLURM_ARRAY_TASK_ID}
TMPDIR="./busco/busco${SLURM_ARRAY_TASK_ID}"
CWDIR=$(pwd)

SEED=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
cd ${TMPDIR}

echo "Begin analysis on ${SEED}"
#removes escape chars and spaces. bug fix for mummer. mummer will not take escape characters and spaces in fasta headers
#echo "Begin removing invalid characters in header on ${SEED}"
ln -sf ../../${SEED}
#cat ${SEED} | sed -r 's/[/ =,\t|]+/_/g' | awk -F "_" '{ if (/^>/) {printf($1"_"$2"\n")} else {print $0} }' > $(basename ${SEED} .fasta)_new.fasta
#QRY=$(basename ${SEED} .fasta)_new.fasta
QRY=${SEED}

echo "Begin quast analysis on ${QRY}"
quastrun="quast.py -t ${SLURM_CPUS_PER_TASK} ${QRY} -o quast_$(basename ${QRY} .fasta)"
echo $quastrun
$quastrun
echo "End quast analysis, cat results and begin busco run"
cat quast_$(basename ${QRY} .fasta)/report.txt > ${CWDIR}/$(basename ${QRY} .fasta)_scoresreport.txt
buscorun="BUSCO.py -c ${SLURM_CPUS_PER_TASK} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} .fasta)_${MYLIB}_${SPTAG} ${OPTIONS} -t ./run_$(basename ${QRY} .fasta)_${MYLIB}_${SPTAG}/tmp"
echo $buscorun
$buscorun
echo "End busco run and cat results"
cat run_$(basename ${QRY} .fasta)_${MYLIB}_${SPTAG}/short*.txt >> ${CWDIR}/$(basename ${QRY} .fasta)_scoresreport.txt
cd ..
#tar czf busco${SLURM_ARRAY_TASK_ID}.tar.gz busco${SLURM_ARRAY_TASK_ID}
#rm -rf busco${SLURM_ARRAY_TASK_ID}
echo "Finished on ${QRY}"

