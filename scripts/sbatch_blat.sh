#!/bin/bash
#SBATCH -J blataf		# jobname
#SBATCH -o blataf.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e blataf.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-1000                  # start and stop of the array start-end
###SBATCH -n 1                  # -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 1                    # -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1                   # queue (partition) -- normal, development, largemem, etc.
#SBATCH --mem-per-cpu=4000
#SBATCH -p p1,gcluster			# queue (partition) -- normal, development, largemem, etc.
###SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes

module load blat

REF=myfastafile.fasta

#create this file in preprocessing step: faToTwoBit ${REF} $(basename ${REF} .fasta).2bit
REF=$(basename ${REF} .fasta).2bit

#JOBFILE contains a list of contig fasta files
JOBFILE=jobfile.txt

SEED=$(head -n $SLURM_ARRAY_TASK_ID $JOBFILE | tail -n 1)
OUTPUT=job${SLURM_ARRAY_TASK_ID}_$(basename ${REF} .fasta)_self_aln.psl
faToTwoBit ${SEED} $(basename ${SEED} .fasta).2bit
QRY=$(basename ${SEED} .fasta).2bit
echo ${SEED}
echo "begin blat ${REF} ${QRY} ${OUTPUT}"
blat ${REF} ${QRY} ${OUTPUT}
echo "job ${SGE_TASK_ID} on ${QRY} alignment finished"
