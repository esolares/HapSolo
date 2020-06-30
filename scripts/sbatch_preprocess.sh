#!/bin/bash
#SBATCH -J aeprep		# jobname
#SBATCH -o aeprep.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e aeprep.e%A.%a	# error file name A is the jobid and a is the arraytaskid
###SBATCH -a 1-1		# start and stop of the array start-end
###SBATCH -n 1			# -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 1			# -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1,gcluster			# queue (partition) -- normal, development, largemem, etc.
#SBATCH --mem-per-cpu=1750
###SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes

# This is the variable that has the total number of cores committed to this job. In this case it is 1
#echo $SLURM_CPUS_PER_TASK

module load blat

REF=myassemblyfile.fasta

#create this file: faToTwoBit ${REF} $(basename ${REF} .fasta).2bit
faToTwoBit ${REF} $(basename ${REF} .fasta).2bit

python preprocessfasta.py -i ${REF}

JOBFILE=jobfile.txt
#ls contigs/*.fasta > ${JOBFILE}
wait
