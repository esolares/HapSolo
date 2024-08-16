<img src=./hapsolo_logo.png width=50%/>

An optimization approach for removing secondary haplotigs during diploid genome assembly and scaffolding.

# Singularity Installation Requirements
# RECOMMENDED
Our singularity image of HapSolo contains all required software including BUSCO, and Augustus. Please skip the local install if you are able to pull the Singularity image.
Running singularity will not require running conda and is the best way to get HapSolo and BUSCO3 working.
### Downloading the image (only do once)
```
# Pull with Singularity
singularity pull --arch amd64 library://esolares/default/hapsolo_busco3:0.01
# Pull by unique ID (reproducible even if tags change)
singularity pull library://esolares/default/hapsolo_busco3:sha256.5070a5b9119b11d7d50da0693c0c5185051a0e690e868a27367f17364ddb31be
```

### Singularity RUN COMMANDS:
This section will be required for running HapSolo and BUSCO3 using singularity initiate singularity and mount current working directory
<br>Augustus 3.2.2 Config directory Link: [config.tar.gz](https://drive.google.com/file/d/1zrHS6mqmIM5f_7n5DW4LrPd1sw32uUrf/view?usp=sharing)
<br>Note: You will need to download this then upload to your work directory saved in HOSTDIR. to see this echo $HOSTDIR
```

## MOUNDIR is the mounting point within singulariy
MOUNTDIR=/data
## HOSTDIR is your current working directory with all your files
## Note: You will need the Augustus config folder in your work directory as well as the contigs folder and the odb9 folder for your taxa
HOSTDIR=$(pwd)

# You will need to download the Augustus 3.2.2 config directory. I have provided a link from gdrive above.
# Note: You will need to download this then upload to your work directory saved in HOSTDIR. to see this echo $HOSTDIR

# I will work on getting this so you can execute curl or wget. For now you will need to download using a browser and upload via scp or winscp or other method
# once it is upload you will need to run: 
tar xzvf config.tar.gz

# Here we will spin up the singularity instance. Please make sure you are running all these commands in your compute node
# Note: Spinning up the instance will keep this loaded on the server until you execute the singularity stop instance://hapsolo command
singularity instance start --bind ${HOSTDIR}:${MOUNTDIR} hapsolo_busco3_0.01.sif hapsolo

# Here we will test BUSCO3 and HapSolo as well as ls your current working directory
singularity exec instance://hapsolo BUSCO.py --help
singularity exec instance://hapsolo hapsolo.py --help
singularity exec instance://hapsolo ls /data/
# you should now see all your files. This means that singularity has access to all your files from the instance. If it does not match please check your HOSTDIR 
# and make sure you are mounting your HOSTDIR to your MOUNTDIR as noted above.
# CONGRATULATIONS! You know have singularity running.

# Here is an example BUSCO3 command
singularity exec instance://hapsolo BUSCO.py -c 24 -m geno -i /data/contigs/tig00000002.fasta -l /data/embryophyta_odb9 -sp arabidopsis -t ./tmp -o testbusco

# Please make sure to append "singularity exec instance://hapsolo" to any command you see in this documentation or scripts to get them to work with singularity. 
# Also please make sure you have hapsolo running as an instance. You can check by running the following:
singularity instance list

# To stop your instance you can run the following
singularity instance stop hapsolo

# Please submit bug reports if you have any issues. Sometimes google will not notify me and you can reach me via my website: edwinsolares.com or on twitter: @edwinsolares10
# If you are having errors with your odb9, make sure you are not linking it (ln -s) and have a hard copy in your working directory i.e. ${HOSTDIR}/embryophyta_odb9
```

# Local Installation Requirements (Local install not required if you are running singularity)
## Please try to use the singularity image instead as there have been issues with the conda version of BUSCO.
HapSolo is compatible with Python 2.7 and requires the PANDAS package be installed. There is support for Python 3, but Python 2.7 runs faster.

To do this please install conda and run:
```
conda create --name HapSolo python=2.7
conda activate HapSolo
conda install -c anaconda pandas
```
###### Note: HapSolo has been tested and run on BUSCO V3 and Blat V36 and minimap2. We do not recommend using newer versions of BUSCO as it will lead to miss-classification of BUSCO's.
# Installation
```
git clone https://github.com/esolares/HapSolo.git
cd HapSolo
HAPSOLO=$(pwd)
```
Make sure to export your path to your ENV variables by doing the following and adding this line to all your scripts that run preprocessfasta.py or hapsolo.py
```
export PATH=$HAPSOLO:$PATH
export PATH=$HAPSOLO/scripts:$PATH
```
###### Note: If you are running a custom environment for python, make sure to replace the first line of preprocessfasta.py and hapsolo.py to reflect the proper python running environment

# How to run HapSolo
HapSolo requires a Blat alignment file and a busco directory that contains a busco output for each of the contigs in your contig assembly. To do this we have included a scripts directory that contains scripts that can be used for preprocessing and postprocessing for your HapSolo run. 
###### Note: For simplicity job submission scripts have been added: sbatch_preprocess.sh and qsub_preprocess.sh for SLURM and SGE respectively that only require that the REF variable be assigned to the name of your contig assembly FASTA file.

We highly recommend you run preprocessfasta.py on your contig assembly first. 

The run syntax of this python script is:

```
preprocessfasta.py -i CONTIGASSEMBLY.fasta -m INTEGERSIZEOFMAXCONTIGFORQUERY

usage: preprocessfasta.py [-h] -i INPUT [-m MAXCONTIG]

Preprocess FASTA file and outputs a clean FASTA and seperates contigs based on
unique headers. Removes special chars

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file
  -m MAXCONTIG, --maxcontig MAXCONTIG
                        Max size of contig in Mb to output in contigs folder.
```

By default the maxcontig size is set to 10000000 or 10Mb and is a not a required parameter.

This script will create a contigs directory in your current working directory that contains individual FASTA files for each contig. This script also removes any illegal characters in the FASTA header that could cause inconsistent results from either BUSCO or MUMmer. A qsub and sbatch file have been included in the scripts folder that also create the required jobfile.txt and buscojobfile.txt files, as well as convert your fasta files to 2bit for Blat alignment.

HapSolo requires three arguments: Your preprocessed contig assembly file, your Blat PSL file (gzipped or uncompressed), and the location of your BUSCO results for each contig fasta file. 

The run syntax is as follows:
```
hapsolo.py -i YOURPREPROCESSEDCONTIGASSEMBLY.fasta --psl YOURPSLFILE.psl -b YOURBUSCOOUTPUTDIRECTORY
or
hapsolo.py -i YOURPREPROCESSEDCONTIGASSEMBLY.fasta --paf YOURPAFFILE.paf -b YOURBUSCOOUTPUTDIRECTORY
hapsolo.py -i contigassembly_new.fasta --paf/psl self_alignment.file -b ./contigs/busco/

usage: hapsolo.py [-h] -i INPUT (-p PSL | -a PAF) -b BUSCOS [-m MAXZEROS] [-t THREADS]
                  [-n NITERATIONS] [-B BESTN] [-S THETAS] [-D THETAD]
                  [-F THETAF] [-M THETAM]

Process alignments and BUSCO"s for selecting reduced assembly candidates

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Fasta file
  -p PSL, --psl PSL     BLAT PSL alignnment file
  -a PAF, --paf PAF     Minimap2 PAF alignnment file. Note. paf file
                        functionality is currently experimental

  -b BUSCOS, --buscos BUSCOS
                        Location BUSCO output directories. i.e. buscoN/
  -m MAXZEROS, --maxzeros MAXZEROS
                        Max # of times cost function delta can consecutively
                        be 0. Default = 10
  -t THREADS, --threads THREADS
                        # of threads. Multiplies iterations by threads.
                        Default = 1
  -n NITERATIONS, --niterations NITERATIONS
                        # of total iterations to run per gradient descent.
                        Default = 1000
  -B BESTN, --Bestn BESTN
                        # of best candidate assemblies to return using
                        gradient descent. Default = 1
  -S THETAS, --thetaS THETAS
                        Weight for single BUSCOs in linear fxn. Default = 1.0
  -D THETAD, --thetaD THETAD
                        Weight for single BUSCOs in linear fxn. Default = 1.0
  -F THETAF, --thetaF THETAF
                        Weight for single BUSCOs in linear fxn. Default = 0.0
  -M THETAM, --thetaM THETAM
                        Weight for single BUSCOs in linear fxn. Default = 1.0

-p/--psl and -a/--paf are mutually exclusive

```
# All-by-All Alignment. Please choose Blat or MiniMap2

# MiniMap2
We recommend using one of the following options below when running MiniMap2:
```
minimap2 -t 36 -k19 -w5 -A1 -B2 -O3,13 -E2,1 -s200 -z200 -N50 --min-occ-floor=100 ${QRY} ${QRY} > $(basename ${QRY} .fasta)_self_align.paf
minimap2 -t 36 -P -k19 -w2 -A1 -B2 -O1,6 -E2,1 -s200 -z200 -N50 --min-occ-floor=100 ${QRY} ${QRY} > $(basename ${QRY} .fasta)_self_align.paf
minimap2 -t 36 -P -G 500k -k19 -w2 -A1 -B2 -O2,4 -E2,1 -s200 -z200 -N50 --max-qlen 10000000 --min-occ-floor=100 --paf-no-hit ${QRY} ${QRY} > $(basename ${QRY} .fasta)_self_align.paf
```
-t <INT> can be set to the number of cores allocated for the job. In SLURM the variable is SLURM_CPUS_PER_TASK and in SGE the variable is NSLOTS.
We also recommend testing different options to see if you get better results. HapSolo also runs faster using MiniMap2 paf files.

# Blat all-by-all alignment 
Here we provide scripts for running Blat on an HPC using SGE or SLURM job schedulers.
To preprocess, you should have run the qsub_preprocess.sh or the sbatch_preprocess.sh scripts above. These files only require that you assign your FASTA file name to the REF variable. Once the preprocessing step has been completed a jobfile.txt file, a new FASTA file containing _new.fasta appended to the name and a contigs directory containing FASTA files for each contig will have been created. This step does require quite a bit of memory, so it is recommended that you consider this prior to running this step.

To submit the Blat all-by-all alignment job, we have provided scripts for either SLURM or SGE, sbatch_blat.sh and qsub_blat.sh.

To run:
```
sbatch sbatch_preprocess.sh
wc -l jobfile.txt
```
Here we recommend you insert the output of wc -l to the array setting in your job submission script. Here you will need to assign the new FASTA file that has been created by the preprocess.py Python script to the REF variable in the blat.sh submission scripts. For all intensive purposes, we recommend that you use this new file for all subsequent runs.

To run:
```
sbatch sbatch_blat.sh
```
###### Note: If the RAM usage is managed by your job scheduler, you will need to alot a sufficient amount of RAM to your job submission script. The contig size will determine how much RAM will be required. We recommend setting this to a minimum of 4GB. Larger files may require 8-16GB of RAM. For mammals, which are repeat rich, require a minimum value of 16GB. 
If all alignment jobs completed successfully, we will then need to concatenate the individual PSL files into a larger aligment file. To do this we have provided the bash_andreaconcatpsl.sh script provided by A. Minio at the Cantu Lab in UC Davis.

To run:
```
bash_andreaconcatpsl.sh myoutput_selfaln.PSL
```
This script will look for all job*.PSL files that have been created by the parallel all-by-all alignment. 
###### Note: You could also use the parallel version of Blat, pblat available here https://github.com/icebert/pblat

# BUSCO run for each contig fasta file
Here we need to extract contigs larger than 10Mb. Make sure to run this after the initial preprocess.py step mentioned above.
To do this we recommend you run the following:
```
bash_buscopreprocess.sh
cd contigs
wc -l buscojobfile.txt
```
We will use the output of the last line and enter it into the job array of either the sbatch or qsub script. This should look like this: 1-N, where N is the number returned by wc -l. You will also need to set the number of cores to run for each array job. 
###### Note: For more information on how to run and configure array jobs, please contact your IT help desk or refer the your job scheduler's manual.
Next you wil need to make sure you have the proper lineage and species selected for your particular sample. You can look this up at https://busco.ezlab.org/ and http://bioinf.uni-greifswald.de/augustus/ respectively. Once these changes have been made to your sbatch_busco.sh or qsub_busco.sh script, you can submit them to your HPC queue.


# run HapSolo
Once you have concatenated all your psl output files or finished minimap2, and finished your BUSCO3 run. You can run HapSolo. Currently HapSolo does not accept absolute paths. This is a known bug and wil be fixed shortly. If your current working directy is:
```
/nfsmount/mydir/hapsolowd
```
and your input files shoud look like this:
```
/nfsmount/mydir/hapsolowd/myalignmentfile
/nfsmount/mydir/hapsolowd/myassembly.fasta
/nfsmount/mydir/hapsolowd/contigs/busco
```
then you should run HapSolo as the following:
```
hapsolo.py -i myassembly.fasta --psl myalignmentfile.psl -b contigs/busco
```
or
```
hapsolo.py -i myassembly.fasta --paf myalignmentfile.paf -b contigs/busco
```
