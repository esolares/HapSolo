# HapSolo
An optimization approach for removing secondary haplotigs during diploid genome assembly and scaffolding.

# Installation requirements

HapSolo is compatible with Python 2.7 and requires the PANDAS package be installed.

To do this run:
```
pip install pandas
```
###### Note: HapSolo has been tested and run on BUSCO V3 and Blat V36. minimap2 support to come soon.
# Installation
```
https://github.com/esolares/HapSolo.git
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

We highly recommend you run preprocessfasta.py on your contig assembly first. The run syntax of this python script is:

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

HapSolo requires three arguments: Your preprocessed contig assembly file, your Blat PSL file (gzipped or uncompressed), and the location of your BUSCO results for each contig fasta file. The run syntax is as follows:
```
hapsolo.py -i YOURCONTIGASSEMBLY.fasta -p YOURPSLFILE.psl -b YOURBUSCOOUTPUTDIRECTORY
hapsolo.py -i contigassembly.fasta -p self_alignment.psl -b ./contigs/busco/

usage: hapsolo.py [-h] -i INPUT -p PSL -b BUSCOS [-m MAXZEROS] [-t THREADS]
                  [-n NITERATIONS] [-B BESTN] [-S THETAS] [-D THETAD]
                  [-F THETAF] [-M THETAM]

Process alignments and BUSCO"s for selecting reduced assembly candidates

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Fasta file
  -p PSL, --psl PSL     BLAT PSL alignnment file
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

```

# BUSCO run for each contig fasta file
To do this we recommend you run the following:
```
cd contigs
bash_buscopreprocess.sh
wc -l buscojobfile.txt
```
We will use the output of the last line and enter it into the job array of either the sbatch or qsub script. This should look like this: 1-N, where N is the number returned by wc -l. You will also need to set the number of cores to run for each array job. 
###### Note: For more information on how to run and configure array jobs, please contact your IT help desk or refer the your job scheduler's manual.
Next you wil need to make sure you have the proper lineage and species selected for your particular sample. You can look this up at https://busco.ezlab.org/ and http://bioinf.uni-greifswald.de/augustus/ respectively. Once these changes have been made to your sbatch_busco.sh or qsub_busco.sh script, you can submit them to your HPC queue.

