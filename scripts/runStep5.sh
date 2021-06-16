#!/bin/bash

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 1-22

## code to run step 5 run it by qsub runStep5.sh
## to allow for superpopulation specific runs we add
## a pop argument

pop=$1

source env/bin/activate

python3 P03_extract_bpoints.py \
	--name=chr"$SGE_TASK_ID" \
	--subset=fourier_ls  \
	--dataset_path="$pop"/ \
	--input_pickle_fname="$pop"/chr"$SGE_TASK_ID"_minima.pickle \
	> "$pop"/chr"$SGE_TASK_ID".bed
