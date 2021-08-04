#!/bin/bash

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 1-22

## code to run step 4 run it by qsub runStep4.sh
## to allow for superpopulation specific analysis we
## added a pop subdir
pop=$1
nsnps=$2

source env/bin/activate

python3 P02_minima_pipeline.py \
	--input_fname="$pop"/chr"$SGE_TASK_ID".vector.txt.gz \
	--chr_name=chr"$SGE_TASK_ID" \
	--dataset_path="$pop"/ \
	--n_snps_bw_bpoints="$nsnps" \
	--out_fname="$pop"/chr"$SGE_TASK_ID"_minima.pickle
