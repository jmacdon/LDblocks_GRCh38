#!/bin/bash

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 1-22

## code to run step 3 run it by qsub runStep3.sh <population dir>
## e.g., qsub ./runStep3.sh EUR
pop=$1

source env/bin/activate

python3 P01_matrix_to_vector_pipeline.py --dataset_path="$pop"/ --name=chr"$SGE_TASK_ID" --out_fname="$pop"/chr"$SGE_TASK_ID".vector.txt.gz
