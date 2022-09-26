## Copyright 2022 University of Washington

## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:

## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#!/bin/bash

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 1-22

## code to run step 5 run it by qsub runStep5.sh
## to allow for superpopulation specific runs we add
## a pop argument

pop=$1

source ../env/bin/activate

python3 ../P03_extract_bpoints.py \
	--name=chr"$SGE_TASK_ID" \
	--subset=fourier_ls  \
	--dataset_path="$pop"/ \
	--input_pickle_fname="$pop"/chr"$SGE_TASK_ID"_minima.pickle \
	> "$pop"/chr"$SGE_TASK_ID".bed
