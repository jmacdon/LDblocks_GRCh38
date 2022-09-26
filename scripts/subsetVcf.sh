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

bcftools view --min-af 0.01 --types snps -o vcfsubsets/chr"$SGE_TASK_ID".vcf.gz --output-type z \
	 /projects/lindstroem/1000Genomes/hg38/ALL.chr"$SGE_TASK_ID".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools index vcfsubsets/chr"$SGE_TASK_ID".vcf.gz
	 
