#!/bin/bash

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 1-22

bcftools view --min-af 0.01 --types snps -o vcfsubsets/chr"$SGE_TASK_ID".vcf.gz --output-type z \
	 /projects/lindstroem/1000Genomes/hg38/ALL.chr"$SGE_TASK_ID".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools index vcfsubsets/chr"$SGE_TASK_ID".vcf.gz
	 
