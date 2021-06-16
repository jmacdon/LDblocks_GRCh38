#!/bin/bash

## subset all the VCFs and index#$ -cwd


#$ -l h_vmem=20G
#$ -q lindstroem.q
#$ -t 1-22
#$ -cwd



zcat ../chr"$SGE_TASK_ID".tab.gz | cut -f 1-2 > chr"$SGE_TASK_ID"

bcftools view -R chr"$SGE_TASK_ID" \
	 /projects/lindstroem/1000Genomes/hg38/ALL.chr"$SGE_TASK_ID".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz \
	 -O z -o chr"$SGE_TASK_ID".vcf.gz
bcftools index chr"$SGE_TASK_ID".vcf.gz 

