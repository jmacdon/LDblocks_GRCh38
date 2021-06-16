#!/bin/bash

#$ -cwd
#$ -q lindstroem.q
#$ -l h_vmem=20G

source ./venv/bin/activate
export PYTHONPATH=/projects/lindstroem/.py-site-packages/

fn=$1
map=$2
chr=$3
vcf=$4
linenum=$SGE_TASK_ID

echo "$fn" "$map" "$chr" "$vcf" "$linenum"


# i=1
# cat $fn | while read start stop || [[ -n $line ]];
# do
#     if [[ $i -eq $linenum ]]
#        tabix -h "$vcf" "$chr":"$start"-"$stop" | python3 P00_01_calc_covariance.py "$map" eurinds.txt 11418 1e-7 "$chr"/"$chr"."$start"."$stop".gz
#     fi
# done
