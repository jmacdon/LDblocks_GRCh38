#!/bin/sh

#$ -q lindstroem.q
#$ -cwd
#$ -l h_vmem=20G
#$ -t 101-105


source ../venv/bin/activate
export PYTHONPATH=/projects/lindstroem/.py-site-packages/

fn=$1
chr=${fn/_partitions/}
chr=${chr/scripts\//}
vcf=$2
map=$3
indfile=$4
pop=$5
popsize=$6

## get the row of the partition file and read into an array
inline=(`awk "NR==$SGE_TASK_ID" $fn`)
start=${inline[0]}
stop=${inline[1]}


tabix -h "$vcf" "$chr":"$start"-"$stop" | python3 ../P00_01_calc_covariance.py "$map" "$indfile" "$popsize" 1e-7 "$pop"/"$chr"/"$chr"."$start"."$stop".gz 
