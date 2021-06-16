#!/bin/bash

## run runCov on all chromosomes
## Modifiec to use different populations so we need
## to pass in the file containing the individuals and the
## output population dir
i=0
indfile=$1
pop=$2
popsize=$3

## to do this by superpopulation we need to move the partitions..

mkdir -p "$pop"
cp -r scripts/ "$pop"/

for f in scripts/*partitions
do
    chrsub=${f/_partitions/}
    chrsub=${chrsub/scripts\//}
    mkdir -p "$pop"/"$chrsub"
    if [[ $i -eq 0 ]]
    then
	nlines=`wc -l $f | awk '{print $1}'`
	sed "s/CHANGEME/$nlines/" runCov.blank > runCov.sh
	#echo "should run qsub  -N covrun${i+1} ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz"
	qsub -N covrun"${i+1}" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
	((i++))
    else
	nextrun=$((i+1))
	nlines=`wc -l $f | awk '{print $1}'`
	sed "s/CHANGEME/$nlines/" runCov.blank > runCov.sh
	#echo "should run qsub -hold_jid covrun$i -N covrun$nextrun ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz"
	qsub -hold_jid covrun$i -N covrun$nextrun ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
	((i++))
    fi
done

    
