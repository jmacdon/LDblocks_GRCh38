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

## run runCov on all chromosomes
## Modified to use different populations so we need
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
    ## most chrs are over 100 partitions; only run <= 100 at a time
    nlines=`wc -l $f | awk '{print $1}'`
    nblocks=$(($nlines/100))
    nblocks=${nblocks/.*}
    if [ $i -eq 0 ]
    then
	if [ $nblocks -eq 0 ]
	   ## there are fewer than 100 blocks so just run it
	then
	    ## Adjust runCov.sh to reflect the number of blocks to run
	    sed "s/START/1/" runCov.blank > runCov.sh
	    sed -i "s/END/$nlines/" runCov.sh
	    qsub -N covrun"$((i+1))" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
	    ((i++))
	else
	    ## there are more than 100 blocks, so run 100 at a time and use hold_jid to make remainder wait
	    for (( j=0; j<=$nblocks; j++ ))
	    do
		if [[ $j -lt $nblocks ]]
		then
		   sed "s/START/$((j*100+1))/" runCov.blank > runCov.sh
		   sed -i "s/END/$(((j+1)*100))/" runCov.sh
		   qsub -N covrun"$((i+1))" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
		   ((i++))
		else
		    nextrun=$((i+1))
		    sed "s/START/$((j*100+1))/" runCov.blank > runCov.sh
		    sed -i "s/END/$nlines/" runCov.sh
		    qsub -hold_jid covrun$i -N covrun"$nextrun" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
		    ((i++))
		fi
	    done
	fi
    else
	## After the first set, we just put in the queue and use hold_jid to make the process wait until the previous one finishes
	if [ $nblocks -eq 0 ]
	then
	    nextrun=$((i+1))
	    sed "s/START/1/" runCov.blank > runCov.sh
	    sed -i "s/END/$nlines/" runCov.sh
	    qsub -hold_jid covrun$i -N covrun"$nextrun" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
	    ((i++))
	else
	    for (( j=0; j<=$nblocks; j++ ))
	    do
		if [[ $j -lt $nblocks ]]
		then
		    nextrun=$((i+1))
		    sed "s/START/$((j*100+1))/" runCov.blank > runCov.sh
		    sed -i "s/END/$(((j+1)*100))/" runCov.sh
		     qsub -hold_jid covrun$i -N covrun"$nextrun" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
		    ((i++))
		else
		    nextrun=$((i+1))
		    sed "s/START/$((j*100+1))/" runCov.blank > runCov.sh
		    sed -i "s/END/$nlines/" runCov.sh
		    qsub -hold_jid covrun$i -N covrun"$nextrun" ./runCov.sh $f vcfsubsets/"$chrsub".vcf.gz "$chrsub".tab.gz "$indfile" "$pop" "$popsize"
		    ((i++))
		fi
	    done
	fi
    fi
done

       
		       
	


    
