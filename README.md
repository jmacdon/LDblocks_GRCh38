# LDblocks_GRCh38 <a name="top"/>

The goal of the code in this repository is to generate approximately
independent LD blocks based on the GRCh38 genome. While there are
existing LD blocks (for example
[here](https://github.com/bogdanlab/RHOGE)), they are based on
GRCh37. Methods to convert genetic loci from one genome build (notably
the UCSC liftOver tool) do not work well for genetic blocks, tending
to fragment the block and often spreading portions across different
chromosomes.

Instead, we use
[LDetect](https://bitbucket.org/nygcresearch/ldetect/src/master/) to
generate new LD blocks, using GRCh38-based [1000Genomes
data](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/)
and new sex-averaged dense recombination rate data from either [deCode
Genetics](https://www.science.org/doi/suppl/10.1126/science.aau1043/suppl_file/aau1043_datas3.gz)
or [pyrho](https://github.com/popgenmethods/pyhro) (Spence and Song [2019](https://www.science.org/doi/10.1126/sciadv.aaw9206)).


### Table of contents
+ [Overview](#justreaditdude)
+ [Pipeline prep](#pipelineprep)
+ [Parse data](#parseit)
+ [Process data](#processit)
+ [Technical addendum](#ughbro)

## Overview <a name=justreaditdude/>

A high level summary of the process that we used is as follows:

+ Install LDetect and bcftools/tabix
+ Download and parse the 1000Genomes VCF files
  + Bash code below
  + subsetVcf.sh script (in scripts dir of this repository)
+ Interpolate recombination rates onto 1000G variants
+ Generate LD blocks using a set of bash scripts
  + Partition chromosomes (bash code below)
  + Calculate covariance matrices (runAllCov.sh in scripts dir)
  + Convert covariance matrices to vector (runStep3.sh in scripts dir)
  + Calculate minima (runStep4.sh in scripts dir)
  + Output bed files (runStep5.sh in scripts dir)
  

## Pipeline prep <a name="pipelineprep"/>

These data were generated on a Linux cluster, as it is beneficial to
parallelize many of the steps. We installed LDetect in a virtualenv
called 'env' following the recommendation to use `pip`

```sh
virtualenv env
source env/bin/activate
pip install ldetect
```

This installs the LDetect software, as well as downloading a set of
example scripts that can be used to perform many of the required
steps. We also used [bcftools](http://www.htslib.org/download/) and
[tabix](http://www.htslib.org/download/) which is part of htslib, to
process the VCF files.


## Parse data <a name="parseit"/>

We downloaded the December 2018 biallelic SNV 1000Genomes GRCh38 VCF
files from the link noted above, as well as the GRCh38 recombination
maps from the Halldorsson Science paper (deCode link above). Please
note that this file has a gz extension, but is not compressed. The
recombination data file contains five columns; the chromosome, the
start and end of the interval (in bases), the recombination rate
within the interval (in cM/Mb) and the recombination rate at the end
of the interval (cM). LDetect expects a three-column file containing
the chromosome, variant position, and recombination rate. To generate
that file we first subsetted the VCF files for each chromosome to only
include SNPs with a MAF>0.01, using the subsetVcf.sh bash
script. *Please note that most of the bash scripts include paths
unique to our compute server. We leave the paths as is, with the
expectation that anybody wanting to replicate will modify to suit.*

The subsetVcf.sh script is run under `qsub` and uses `bcftools` to
subset the VCF and index, parallelized over all chromosomes. 

```sh

qsub ./subsetVcf.sh


```


After generating subsetted VCF files, we interpolated data from deCode
and pyhro genetic map files using an R script (interpolate.R or
interpolate_pyrho.R), which generates the expected format for
LDetect. This file is self-contained and can be called at the command
line, using `R --vanilla < interpolate.R` (or interpolate_pyhro.R,
depending on the source data) to output gzipped genetic map files for
each chromosome. If the Bioconductor `GenomicFeatures` package is not
installed, it will be automatically installed. If the
recombination map data are not in the working directory, they will be
downloaded automatically in the case of deCODE, or an error will be
issued with a link to manually download and uncompress the pyhro data.



## Process data <a name="processit"/>

To generate the LD blocks we followed the general instructions
provided at the [LDetect bitbucket
repository](https://bitbucket.org/nygcresearch/ldetect/src/master/)
with some small changes required by the particulars of our data. The
first step is to generate partitions of the genome that can then be
processed in parallel.

```sh
python3 P00_00_partition_chromosome.py <genetic_map> <n_individuals> <output_file>
```

Because each step is fast, we did it sequentially. The number of
individuals is the number of individuals in the reference panel, which
in our case is 417.

```sh
for f in chr{1..22}.tab.gz; do python3 P00_00_partition_chromosome.py $f 417 scripts/${f/.tab.gz/}_partitions; done
```

Note that the LDetect Python scripts expect things to be in a `scripts/`
directory, so we followed that convention. The next step is to use
these partitions to compute the covariance matrix. In other words, we
want the covariance of the variants within a chromosomal region, for
those individuals of a given ancestry. The VCF files we have include
all of the 1000Genomes subjects, and we only want subjects with
a given genetic ancestry. For the deCODE data we selected individuals from the
following sub-populations (TSI, IBS, CEU, GBR) using a bash script:

```sh
for dir in TSI IBS CEU GBR 
 do
     curl -s -L ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/"$dir"/ \
 | awk '{print $NF}' >> eurinds.txt
 done
```
Which collects all the 1000Genomes subject IDs from those
sub-populations in a text file. There are more subjects in
this text file than we have in the VCF files, so we used a simple
script to get the intersection:

```sh

zcat vcfsubsets/chr9.vcf.gz | head -n 30 | awk '$1 ~/#CHR/ {print $0}' \
| cut -f 10- > subjects.txt

cat subjects.txt eurinds.txt | sort | uniq -d > tmp.txt; mv tmp.txt eurinds.txt

```

The first line simply captures the subject IDs from one of the VCF
files (they all have the same IDs) and the second returns all the
European subject IDs that are found in both the VCFs and the IDs we
got from 1000Genomes.

We used similar scripts for the pyhro recombination maps, depending on
the goal. For sub-population specific LD blocks, we only selected
those subjects from a given sub-population. For super-population maps,
we similarly chose all sub-populations within a super-population to
compute correlations. But do note that we used a single sub-population
for the genetic map (e.g., for the AFR super-population, we used the
GWD sub-population genetic map).

For EUR, we omitted FIN, due to large differences between FIN and
other European ancestries. For AFR we omitted both ASW and ACB due to
admixture. Block statistics for the four super-populations using pyrho
recombination maps can be found [here](data/popTables.html).

We then computed all the correlation values in parallel using
`runAllCov.sh` which can be found in the scripts directory.

```sh
runAllCov.sh eurinds.txt EUR 11418
```

The final argument for that script (11418) is the effective population
size for Europeans. This script is meant to limit the number of
processes sent to the compute nodes to 100 or less. The remaining
processes are held in the queue and sent for processing only after the
preceding run has finished. 

This step is the most computationally expensive, taking several days
to finish.

There are three more steps; convert the covariance matrices into
vectors, calculate the minima across the covariance matrices, and then
extract the minima (which are output as python .pickle files) into .bed
files. We used three scripts to parallelize these steps (which are
steps 3-5 in the LDetect example).

```sh

qsub ./runStep3.sh EUR
qsub ./runStep4.sh EUR
qsub ./runStep5.sh EUR

``` 

## Technical addendum <a name=ughbro/>

This repository is mainly meant as a way for people to get the LD
blocks for their own use, as well as to document exactly how the LD
blocks were generated. If you just want the LD blocks, see the data
directory of this repository. However, it is not inconceivable that
someone might want to use the scripts to generate their own LD
blocks. To that end, we provide some technical pointers.

First, please note that the bash scripts provided are intended to be
used on a cluster environment, which usually means that `qsub` is
available. They will not work 'out of the box', as we have hard-coded
the queue that we used, as well as some of the paths to our
data. These will have to be adjusted to correspond to your own cluster
in order for them to work correctly. For example, here is
`runStep3.sh`:

```sh
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

```
We used our internal queue (the `#$ -q lindstroem.q` directive). This
will have to be modified to use your own queue.

Second, we installed LDetect in a `virtualenv` called 'env', and the
bash scripts (notably runStep3.sh, runStep4.sh and runStep5.sh) have a
line that activates that `virtualenv` prior to running the python code
(the line `source env/bin/activate`). The easiest thing to do would be
to copy that paradigm.

Third, there appears to be a bug in the LDetect codebase that causes
a problem in step 4, where the minima are computed. This only affects
the uniform local minima, which might not matter (we used the
fourier-ls breakpoints which correspond to the low-pass filter with
local search). There is a cryptic line in the README for LDetect that
says

> This file (P02_minima_pipeline.py) can be tweaked to remove all but
> the low-pass filter with local search algorithm in order to reduce
> total runtime.

What we had to do to 'fix' the problem is complicated and boring and
not worth going into. Suffice it to say that lines 103-105 in
`P02_minima_pipeline.py` look like this:

```python

metric_out_uniform_local_search = apply_metric(chr_name, begin, end, config, breakpoint_loci_uniform_local_search['loci'])
flat.print_log_msg('Global metric:')
print_metric(metric_out_uniform_local_search)

```
And commenting out those lines will bypass the bug. To follow the
recommendation from the LDetect authors ('remove all but the low-pass
filter...'), one would also comment out lines 67-78:

```python

# METRIC FOR UNIFORM BREAKPOINTS
flat.print_log_msg('* Calculating metric for uniform breakpoints...')
# step = int((end-begin)/(len(breakpoint_loci)+1))
# breakpoint_loci_uniform = [l for l in range(begin+step, end-step+1, step)] 
step = int(len(init_array_x)/(len(breakpoint_loci)+1))
breakpoint_loci_uniform = [init_array_x[i] for i in range(step, len(init_array_x)-step+1, step)]

# metric_out_uniform = apply_metric(chr_name, begin, end, cnst.const[dataset], breakpoint_loci_uniform)
metric_out_uniform = apply_metric(chr_name, begin, end, config, breakpoint_loci_uniform)
flat.print_log_msg('Global metric:')
print_metric(metric_out_uniform)

```


