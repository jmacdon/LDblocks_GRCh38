## code to generate  genetic map by interpolating between DeCode positoins

## we use chr22 as a test for this, and for a run of LDetect
## for input we have DeCode data that looks like

## Chr	Begin	End	cMperMb	cM
## chr1	1431813	1515567	0.032449236802513395	0.0027177533791577072
## chr1	1515567	1530002	0.18959655991987387	0.0054545797216010855
## chr1	1530002	1534402	0.004630788906729801	0.005474955192790697
## chr1	1534402	1538787	0.019864562217819185	0.005562061298115834
## chr1	1538787	1541864	0.0003320830655916105	0.00556308311770866

## and the goal is to interpolate the cM for the SNPs in the 1000G VCF files, which in the
## end should look like

## chr22	15287922	1.457757
## chr22	16370978	1.585395
## chr22	16372046	1.585931
## chr22	16373044	1.586426
## chr22	16373740	1.586781
## chr22	16374612	1.587307

## where the second column is  physical position and the
## third column is the interpolated cM from DeCode

## basically we have the beginning and end of each interval from DeCODE, and
## also the cM/Mb, so assuming linearity, that can be converted to cM/Bp.
## We also have the distance from the start of the interval and the SNP position,
## so the cM for a given position is just starting cM of interval + (distance from start)*cM/Bp

## Any position < first DeCODE position gets a 0, and any position > last DeCODE position gets the max cM value


library(GenomicRanges)
decode <- read.table("aau1043_datas3.gz", header = TRUE)

getPos <- function(chr){
    file <- paste0("vcfsubsets/", chr, ".vcf.gz")
    pos <- system(paste("zcat", file, "| awk '$1 !~ /#/ {print $2}'"), intern = TRUE)
    gr <- GRanges(rep(chr, length(pos)), IRanges(as.numeric(pos), width = 1))
    gr
}


decodeGr <- function(decodedat, chr, chrgr, try.a.sub = FALSE){
    tmp <- decodedat[decodedat$Chr == chr,]
    if(try.a.sub) tmp <- tmp[1:10,]
    ## adjust to have decode data start at zero and end at the max extent of the data
    tmp[1,2] <- 0
    if(start(chrgr)[length(chrgr)] > tmp[nrow(tmp),3])
       tmp[nrow(tmp),3] <- start(chrgr)[length(chrgr)] + 1
    gr <- GRanges(tmp[,1], IRanges(tmp[,2], width = 1), cMperMb = tmp[,4], cM = tmp[,5])
    gr
}



interpThat <- function(inds, snpgr, decodeblocks) {
    blocklst <- split(start(snpgr), inds)
    blocks <- unique(inds)
    decodeblocks <- decodeblocks[blocks]
    cms <- do.call(c, lapply(seq(along = blocklst), function(x) {
                          cMperBp <- decodeblocks$cMperMb[x]/1e6
                          if(x == 1L)
                              startcM <- 0
                          else
                              startcM <- decodeblocks$cM[x-1]
                          out <- startcM + (blocklst[[x]] - start(decodeblocks[x])) * cMperBp
                      }))
    cms
}

writeOut <- function(dat, fname){
    f <-  gzfile(fname, "w")
    write.table(dat, f, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    close(f)
}


## master function

runIt <- function(chr){
    fname <- paste0(chr, ".tab.gz")
    gr <- getPos(chr)
    decgr <- decodeGr(decode, chr, gr)
    inds <- follow(gr, decgr)
    cms <- interpThat(inds, gr, decgr)
    gr$cM <- cms
    out <- as(gr, "data.frame")[,c(1,2,6)]
    writeOut(out, fname)
}


## run them all
for(i in paste0("chr", 1:22)) runIt(i)
