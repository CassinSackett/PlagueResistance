#following the vignette here: https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
library(vcfR)
setwd("/Users/Loren/pdog_genomics/BTPD_resistant_not/")

# file should contain no missing data and only biallelic SNPs are allowed (MNPs are read as missing)
btpd13vcf <- read.vcfR("btpd_SNP_QC_vcftools_complete.fltr_mean2sd_noMNP13indMAC2.contig.vcf") #bigsnpr package not working for me so I'm going to start with fewer SNPs and run outflank with all
dim(btpd13vcf)
btpd13vcf[1:14,1:14]

bt13geno <- extract.gt(btpd13vcf) # Character matrix containing the genotypes
bt13position <- getPOS(btpd13vcf) # Positions in bp
bt13chromosome <- getCHROM(btpd13vcf) # Chromosome information for each SNP in file (i.e., each scaffold repeated many times)
length(bt13chromosome)

#Paste the SNPs to the scaffolds - this is important if you have SNPs at the same position on different scaffolds
SNPnames <- paste(bt13chromosome,bt13position, sep=":")

btmat13 <- matrix(NA, nrow = nrow(bt13geno), ncol = ncol(bt13geno)) #ncol=14, nrow=7925228
str(btmat13)

newbtmat13 <- t(btmat13) #transpose to get indivs in rows and SNPs in columns for OutFLANK

newbtmat13[bt13geno %in% c("0/0", "0|0")] <- 0
newbtmat13[bt13geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
newbtmat13[bt13geno %in% c("1/1", "1|1")] <- 2

newbtmat13[is.na(newbtmat13)] <- 9
sum(is.na(newbtmat13)) #should be 0.
levels(factor(as.numeric(unlist(newbtmat13)))) #should have 0,1,2,9
str(newbtmat13)
newbtmat13[1:13,1:8]

table(as.vector(newbtmat13)) #MNPs are coded as missing

popmembership <- c("fatality","fatality","survived","survived","fatality","fatality","fatality","survived","survived","survived","fatality","fatality","survived") #with the antibody individual removed


##### calculate FST ######
#example syntax to generate the FST matrix first line below (from https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf)
#bt_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop). This first step takes forever! a few hours with my 7M SNP dataset (15:26 - 16:56 for 6.5M dataset)
library(OutFLANK)
bt_fst <- MakeDiploidFSTMat(newbtmat13, locusNames = SNPnames, popNames = popmembership) #the only way I found to get scaffolds in the results is to add the paste statement above. 

#first check He of data to determine what threshold to set for exclusion
plot(bt_fst$He, bt_fst$FST)

plot(bt_fst$FST, bt_fst$FSTNoCorr)
abline(0,1, col="blue", lwd=2)


#Running outflank. NumSamples is the number of pops. Do this with thinned/ trimmed dataset to estimate the neutral FST distribution
#OutFLANK(FstDataFrame=bt_fst, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.05, NumberOfSamples=3, qthreshold=0.05) #no outliers with these defaults

# Before running the OutFLANK() function to estimate the parameters on the neutral FST distribution, you will want to identify a quasi-independent set of SNPs to calculate FSTbar and df. A common way of obtaining these SNPs is to thin for linkage disequilibrium (SNP thinning).

#vcftools thinning doesn't keep the list subsetted in a way that OutFLANK needs (or I didn't figure out how), so I'll thin here
library(bigstatsr)
library(bigsnpr)

thinmat <- add_code256(big_copy(newbtmat13, type="raw"), code=bigsnpr:::CODE_012)
chromnums <- gsub("contig", "", bt13chromosome) #omit the letters 'contig'
chroms <- as.integer(chromnums) #force contig #s to be integers
bt13thinned <- snp_autoSVD(G=thinmat, infos.chr = chroms, infos.pos = bt13position, roll.size=0, min.mac=3, thr.r2=0.5) #default MAC=10 but with my tiny sample size I need it smaller. size is in kb if infos.pos is provided.  I keep getting the error "Error: Problem with columns [2, 3, 4, 5, 6, 7, 8, 9, 10] of 'U'." but I don't know what U is, and that seems to be all of the columns at least within some range. I found nothing online when googling. I restarted R and this message disappeared but now i get an error "Error in bigutilsr::rollmean(S.col[ind], roll.size) : Parameter 'size' is too large." even when size is 1. It works with roll.size=0 and r2=0.7 but...


bt13thindex <- attr(bt13thinned, which="subset") #indexes of remaining SNPs after pruning
length(bt13thindex) #3945629 with r2=0.7 & mac=2 || 1900509 w/ r2=0.5 & mac=3



#should be with thinned fst 
#outLtest <- OutFLANK(FstDataFrame=bt_fst, LeftTrimFraction=0.01, RightTrimFraction=0.01, Hmin=0.1, NumberOfSamples=2, qthreshold=0.08)
#outLtest$numberHighFstOutliers
#str(outLtest)

out_trimmed <- OutFLANK(FstDataFrame=bt_fst[bt13thindex,], LeftTrimFraction=0.01, RightTrimFraction=0.01, Hmin=0.1, NumberOfSamples=2, qthreshold=0.0775)
out_trimmed$numberHighFstOutliers
str(out_trimmed)

#Plotting neutral distribution of FST:
#We recommend that you plot the distribution of FSTNoCorr for your loci against the inferred distribution. The included function OutFLANKResultsPlotter uses the output of OutFLANK to plot the distribution of FSTNoCorr:

OutFLANKResultsPlotter(out_trimmed, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=FALSE, RightZoomFraction = 0.05, titletext = "Neutral FST distribution") #Zoom=TRUE plots only loci in the highest quantiles of FST, with the proportion determined by RightZoomFraction

hist(out_trimmed$results$pvaluesRightTail)


#now calculate p values for all loci in the dataset 
#Now that we’ve estimated neutral mean FST and df to a quasi-independent set of SNPs, we can go back and calculate P-values and q-values for all the loci in our dataset.
#Note that it is important to run this code with the uncorrected FSTs

pvals <- pOutlierFinderChiSqNoCorr(bt_fst, Fstbar = out_trimmed$FSTNoCorrbar, dfInferred = out_trimmed$dfInferred, qthreshold=0.1, Hmin=0.1)
head(pvals)
pvals[724392,]

quasiouts <- pvals$OutlierFlag==TRUE #for plotting
#quasioutliers <- as.matrix(pvals$OutlierFlag==TRUE)

plot(pvals$He, pvals$FST, pch=19, col=rgb(0,0,0,0.1))
points(pvals$He[quasiouts], pvals$FST[quasiouts], col="seagreen")

#quasioutliers <- pvals$OutlierFlag==TRUE 
#write.table(quasiouts, file="outflank_results_q08.txt", sep="\t") #a list of trues and falses with their index in this file (not that useful)
write.table(pvals, file="outflank_thinned_resultsq08mac3.txt", sep="\t") #gigantic file - 1.3GB for my 6.5M SNPs

#Highlight outliers on Manhattan Plot
#For publication, we want to show the accurate estimate of FST, not the uncorrected estimate. Remember to exclude those low H loci!

plot(pvals$LocusName[pvals$He>0.1], pvals$FST[pvals$He>0.1], xlab="Position", ylab="FST", col=rgb(0,0,0,0.2)) #could do this with SNP name but not scaffold:SNP name. In an amazing twist, they are currently in genomic order
points(pvals$LocusName[quasiouts], pvals$FST[quasiouts], col="seagreen", pch=20)

#Too memory intensive to plot on laptop - should just do this by scaffold
plot(row.names(pvals)[pvals$He>0.1], pvals$FST[pvals$He>0.1], xlab="Position", ylab="FST", col=rgb(0,0,0,0.2)) #could do this with SNP name but not scaffold:SNP name. In an amazing twist, they are currently in genomic order so I can use row names as the x axis
points(pvals$LocusName[quasiouts], pvals$FST[quasiouts], col="seagreen", pch=20)

#Now I'd like to plot by scaffold



