#following the vignette here: https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
library(vcfR)
setwd("/Users/Loren/pdog_genomics/BTPD_resistant_not/")

btpd13vcf <- read.vcfR("btpd_SNP_QC_vcftools_complete.fltr_mean2sd_6surv7fatal.recode.vcf.gz")
dim(btpd13vcf)
btpd13vcf[1:14,1:14]

bt13geno <- extract.gt(btpd13vcf) # Character matrix containing the genotypes
bt13position <- getPOS(btpd13vcf) # Positions in bp
bt13chromosome <- getCHROM(btpd13vcf) # Chromosome information for each SNP in file (i.e., each scaffold repeated many times)

btmat13 <- matrix(NA, nrow = nrow(bt13geno), ncol = ncol(bt13geno)) #ncol=14, nrow=7925228

newbtmat13 <- t(btmat13) #transpose to get indivs in rows and SNPs in columns for OutFLANK

newbtmat13[bt13geno %in% c("0/0", "0|0")] <- 0
newbtmat13[bt13geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
newbtmat13[bt13geno %in% c("1/1", "1|1")] <- 2

newbtmat13[is.na(newbtmat13)] <- 9
sum(is.na(newbtmat13)) #should be 0.
levels(factor(as.numeric(unlist(newbtmat13)))) #should have 0,1,2,9


table(as.vector(newbtmat13)) #I wonder if the 9s are MNPs

#popmembership <- c("fatality","fatality","survived","survived","fatality","antibodies","fatality","fatality","survived","survived","survived","fatality","fatality","survived")
popmembership <- c("fatality","fatality","survived","survived","fatality","fatality","fatality","survived","survived","survived","fatality","fatality","survived") #with the antibody individual removed


######### generate a file of pruned SNPs ############
#btpd13thinvcf <- read.vcfR("btpd_SNP_QC_vcftools_complete.fltr_mean2sd_6surv7fatal.thin5k.vcf") #vcftools thinning doesn't keep the list it subsetted in a way that OutFLANK needs :(
library(bigstatsr)
library(bigsnpr)

thinmat <- add_code256(big_copy(t(btmat13$G), type="raw"), code=bigsnpr:::CODE_012)
bt13thinned <- snp_autoSVD(G=thinmat, infos.chr = btmat13$chromosome, infos.pos = btmat13$position)
bt13thindex <- attr(bt13thinned, which="subset") #indexes of remaining SNPs after pruning
length(bt13thindex)


#example syntax to generate the FST matrix first line below (from https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf)
#bt_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop). Takes forever!
library(OutFLANK)
bt_fst <- MakeDiploidFSTMat(newbtmat13, locusNames = bt13position, popNames = popmembership) 

thinned_fst <- MakeDiploidFSTMat(btmat13thin, locusNames = bt13positionthin, popNames = popmembership) 


#Running outflank. NumSamples is the number of pops. I'd like to just get rid of the antibodies indiv but it's in the vcf file. First run below is with default settings
#OutFLANK(FstDataFrame=bt_fst, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.05, NumberOfSamples=3, qthreshold=0.05) #no outliers with these defaults

#first check He of data to determine what threshold to set for exclusion
plot(thinned_fst$He, thinned_fst$FST)

plot(bt_fst$FST, bt_fst$FSTNoCorr)
abline(0,1, col="blue", lwd=2)

## ********************************************************* ##
# Before running the OutFLANK() function to estimate the parameters on the neutral FST distribution, you will want to identify a quasi-independent set of SNPs to calculate FSTbar and df. A common way of obtaining these SNPs is to thin for linkage disequilibrium (SNP thinning).

outLtest <- OutFLANK(FstDataFrame=thinned_fst, LeftTrimFraction=0.01, RightTrimFraction=0.1, Hmin=0.1, NumberOfSamples=2, qthreshold=0.2)
str(outLtest)
outLtest$numberHighFstOutliers

#Plotting neutral distribution of FST:
#We recommend that you plot the distribution of FSTNoCorr for your loci against the inferred distribution. The included function OutFLANKResultsPlotter uses the output of OutFLANK to plot the distribution of FSTNoCorr:

OutFLANKResultsPlotter(outLtest, withOutliers=TRUE, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=FALSE, RightZoomFraction = 0.05, titletext = "Neutral FST distribution") #Zoom=TRUE plots only loci in the highest quantiles of FST, with the proportion determined by RightZoomFraction

hist(outLtest$results$pvaluesRightTail)


#now calculate p values for all loci in the dataset 
#Now that we’ve estimated neutral mean FST and df to a quasi-independent set of SNPs, we can go back and calculate P-values and q-values for all the loci in our dataset.
#Note that it is important to run this code with the uncorrected FSTs

pvals <- pOutlierFinderChiSqNoCorr(bt_fst, Fstbar = outLtest$FSTNoCorrbar, dfInferred = outLtest$dfInferred, qthreshold=0.14, Hmin=0.15)

plot(pvals$He, pvals$FST, pch=19, col=rgb(0,0,0,0.1))
points(pvals$He[quasioutliers], pvals$FST[quasioutliers], col="seagreen")

#quasioutliers <- pvals$OutlierFlag==TRUE 
#write.table(quasioutliers, file="outflank_top2127.txt", row.names=FALSE, quote=FALSE, sep="\t")

#Highlight outliers on Manhattan Plot
#For publication, we want to show the accurate estimate of FST, not the uncorrected estimate. Remember to exclude those low H loci!

# need to change alllloci to pvals
plot(allloci$LocusName[allloci$He>0.1], allloci$FST[allloci$He>0.1], xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(allloci$LocusName[quasioutliers], allloci$FST[quasioutliers], col="seagreen", pch=20)






