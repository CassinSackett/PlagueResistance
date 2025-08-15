[![hackmd-github-sync-badge](https://hackmd.io/g2NaGk4fQx2kg8yPS1Icbw/badge)](https://hackmd.io/g2NaGk4fQx2kg8yPS1Icbw)
###### tags: `selection`, `GWAS`, `outlier`
:ear_of_rice: :dog2: 

# Tests for selection on plague resistant BTPD :chipmunk: 

>The goal of this pipeline is to detect evolutionary signatures of selection that will enable the determination of the genetic basis of a phenotype (e.g., plague resistance in prairie dogs). The pipeline starts with a [filtered vcf file, which we generated previously](https://hackmd.io/@sacketlc/HkEAgfbmK), and conducts tests for selection within genomes in a small sample.

- Table of Contents
[ToC]


## I. GWAS in plink
:::warning
This code is for plink v1.9. I will update the pipeline in the future for plink2.
:::


### 1. Convert vcf to plink format
First, we need to convert our .vcf file to plink2 format. This is surprisingly fast :woman-running: 

```
./plink2 --vcf ../btpd_SNP_QC_vcftools_complete.fltr_mean2sd_MAC2.recode.vcf --out btpd_SNPqc_complete_fltr_mean2sd_MAC2 --allow-extra-chr
``` 

This will generate several output files with the same prefix; you will supply plink2 with the prefix and it will use all of these files.


### 2. Generate phenotype file for plink
You need to tell plink2 which individuals to group in each phenotype group. You will generate a tab-delimited file that looks like this:

``` 
S_19A26070701_DSW64907-V	S_19A26070701_DSW64907-V	1
S_47A13080701_DSW64908-V	S_47A13080701_DSW64908-V	1
S_47A15080703_DSW64910-V	S_47A15080703_DSW64910-V	1
D_17A25060701_DSW64899-V	D_17A25060701_DSW64899-V	2
D_17A14070804_DSW64904-V	D_17A14070804_DSW64904-V	2
D_17a29060701_DSW64902-V	D_17a29060701_DSW64902-V	2
``` 
You can also have a header, which can be useful if you have multiple phenotypes, and you just have to include the ```--pheno-name``` flag in your plink command. The fields in my file are FAMILY_ID, WITHIN_FAMILY_ID, and PHENOTYPE.

In my case, you can see that each individual is written twice, as its name is read in the vcf file, and then its phenotype follows in the third column. I did this because I was not interested in among-family variation, but this could be useful to include in stratified populations. In this example, the individuals whose names begin with 'S' survived :smiley: (I assigned phenotype 1) and those with names beginning with 'D' died :skull_and_crossbones: (I assigned phenotype 2). In plink's syntax, 1 is a control and 2 is a case.

### 3. Edit your .map file to get correct scaffold names
This is an optional step, but it will make downstream analyses easier. Since I do not have a chromosome-level assembly, my .map file looked like this:

``` 
0	scaffold1:11948	0	11948
0	scaffold1:12214	0	12214
0	scaffold1:13653	0	13653
0	scaffold1:14248	0	14248
``` 

where in the first column, chromosome number was replaced by the conversion with a zero because plink's conversion did not recognize the chromosome names. To fix this, I did the following:

``` 
cut -f2- original.map > cut.map
``` 
followed by


``` 
cat cut.map | awk -F '\t' '{split($1, a, ":"); print a[1] "\t" a[2] "\t" $2 "\t" $3}' > fixed.map
``` 
This will use awk to split column 1 on the ```:``` and then print both halves of that column as new separate columns, followed by all the other columns that were in the file.

The new file looks like this:
``` 
scaffold1	11948	0	11948
scaffold1	12214	0	12214
scaffold1	13653	0	13653
scaffold1	14248	0	14248
``` 
where the columns are chrom, identifier, distance in cM (0 as default), position.


### 4. Decide which test to use in plink

:eyes: the syntax in plink2 has changed and I need to update this pipeline to reflect  that. :eyes:

#### plink v1.9:
You need to decide which test to perform (e.g., allele only or genotype). Be aware that, for instance --assoc and --model will produce different results. You also need to decide which method to use to statistically correct for multiple tests. I tend to use Genomic Control for additive models. I use ```--assoc``` , along with their 95% (or your chosen level) confidence intervals ```--ci 0.95``` and a flag to calculate adjusted p values ```--assoc``` like this:

``` 
/Users/Loren/plink_mac_20210606/plink --file btpd_SNP_QC_vcftools_complete.fltr_mean2sd --allow-extra-chr --allow-no-sex --assoc --adjust --pheno ../btpd_SNP_QC_vcftools_complete.fltr_mean2sd_phenos.txt --ci 0.95 --out btpd_SNP_QC_vcftools_complete.fltr_mean2sd_assocCI95
``` 

You may need to remove the ```--adjust``` to get the full model output with the following fields:
```
 CHR    SNP    BP    A1    F_A    F_U    A2     CHISQ    P    OR    SE    L95    U95 

```
Note: A1 is the name of the minor allele *in the sample* (i.e., the specific allele could change with different population membership or adding/removing individuals) and A2 is the name of the major allele in the sample. F_A and F_U are the frequency of the minor allele in cases (affected individuals) and controls (unaffected individuals), respectively. OR is the estimated odds ratio for A1 (i.e., assumes A2 is the  reference allele).

It is also important to look at the odds ratios, as well as loci that have certain genotypes strongly associated with the phenotype of interest but that cannot have OR calculated because the allele frequency of the minor allele is 0 in the unaffected / control group. (This will be less likely in larger samples.)

#### plink2:
 :construction: THIS IS STILL IN DEVELOPMENT :construction: 
 
First, convert plink1 files (.ped + .map) to plink2 format (.pgen + .psam + .pvar) (this is fast)
```
/path/to/plink2 --pedmap /path/to/prefix_of_ped_map-files --allow-extra-chr --out output_prefix
```

:one:  The first step is to calculate allele frequencies in each population so that we can account for them in the GWAS (I did not do this for the natural epizootic because I had only 4 samples per population).
```
/path/to/plink2 --pfile /path/to/prefix/file --allow-extra-chr --pheno /path/to/phenotype/file.txt --freq --out /path/to/prefix/file

```
This will create a table of allele frequencies output as a .afreq file.

:two:  Next, we run a PCA to account for population structure when running the GWAS 

```
/path/to/plink2 --pfile /path/to/prefix/file --allow-extra-chr --pheno /path/to/phenotype-file.txt --read-freq /path/to/file.acount --pca --out output
```

:three:  Now we are ready to run the association test! We need to include our top PCs in the association test as covariates using ```--covar``` and specify the number of PCs with ```--parameters``` to account for population structure. It may be a good idea to run a PCA in a program that produces a scree plot. We will also include our allele frequency file.

```
/path/to/plink2 --pfile /path/to/file-prefix --allow-extra-chr --pheno /path/to/phenotype-file.txt --read-freq /path/to/file.afreq --glm --covar /path/to/file.eigenvec  --parameters 1-4 --out /path/to/prefix/file
```

If you are interested in the genotypic model 
``` 
plink2 --pfile mydata \
       --glm genotypic \
       --parameters 1-2

``` 
should induce
Additive effect ('ADD') (parameter 1) and
Dominance deviation ('DOMDEV') (parameter 2) (or do I need to specify --glm genotypic interaction for this to work?).

You could replace "genotypic" with dominant, recessive, hetonly and interaction.


### 5. Manage the plink output
I like to edit the original output from plink to get rid of the 0s it assigns my chromosomes, remove the SE column, etc.

First, change to tab-delimited using translate (tr) and make just one tab between columns (s for squeeze all consecutive delimiters into one):

```
tr -s ' ' '\t' < plink.out > my.out
```


I now have output files for the sites with the highest OR, the lowest OR, the most significant associations, and the highest chi-square. I want to combine those into a single file and look for overlap. 

First, I tag each file by adding a column at the end for the specific analysis:

```
sed 's/$/ \tORmas40/' btpd_ORmas40.txt > btpd_ORmas40.annotated.txt 
```
and then I can concatenate the files together with ``` cat btpd_ORmas40.annotated.txt btpd_pless0.001.annotated.txt btpd_x2mas18.annotated.txt > btpd_allcriteria_topcandidates.txt ```

### 6. Generate random GWAS data to provide null expectations

I am sure this could be done much more elegantly in [SLiM](https://messerlab.org/slim/). However, my approach was to randomly generate case/ control groups 1000 times by randomizing the phenotypes among individuals, and then running a GWAS on those groups, and pulling out the significant loci. This way, I can get an idea of how many false positives (significant loci in random groupings) to expect, and can use these null-model candidate loci downstream.

*Edit: on a Mac OSX, the 'sort -R' did not work* 

For plink2:
First, generate random groupings of individuals:
```
for i in {1..1000}; do sort -R ../btpd_SNP_QC_vcftools_complete.fltr_mean2sd_bin.fam > btpd_SNP_QC_vcftools_complete.fltr_mean2sd_binv$i.fam; head -n 7 btpd_SNP_QC_vcftools_complete.fltr_mean2sd_binv$i.fam > pop1v$i.txt; tail -n 7 btpd_SNP_QC_vcftools_complete.fltr_mean2sd_binv$i.fam > pop2v$i.txt; done
```
then you can do this:
```
plink1.9.exe --bfile file --keep set1.fam --make-bed --out subset1

plink1.9.exe --bfile file --keep set2.fam --make-bed --out subset2
```

and continue with analyses.

For plink1.9, I want to randomize the phenotype in column 3 of my .pheno file, which looks like this
R_47A13080701_DSW64908-V        R_47A13080701_DSW64908-V        1
R_47A15080703_DSW64910-V        R_47A15080703_DSW64910-V        1
R_17A25060701_DSW64899-V        R_17A25060701_DSW64899-V        2
R_17A14070804_DSW64904-V        R_17A14070804_DSW64904-V        2
...






## II. Outlier FST
I am more confident that my candidates are real when using several different programs to detect outliers. Two methods I like based on ease and widespread use are 1) site-by-site and windowed FST in vcftools and 2) BayPass.

### 1. Calculate FST at each site in  [VCFtools](https://vcftools.github.io/man_latest.html) 
to calculate Fst at each site so we can determine which loci fall outside the top 1% (or appropriate cutoff) of genome-wide FST values. VCFtools uses the Weir & Cockerham (1984) estimator. First, we create two popfiles, one with each population. These popfiles are simply text files with the list of individuals, and I use the ```--missing-indv``` flag in VCFtools to get the sample names as they are understood by the software (e.g., in this case, I omitted the .sorted_markdup.g.vcf from each individual). 

```
/sackettl/software/vcftools/bin/vcftools --vcf btpd_SNPQC_complete_DPmean2sd.vcf --missing-indv
``` 

VCFtools uses Weir & Cockerham's estimate, controlling for sample size, so many datasets have negative estimates. I change these to zero immediately (but you could do so later in R) with a python script '[convert_negativeFSTs_to_0.py](https://github.com/CassinSackett/SNP_capture/blob/master/convert_negative-fst_to_zero.py)'. 

Download the python script and edit the output filename (line 3), input file name (line 41), and the column number that you want to adjust (lines 19, 23, 29). Python is zero-indexed, so if you want to change column 3, there should be a [2] in between the brackets. Finally, you may need to change the minimum line length (line 15) if you have fewer columns (len refers to the number of columns).

Next, I sort from highest to lowest FST with ```sort -n -k3,3 surv-fat.weir.fst > surv-fat_sorted.fst``` and scroll through to look for obvious breaks in the FST values, and I plot the distribution of FST to confirm this.


### 2. Calculate FST in windows using  [VCFtools](https://vcftools.github.io/man_latest.html) 

First, decide on the window size you wish to use for calculated FST. This may depend on how your data will look when plotting (e.g., plotting FST across a large number of scaffolds may be easier to visualize if FST windows are larger) or on some characteristic of your taxon (e.g., size of LD blocks).



### 3. Test for outliers in BayeScan
:warning: My BayeScan pilot runs kept getting stuck at ~14% no matter how much memory I allocated. I am not sure if this was a specific problem with our cluster, but I ultimately decided not to use this analysis, so this code is not final.  :warning: 

First, we need to convert our vcf file to a BayeScan input.  We will doing this using a [perl script written by Santiago Sanchez-Ramirez](https://github.com/santiagosnchez/vcf2bayescan). 

```
#SBATCH lines

cd /sackettl/BTPD_7resistant/

perl vcf2bayescan/vcf2bayescan.pl -p popfile_new.txt -v btpd_SNP_QC_vcftools_complete.fltr_mean2sd.recode.vcf  

``` 
I included a population definition file here, using the same individual names as above. The file looks like this:

```
R_47A15080701_DSW64911-V   fatality
R_47A15080702_DSW64909-V   fatality
R_17A26060701_DSW64900-V   survivor
R_17a28060701_DSW64901-V   survivor
...
```
You have to make sure there are no extra lines in the pop file, or they will be read as an extra population with no individuals.


### 4. Detect outliers with OutFLANK
:warning: My outflank output produced ~1000 inferred outliers, but when I look at allele frequencies at those SNPs in the two groups, they were indistinguishable. I could not figure out what was going wrong and decided not to use this analysis.   :warning: 

[OutFLANK](https://github.com/whitlock/OutFLANK/) is an R package that calculates its own FST values because it does not allow negative values (which are generated by estimators that make adjustments for sample size), and it removes homozygous loci when creating the null distribution. This analysis uses an R script found [here](https://github.com/CassinSackett/PlagueResistance/blob/main/outflank_test.R).

After you have installed the OutFLANK package in R, run these steps to read in your file, convert it, and run OutFLANK:

```   
library(OutFLANK)
library(vcfR)
setwd("/Users/Loren/BTPD_resistant_not/")

btpd14vcf <- read.vcfR("btpd_SNP_QC_vcftools_complete.fltr_mean2sd.recode.vcf.gz")

bt14geno <- extract.gt(btpd14vcf) # Character matrix containing the genotypes
bt14position <- getPOS(btpd14vcf) # Positions in bp
bt14chromosome <- getCHROM(btpd14vcf) # Chromosome/scaffold information

#If you have any instances of two SNPs at the same position on different scaffolds, you need to link the scaffolds and position so that you can figure out which scaffolds your outliers are on:
SNPnames <- paste(bt13chromosome,bt13position, sep=":")


btmat <- matrix(NA, nrow = nrow(bt14geno), ncol = ncol(bt14geno))

sum(is.na(btmat)) #should be 0. If not, need to replace NAs with 9s so OutFLANK can read the data

newbtmat <- t(btmat) #transpose because we need SNPs in columns and individuals in rows

newbtmat[is.na(newbtmat)] <- 9 #change NAs to 9s
sum(is.na(newbtmat)) #verify this is 0 - worked correctly
levels(factor(as.numeric(unlist(newbtmat))))

newbtmat[bt14geno %in% c("0/0", "0|0")] <- 0
newbtmat[bt14geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
newbtmat[bt14geno %in% c("1/1", "1|1")] <- 2

table(as.vector(newbtmat))

popmembership <- c("fatality","fatality","survived","survived","fatality","antibodies","fatality","fatality","survived","survived","survived","fatality","fatality","survived")

bt_fst <- MakeDiploidFSTMat(newbtmat, locusNames = bt14position, popNames = popmembership)

#run OutFLANK
OutFLANK(FstDataFrame=bt_fst, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.05, NumberOfSamples=3, qthreshold=0.05)
```


### 5. Test for loci under selection in BayPass

After you have the BayeScan file for #3 above, use the python script [geste2baypass.py](https://github.com/CoBiG2/RAD_Tools/blob/master/geste2baypass.py) to convert this file to BayPass format.

```
python3 geste2baypass.py BayeScan/btpd_QC_complete-mean2sd_popnums.txt btpd_QC_complete_mean2sd_BP.txt 
```

Then, if you are using the default parameter values, you can run Baypass in just one line:

```
g_baypass -gfile btpd_QC_complete_mean2sd_BP.txt -outprefix btpd13_default2 -nthreads 44
```

The output will have a list of SNP names as 1 through Nsnps, so you'll need to pair them with their identities again. You can get a list of sites in the vcf file and unless you have manipulated the vcf file, the BayeScan input file or the BayPass input file, the order will be the same.
```
export PERL5LIB=/project/sackettl/vcftools/src/perl

/project/sackettl/vcftools/bin/vcf-validator btpd_QC_complete_mean2sd.recode.vcf  

```
and then
```
paste btpd_QC_complete_mean2sd.sites btpd_QC_complete_mean2sd_deflt_summary_pi_xtx.txt > btpd_QC_complete_mean2sd_deflt_pi_xtx_sites.txt
```




## III. Estimates of site-by-site heterozygosity
You may theoretically be able to get this in an R package, but the strategy of using brute force :weight_lifter:  in bash gives consistent, manipulable results. But mostly, I'm better at coding in bash than in R :grimacing: 

For large datasets, I use an interactive job instead of running on my laptop. You don't want to request too many resources or you'll be in the queue for a long time.

```
srun --time=2:00:00 --ntasks=8 --nodes=1 --account=loni_xx --partition=single --pty /bin/bash
```
where ```loni_xx``` is the name of the computing allocation. 


#### 1. Generate a file of positions and genotypes for all samples using vcf-query. GT is GenoType, and the brackets loop over all samples. This also works on .vcf.gz files.
```
/vcftools/vcf-query –f '%CHROM:%POS\t%REF[\t%GT]\n' BTPD_snponly.vcf > BTPD_genotypes.txt 
```

#### 2. Generate a list of individuals in the same order as the vcf file
```
vcf-query –l BTPD_snponly.vcf > BTPD-list.txt 
```

#### 3. We will use that list to create a header file that we will then concatenate the genotypes onto. 
First, insert 2 lines at the top of the BTPD-list.txt file with GHROM:POS on the first line and REF_ALELLE on the second line. Your new file will look like this:
```
CHROM:POS
REF_ALLELE
indiv1
indiv2
...
```

Next, transpose (i.e., replace new lines with tabs) to create a tab-delimited header file 
```
tr '\n' '\t' < BTPD-list.txt > BTPD-header.txt
```
:::danger
Using a text editor, add a new line (\n) to the end of this line so that concatenation will occur correctly.
:::

Now you can concatenate the genotypes file onto the header file so that you know which individual each genotype corresponds to:
```
cat BTPD-header.txt BTPD_genotypes.txt > BTPD_genotypes_withheader.txt
```

:::info
One of my datasets for some reason had a mix of genotypes separated by the classic '/' character, while others were separated by the '|'. I fixed this by using sed to replace the pipes with the slashes, and learned a fun fact: sed's substitute can use any character as a delimiter; it doesn't have to be the classic ```'s/old/new/g' ```. I used an exclamation point:
```
sed 's!|!/!g' BTPD_genotypes_withheader.txt > BTPD_genotypes_withheader_edit.txt
```
:::


#### 4.  Generate files of all heterozygous and homozygous sites for each individual

We will do this by looping through all samples using an awk script, naming the files according to the column/header name. First, create the awk script (this was written by Taylor Callicrate): 
```
BEGIN{
        FS="\t"
}
NR==1 {
        name=$x;
        print (name) > (name ".homs");
        next
}
{
split($x, a, "/"); if (a[1] == a[2]) print ($1 "\t" $2 "\t" $x) >> (name ".homs")
        }
END {
}

```
and then loop through all your samples--starting in column 3 until the end of the file, usually--to execute this script:
```
for i in {3..146}; do awk -f loren2.awk -v x=$i BTPD_N144_gtypes.txt; done
```


#### 5. Generate a file from your reference genome of windows of a certain size

First, we'll need to extract the scaffold lengths from our reference genome.
```

```


Next, we can make windows corresponding to each scaffold. This allows us to use, for instance, 10000 base pair windows without overshooting beyond the end of a chromosome. Using scaffold lengths as an input will cause the last window on each scaffold to be shorter than 10000bp according to each scaffold's true length.

```
bedtools2/bin/bedtools makewindows -g ../BTPD_firstgenome/scaffold_lengths.txt -w 100000 > btpdgenome_100kwind.bed
```

#### 6. Use coverage bed to examine the coverage of both heterozygotes and homozygotes in each window
```
```
Once you have the first .hets and .homs files after loren2.awk, remove the header since the filename already contains the individual name. 







## IV. Long Runs of Homozygosity (LROH)
Typically, LROH are used to detect inbreeding in certain individuals. However, if multiple outbred individuals exhibit LROH in the same region, this can be indicative of a selective sweep.

First, get estimates of inbreeding coefficients for each individual in vcftools. This needs to be done on a per-scaffold basis, so if you have tons of scaffolds, you may want to consider using only the longest ones, or only for the scaffolds that have candidate SNPs from the other analyses.

```
for i in {1..147}; do ../../vcftools/src/cpp/vcftools --vcf btpd_complete.fltr_mean2sd.recode.vcf --LROH --chr scaffold$i --out btpd_complete.fltr_mean2sd_$i; done
```




## V. Calculate LD in regions of interest
I do this in plink (in vcftools I was not getting the pairwise estimates I wanted). This syntax will output results in a table:


```
~/plink_mac_20210606/plink --file ../plink_analyses/btpd_SNP_QC_vcftools_complete.fltr_mean2sd --allow-extra-chr --allow-no-sex --chr scaffold214 --r2 --out btpd_SNP_QC_vcftools_complete.fltr_mean2sd_scaff214LD
```



## VI. Follow up on candidate loci
### 6a. Allele frequency differences
Since I am interested in adaptation to a novel selection pressure (evolved resistance to an introduced pathogen), I wanted to look at all of my candidate SNPs that had a different allele in survivors from the allele in the reference genome (which I am assuming to be a non-resistance allele).

First, get your list of candidate sites. Mine looked something like this:
```
scaffold2693	50281
scaffold2693	50354
scaffold2693	50406
...
```
I am going to make a .vcf file with just my candidate sites (e.g., to compare a PCA with all sites and candidate sites):

```
~/vcftools/src/cpp/vcftools --gzvcf btpd_SNP_QC_complete.fltr_mean2sd.recode.vcf.gz --positions moderate-candidates_1334sites.txt --out moderate-1334candidates --recode
```

I used the .012 output from vcftools to count the number of alternate alleles in survivors at each SNP. To do so:

Get list of positions .012.pos and put them all in one column ```tr "\t" ":" < .012.pos > .positions```

Transpose columns into a row to create a header file ```tr "\n" "\t" < .positions > genotypes.header```

and then opened the header in BBEdit to add a new line to the end of the file (I need to figure out how to do this in bash!)

Then I separated my .012 file into two files with the two groups I am comparing (survivors vs fatalities). In each file, I summed the number of non-reference alleles at each SNP across all survivors:

```
awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print "for column "i" is " sum[i];}' btpd_moderate1300candidates.genotypes.txt > btpd_moderate1300candidates.num_nonref-alleles.txt
```
and then concatenated that on the end of the genotypes file

```
cat header genotypes alternate_allelecount > btpd_7survivors_gtypes_altallelecounts.txt
```

Next, I want to sort the loci in order of the MOST alternate alleles in the survivors group. There are ways to do this directly in bash, but I found it easier to transpose the file and then sort by columns.

To transpose, I used the following script which I found [here](https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash).

```
BEGIN { FS=OFS="\t" }
{
    for (rowNr=1;rowNr<=NF;rowNr++) {
        cell[rowNr,NR] = $rowNr
    }
    maxRows = (NF > maxRows ? NF : maxRows)
    maxCols = NR
}
END {
    for (rowNr=1;rowNr<=maxRows;rowNr++) {
        for (colNr=1;colNr<=maxCols;colNr++) {
            printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
        }
    }
}
```
and then you call that script with

``` 
awk -f transpose_file.awk btpd_7survivors_gtypes_altallelecounts.txt > btpd_7survivors_gtypes_altallelecounts_transpose.txt 
```

Now we have two files that looks like this:

```
indiv_name	R_17a29060703_antibody	R_17A26060701_survivor	R_17a28060701_survivor	R_19A24070805_survivor	R_19A26070701_survivor	R_47A13080701_survivor	R_47A15080703_survivor	number_alternate_alleles
scaffold3:3018055	2	2	2	2	1	2	2	13
scaffold3:3019133	2	2	2	2	1	2	2	13
```

I want to calculate the difference in number of alternate alleles between my two groups, so I first pasted the files from each group next to each other:

```
paste btpd_7survivors_gtypes_alleletallies_transpose.txt btpd_7fatalities_gtypes_alleletallies_transpose.txt > btpd14_gtypes_alt-alleletallies.txt
```
You could cut the duplicate scaffold:pos field, or keep it as a check that your loci are in the same order in the two groups.

Finally, you can use awk to calculate the absolute value of the difference in allele count between groups, where the $9 refers to column 9 (the first column with my allele tally):
```
awk 'BEGIN{FS=OFS="\t"} {print $0, ($9>$18 ? $9-$18 : $18-$9)}' btpd14_gtypes_alt-alleletallies.txt > btpd14_gtypes_alt-allele-diffs.txt
```
Now you can follow up on these loci, such as sorting by magnitude of differences, searching for fixed differences, or finding alleles that are present only in one group (e.g., survivors)

:eyes: A note: If you suspect adaptation was conferred by a recent mutation, or if you expect dominance of a mutation, you may also want to look at loci that have one alternate allele per individual in the "survived" group that is absent in the "fatalities" group. In these data, there would be 7 alternate alleles in the survivors group and 0 in the fatalities.


### 6b. Inferring gene location of candidate SNPs
Right now, we have a list of candidate SNPs on different scaffolds. But unless the genome is very well annotated, we may want to blast the genomic region containing our SNP to figure out what gene or genomic region the SNP is in. Suppose our list ```candidates-list.bed``` looks something like this:

```
scaffold12    499208
scaffold496    292017
```

First, you need to figure out what size region you want to blast, and make a bed file with the start and end position. Here, we'll blast 500 bp upstream (note the '-' prior to the equal sign in var1) and 1000 bp downstream (note the '-' prior to the equal sign in var1)

```
cat candidates-list.bed | awk '{var1=$2; var1-=500; var2=$2; var2+=1000; print $1 "\t" var1 "\t" var2}' > candidates_500bpup1kbdown.bed
```

Now we can use [bedtools](https://github.com/arq5x/bedtools2) to extract that region from the reference genome so we have something to blast:

```
bedtools getfasta -fi HcZfUnix_ref/HcZfUnix.fasta -bed candidates_500bpup1kbdown.bed -fo candidates_500bpup1kbdown.fa
```
where '-fi' refers to the reference genome fasta, and '-fo' refers to the name of the fasta sequence you are generating.

Finally, you can blast against a custom or existing database: 

```
/ncbi-blast-2.2.30+/bin/blastn -query ExpInf_400bflanks.fa -out ExpInf_400bflanks_blastout.txt -db nr -remote -evalue 1e-3 -outfmt 7 
```

```-outfmt 7``` puts the output in a nice table w/ comment lines; ```-outfmt 11``` is _complete_ and can be later converted to other formats with an NCBI (?) script blast_formatter.). If there is no outfmt flag, you can add ```-line_length 80``` (or whatever you want) to lengthen the description. ```-db``` asks for the database name: nr is the non-redundant protein sequences from several sources; nt is the partially non-redundant nucleotide sequences.





