DADA2_FALL21_EEMB144L
================
Hope Hahn
11/22/2021

# Load packages

``` r
library(tidyverse)
library(dada2)
library(ShortRead)
```

# Import data

``` r
#save the path to the directory with a COPY of your unzipped fastq files that you WILL work with. MAKE SURE YOU HAVE ANOTHER DIRECTORY WITH THE FILES THAT YOU WILL NEVER DIRECTLY WORK WITH. 

path <- "~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq" #make sure there is no / at the end of the path
#also make sure there are no unzipped files in this directory

#store the names of fwd and rv files as lists
fnFs <- list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)
fnRs <- list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)
```

# Retrieve orientation of primers

``` r
FWD = "GTGYCAGCMGCCGCGGTAA"
REV = "GGACTACNVGGGTWTCTAAT"

#now store all orientations of fwd and rvs primers 
allOrients <- function(primer) {
  # Biostrings works w/ DNAString objects rather than character vectors
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
              RevComp = reverseComplement(dna))
  #Convert back to character vector
  return(sapply(orients, toString))
  }
  
#Store the fwd and rvs orientations separately 
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#view the orientations of the primers
FWD.orients
```

    ##               Forward            Complement               Reverse 
    ## "GTGYCAGCMGCCGCGGTAA" "CACRGTCGKCGGCGCCATT" "AATGGCGCCGMCGACYGTG" 
    ##               RevComp 
    ## "TTACCGCGGCKGCTGRCAC"

``` r
REV.orients
```

    ##                Forward             Complement                Reverse 
    ## "GGACTACNVGGGTWTCTAAT" "CCTGATGNBCCCAWAGATTA" "TAATCTWTGGGVNCATCAGG" 
    ##                RevComp 
    ## "ATTAGAWACCCBNGTAGTCC"

#search for Primers

``` r
primerHits <- function(primer, fn) {
  #Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits >0))
  }
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```

    ##                  Forward Complement Reverse RevComp
    ## FWD.ForwardReads       0          0       0       0
    ## FWD.ReverseReads       0          0       0       1
    ## REV.ForwardReads       0          0       0       1
    ## REV.ReverseReads       0          0       0       0

# Insepct read quality profiles

You should look at least some of the quality profiles to assess the
quality of the sequencing run.

## Forward reads

``` r
plotQualityProfile(fnFs[1:12])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](DADA2_FALL21_EEMB144L_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Reverse reads

``` r
plotQualityProfile(fnRs[1:12])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](DADA2_FALL21_EEMB144L_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Filtering and Trimming

``` r
#Get the sample names
#define the basename of the fnFs as the first part of each fastq file name until "_L"
#apply this to all samples
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
sample.names
```

    ##  [1] "144_A0_S1"    "144_A4_S2"    "144_A8_S3"    "144_B0_S4"    "144_B4_S5"   
    ##  [6] "144_B8_S6"    "144_C0_S7"    "144_C4_2_S33" "144_D8_2_S34" "144_E0_2_S35"
    ## [11] "144_G0_S24"   "144_G4_S25"   "144_G8_S26"   "144_H0_S27"   "144_H4_S28"  
    ## [16] "144_H8_S29"

``` r
#created a "filtered" folder in the working directory as a place to put all the new filtered fastQ files. 
filt_path <- file.path(path, "filtered")
#add the appropriate designation string to any new files made that will be put into the "filtered" folder
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 150), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)
#look at this output. This tells you how many reads were removed. 
out
```

    ##                                   reads.in reads.out
    ## 144_A0_S1_L001_R1_001.fastq.gz       10955      9112
    ## 144_A4_S2_L001_R1_001.fastq.gz       10502      8124
    ## 144_A8_S3_L001_R1_001.fastq.gz        4506      3220
    ## 144_B0_S4_L001_R1_001.fastq.gz       17862     13583
    ## 144_B4_S5_L001_R1_001.fastq.gz       12206      9589
    ## 144_B8_S6_L001_R1_001.fastq.gz       12477     10322
    ## 144_C0_S7_L001_R1_001.fastq.gz        1697      1166
    ## 144_C4_2_S33_L001_R1_001.fastq.gz   326236    311641
    ## 144_D8_2_S34_L001_R1_001.fastq.gz   298669    282860
    ## 144_E0_2_S35_L001_R1_001.fastq.gz   329028    314579
    ## 144_G0_S24_L001_R1_001.fastq.gz      40935     36168
    ## 144_G4_S25_L001_R1_001.fastq.gz      40109     35236
    ## 144_G8_S26_L001_R1_001.fastq.gz      35610     31788
    ## 144_H0_S27_L001_R1_001.fastq.gz      63711     57388
    ## 144_H4_S28_L001_R1_001.fastq.gz      27892     24291
    ## 144_H8_S29_L001_R1_001.fastq.gz      36860     32338

# Learn the error rates

``` r
errF <- learnErrors(filtFs, multithread = TRUE)
```

    ## 155908080 total bases in 649617 reads from 9 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread = TRUE)
```

    ## 144629400 total bases in 964196 reads from 10 samples will be used for learning the error rates.

<img src="DADA2_FALL21_EEMB144L_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

<img src="DADA2_FALL21_EEMB144L_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

# Dereplication

``` r
derepFs <- derepFastq(filtFs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A0_S1_F_filt.fastq

    ## Encountered 6849 unique sequences from 9112 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A4_S2_F_filt.fastq

    ## Encountered 6442 unique sequences from 8124 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A8_S3_F_filt.fastq

    ## Encountered 2718 unique sequences from 3220 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B0_S4_F_filt.fastq

    ## Encountered 10336 unique sequences from 13583 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B4_S5_F_filt.fastq

    ## Encountered 6917 unique sequences from 9589 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B8_S6_F_filt.fastq

    ## Encountered 6872 unique sequences from 10322 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_C0_S7_F_filt.fastq

    ## Encountered 885 unique sequences from 1166 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_C4_2_S33_F_filt.fastq

    ## Encountered 75469 unique sequences from 311641 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_D8_2_S34_F_filt.fastq

    ## Encountered 66946 unique sequences from 282860 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_E0_2_S35_F_filt.fastq

    ## Encountered 71204 unique sequences from 314579 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G0_S24_F_filt.fastq

    ## Encountered 10822 unique sequences from 36168 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G4_S25_F_filt.fastq

    ## Encountered 8483 unique sequences from 35236 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G8_S26_F_filt.fastq

    ## Encountered 8503 unique sequences from 31788 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H0_S27_F_filt.fastq

    ## Encountered 15044 unique sequences from 57388 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H4_S28_F_filt.fastq

    ## Encountered 5919 unique sequences from 24291 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H8_S29_F_filt.fastq

    ## Encountered 9702 unique sequences from 32338 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A0_S1_R_filt.fastq

    ## Encountered 7544 unique sequences from 9112 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A4_S2_R_filt.fastq

    ## Encountered 6859 unique sequences from 8124 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_A8_S3_R_filt.fastq

    ## Encountered 2941 unique sequences from 3220 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B0_S4_R_filt.fastq

    ## Encountered 11786 unique sequences from 13583 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B4_S5_R_filt.fastq

    ## Encountered 7882 unique sequences from 9589 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_B8_S6_R_filt.fastq

    ## Encountered 8371 unique sequences from 10322 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_C0_S7_R_filt.fastq

    ## Encountered 882 unique sequences from 1166 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_C4_2_S33_R_filt.fastq

    ## Encountered 83078 unique sequences from 311641 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_D8_2_S34_R_filt.fastq

    ## Encountered 70371 unique sequences from 282860 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_E0_2_S35_R_filt.fastq

    ## Encountered 74371 unique sequences from 314579 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G0_S24_R_filt.fastq

    ## Encountered 13092 unique sequences from 36168 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G4_S25_R_filt.fastq

    ## Encountered 11615 unique sequences from 35236 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_G8_S26_R_filt.fastq

    ## Encountered 10202 unique sequences from 31788 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H0_S27_R_filt.fastq

    ## Encountered 19510 unique sequences from 57388 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H4_S28_R_filt.fastq

    ## Encountered 9788 unique sequences from 24291 total sequences read.

    ## Dereplicating sequence entries in Fastq file: ~/Desktop/github/144l_students_2021/Input_Data/week9/Kelp_Remin_fastq/filtered/144_H8_S29_R_filt.fastq

    ## Encountered 11301 unique sequences from 32338 total sequences read.

``` r
# Name the derep-class objects by the sample names 
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

# Infer sequence variants

``` r
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
```

    ## Sample 1 - 9112 reads in 6849 unique sequences.
    ## Sample 2 - 8124 reads in 6442 unique sequences.
    ## Sample 3 - 3220 reads in 2718 unique sequences.
    ## Sample 4 - 13583 reads in 10336 unique sequences.
    ## Sample 5 - 9589 reads in 6917 unique sequences.
    ## Sample 6 - 10322 reads in 6872 unique sequences.
    ## Sample 7 - 1166 reads in 885 unique sequences.
    ## Sample 8 - 311641 reads in 75469 unique sequences.
    ## Sample 9 - 282860 reads in 66946 unique sequences.
    ## Sample 10 - 314579 reads in 71204 unique sequences.
    ## Sample 11 - 36168 reads in 10822 unique sequences.
    ## Sample 12 - 35236 reads in 8483 unique sequences.
    ## Sample 13 - 31788 reads in 8503 unique sequences.
    ## Sample 14 - 57388 reads in 15044 unique sequences.
    ## Sample 15 - 24291 reads in 5919 unique sequences.
    ## Sample 16 - 32338 reads in 9702 unique sequences.

``` r
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
```

    ## Sample 1 - 9112 reads in 7544 unique sequences.
    ## Sample 2 - 8124 reads in 6859 unique sequences.
    ## Sample 3 - 3220 reads in 2941 unique sequences.
    ## Sample 4 - 13583 reads in 11786 unique sequences.
    ## Sample 5 - 9589 reads in 7882 unique sequences.
    ## Sample 6 - 10322 reads in 8371 unique sequences.
    ## Sample 7 - 1166 reads in 882 unique sequences.
    ## Sample 8 - 311641 reads in 83078 unique sequences.
    ## Sample 9 - 282860 reads in 70371 unique sequences.
    ## Sample 10 - 314579 reads in 74371 unique sequences.
    ## Sample 11 - 36168 reads in 13092 unique sequences.
    ## Sample 12 - 35236 reads in 11615 unique sequences.
    ## Sample 13 - 31788 reads in 10202 unique sequences.
    ## Sample 14 - 57388 reads in 19510 unique sequences.
    ## Sample 15 - 24291 reads in 9788 unique sequences.
    ## Sample 16 - 32338 reads in 11301 unique sequences.

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = T)
```

    ## 7328 paired-reads (in 22 unique pairings) successfully merged out of 8217 (in 225 pairings) input.

    ## 6547 paired-reads (in 22 unique pairings) successfully merged out of 7175 (in 219 pairings) input.

    ## 1744 paired-reads (in 10 unique pairings) successfully merged out of 2471 (in 126 pairings) input.

    ## 8482 paired-reads (in 65 unique pairings) successfully merged out of 11735 (in 1115 pairings) input.

    ## 6666 paired-reads (in 59 unique pairings) successfully merged out of 8503 (in 668 pairings) input.

    ## 7145 paired-reads (in 51 unique pairings) successfully merged out of 9263 (in 743 pairings) input.

    ## 530 paired-reads (in 37 unique pairings) successfully merged out of 813 (in 93 pairings) input.

    ## 291343 paired-reads (in 547 unique pairings) successfully merged out of 308334 (in 2587 pairings) input.

    ## 267782 paired-reads (in 400 unique pairings) successfully merged out of 281155 (in 1558 pairings) input.

    ## 292389 paired-reads (in 509 unique pairings) successfully merged out of 312279 (in 2010 pairings) input.

    ## 30153 paired-reads (in 189 unique pairings) successfully merged out of 35599 (in 405 pairings) input.

    ## 32819 paired-reads (in 124 unique pairings) successfully merged out of 34822 (in 339 pairings) input.

    ## 29087 paired-reads (in 142 unique pairings) successfully merged out of 31490 (in 296 pairings) input.

    ## 40598 paired-reads (in 238 unique pairings) successfully merged out of 56677 (in 597 pairings) input.

    ## 22804 paired-reads (in 108 unique pairings) successfully merged out of 23905 (in 277 pairings) input.

    ## 28274 paired-reads (in 147 unique pairings) successfully merged out of 31973 (in 303 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                        sequence
    ## 1 TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGG
    ## 2 TACGGAGGGTCCGAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGCAGGCGGTCAATTAAGTCAGAGGTGAAATCCCATAGCTCAACTATGGAACTGCCTTTGATACTGGTTGACTTGAGTCATATGGAAGTAGATAGAATGTGTAGTGTAGCGGTGAAATGCATAGATATTACACAGAATACCGATTGCGAAGGCAGTCTACTACGTATGTACTGACGCTGAGGGACGAAAGCGTGGGGAGCGAACAGG
    ## 3 TACGGAGGATCCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGCGGTCTAATAAGTCAGAGGTGAAATCCTACAGCTCAACTGTAGCATTGCCTTTGATACTGTTAGACTTGAGTTATTGTGAAGTAGTTAGAATGTGTAGTGTAGCGGTGAAATGCATAGATATTACACAGAATACCGATTGCGAAGGCAGATTACTAACAATATACTGACGCTGAGGGACGAAAGCGTGGGTAGCGAACGGG
    ## 4 TACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGTACGTAGGCGGATTAATAAGTTAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTAAAACTGTTAGTCTTGAGATCGAGAGAGGTGAGTGGAATTCCAAGTGTAGAGGTGAAATTCGTAGATATTTGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTACGAAAGTGTGGGGAGCAAACAGG
    ## 5 CACGGAAGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAACTCACTAGAGTACGATAGAGGGAGGTAGAATTCATAGTGTAGCGGTGGAATGCGTAGATATTATGAAGAATACCAGTGGCGAAGGCGGCCTCCTGGATCTGCACTGACACTGAGGTGCGAAAGCGTGGGTAGCGAACAGG
    ## 6 TACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATTAGTAAGTTAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTAATACTGCTAGTCTTGAGTTCGAGAGAGGTAAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      1233       1       2    137         0      0      1   TRUE
    ## 2      1204       3       3    137         0      0      1   TRUE
    ## 3      1142       2       1    137         0      0      1   TRUE
    ## 4       768       4       4    137         0      0      1   TRUE
    ## 5       430       5       5    137         0      0      1   TRUE
    ## 6       323       7       8    137         0      0      1   TRUE

``` r
saveRDS(mergers, "~/Desktop/github/144l_students_2021/Output_Data/week9/dada_merged.rds")
saveRDS(mergers, "~/Desktop/github//144l_students_2021/Input_Data/week9/dada_merged.rds")
```

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # samples by unique sequence
```

    ## [1]   16 1099

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  251  252  253  254  255  256  257  258  265  266  270  280 
    ##    1   16 1009   51    6    6    4    1    1    1    2    1

# Remove the Chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
```

    ## Identified 221 bimeras out of 1099 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  16 878

check the proportion of sequences that are not chimeras

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9877861

# Assign taxonomy using a reference database

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/github/144l_students_2021/Input_Data/week9/Reference_Database/silva_nr_v138_train_set.fa", multithread = TRUE)
```

``` r
saveRDS(t(seqtab.nochim), "~/Desktop/github/144l_students_2021/Output_Data/week9/seqtab-nochimtaxa.rds")
saveRDS(taxa, "~/Desktop/github/144l_students_2021/Output_Data/week9/taxa.rds")


saveRDS(t(seqtab.nochim), "~/Desktop/github/144l_students_2021/Input_Data/week9/seqtab-nochimtaxa.rds")
saveRDS(taxa, "~/Desktop/github/144l_students_2021/Input_Data/week9/taxa.rds")
```
