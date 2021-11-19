Phyloseq
================
Hope Hahn
11/10/2021

# Intro

We explore the processed ACIDD 16S sequences using
[phyloseq](https://joey711.github.io/phyloseq/)

# Install phyloseq

``` r
#BiocManager::install("phyloseq")
```

``` r
library(tidyverse) 
library(phyloseq)
library(RColorBrewer)
library(readxl)
library(lubridate)
```

# Import Data

``` r
count.tab <- read_rds("~/Desktop/github/144l_students_2021/Input_Data/week7/seqtab-nochimtaxa.rds") #table of counts for each sequence in each sample
tax.tab <- read_rds("~/Desktop/github/144l_students_2021/Input_Data/week7/taxa.rds") #table that matches ASV to sequence

#you will need to download the ACIDD_Exp_Processed_DOC_BGE.rds and the ACIDD_Exp_BactAbund.xlsx files from Justine's Github repository (located in Input_Data/week7/) and put them into your own local Input Data folder
metadata_OG <- read_excel("~/Desktop/github/144l_students_2021/Input_Data/week7/ACIDD_Exp_BactAbund.xlsx", sheet = "Metadata")
glimpse(metadata_OG)
```

    ## Rows: 84
    ## Columns: 18
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <chr> "2017-12-16T21:30", "2017-12-17T10:00", "2017-…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, …
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALS…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171-…

``` r
metadata <- metadata_OG %>%
  mutate(Datetime = ymd_hm(Datetime))
glimpse(metadata)
```

    ## Rows: 84
    ## Columns: 18
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <dttm> 2017-12-16 21:30:00, 2017-12-17 10:00:00, 201…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, …
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALS…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171-…

``` r
metadata_doc <- read_rds("~/Desktop/github/144l_students_2021/Input_Data/week7/ACIDD_Exp_Processed_DOC_BGE.rds")
glimpse(metadata_doc)
```

    ## Rows: 84
    ## Columns: 64
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <dttm> 2017-12-16 21:30:00, 2017-12-17 10:00:00, 201…
    ## $ hours                   <dbl> 0.0, 12.5, 25.5, 48.0, 71.5, 95.0, 118.5, 122.…
    ## $ days                    <dbl> 0.0000000, 0.5208333, 1.0625000, 2.0000000, 2.…
    ## $ TOC                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_TOC                  <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ PTOC                    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_PTOC                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ cells                   <dbl> 1.30e+08, 1.34e+08, 1.28e+08, 1.55e+08, 1.55e+…
    ## $ sd_cells                <dbl> 20900000, 27600000, 22200000, 25200000, 319000…
    ## $ ln_cells                <dbl> 18.68305, 18.71335, 18.66754, 18.85894, 18.858…
    ## $ diff_ln_cells           <dbl> 0.000000000, 0.030305349, -0.045809536, 0.1913…
    ## $ bc                      <dbl> 0.32500, 0.33500, 0.32000, 0.38750, 0.38750, 0…
    ## $ ave_bc                  <dbl> 0.308750, 0.330000, 0.317500, 0.358750, 0.4100…
    ## $ sd_bc                   <dbl> 0.022980970, 0.007071068, 0.003535534, 0.04065…
    ## $ exp_start               <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5…
    ## $ exp_end                 <dbl> 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6…
    ## $ ln_cells_exp_start      <dbl> 18.85894, 18.85894, 18.85894, 18.85894, 18.858…
    ## $ ln_cells_exp_end        <dbl> 19.74776, 19.74776, 19.74776, 19.74776, 19.747…
    ## $ cells_exp_start         <dbl> 1.55e+08, 1.55e+08, 1.55e+08, 1.55e+08, 1.55e+…
    ## $ cells_exp_end           <dbl> 3.77e+08, 3.77e+08, 3.77e+08, 3.77e+08, 3.77e+…
    ## $ days_exp_start          <dbl> 2.979167, 2.979167, 2.979167, 2.979167, 2.9791…
    ## $ days_exp_end            <dbl> 4.937500, 4.937500, 4.937500, 4.937500, 4.9375…
    ## $ mew                     <dbl> 0.4538656, 0.4538656, 0.4538656, 0.4538656, 0.…
    ## $ doubling                <dbl> 1.5272081, 1.5272081, 1.5272081, 1.5272081, 1.…
    ## $ delta_cells             <dbl> 247000000, 247000000, 247000000, 247000000, 24…
    ## $ delta_bc                <dbl> 0.61750, 0.61750, 0.61750, 0.61750, 0.61750, 0…
    ## $ ave_mew                 <dbl> 0.5441279, 0.5441279, 0.5441279, 0.5441279, 0.…
    ## $ sd_mew                  <dbl> 0.09366960, 0.09366960, 0.09366960, 0.09366960…
    ## $ ave_doubling            <dbl> 1.3099139, 1.3099139, 1.3099139, 1.3099139, 1.…
    ## $ sd_doubling             <dbl> 0.22549685, 0.22549685, 0.22549685, 0.22549685…
    ## $ ave_delta_cells         <dbl> 232500000, 232500000, 232500000, 232500000, 23…
    ## $ sd_delta_cells          <dbl> 15047361, 15047361, 15047361, 15047361, 150473…
    ## $ ave_delta_bc            <dbl> 0.581250, 0.581250, 0.581250, 0.581250, 0.5812…
    ## $ sd_delta_bc             <dbl> 0.037618403, 0.037618403, 0.037618403, 0.03761…
    ## $ ave_lag                 <dbl> 3.46875, 3.46875, 3.46875, 3.46875, 3.46875, 3…
    ## $ sd_lag                  <dbl> 0.5080646, 0.5080646, 0.5080646, 0.5080646, 0.…
    ## $ interp_toc              <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ interp_bc               <dbl> 0.3, 0.3, 0.3, 0.4, 0.4, 0.5, 0.9, NA, NA, NA,…
    ## $ doc                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ bioav_doc               <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ doc_exp_end             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ delta_doc               <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ tdelta_doc              <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ bge                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_toc                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_toc                  <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_bioav_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_bioav_doc            <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_delta_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_delta_doc            <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_tdelta_doc          <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_tdelta_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_bge                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_bge                  <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…

``` r
join <- left_join(metadata_doc, metadata)
sample.tab <- join %>%
  drop_na(DNA_SampleID) %>% 
  column_to_rownames(var = "DNA_SampleID") 
```

# Phyloseq Object

We need to create a phyloseq object that merges all three datasets.
Sometimes this doesn’t work beacuse of the format of the data files.
Make sure all the sample names between the sampleinfo.txt and
seqtab-nochimtaxa.txt are the same

``` r
OTU = otu_table(count.tab, taxa_are_rows = TRUE) 
TAX = tax_table(tax.tab)
SAM = sample_data(sample.tab)
ps = phyloseq(OTU,TAX,SAM) 
```

# Filter sequences

We will filter out chloroplasts and mitochondria, because we only
intended to amplify bacterial sequences. It’s good to check you don’t
have anything lurking in the taxonomy table.

``` r
sub_ps <- ps %>%
  # subset_samples(Experiment == "ASH172") %>%  #use this function if you want to only include some subset of your sample set in the subsequent analysis
  subset_taxa(Family  != "mitochondria" & Order  != "Chloroplast")
```

Q1: What did we do in the code chunk above? (What is the difference
between phyloseq object ps and phyloseq object sub_ps?) Why did we do
this?

A1: Because we only want to look at bacterial sequences, we have to take
out the observations that are mitochondria or chloroplasts.We use the
function “subset_taxa” to remove the sequences that belong to
mitochomdria and chloroplasts. sub_ps is derived from ps, but has the
chloroplasts and mitochondria removed.

# Sample Summary

As a first analysis, we will look at the distribution of read counts
from our samples

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />
Q2: Describe what “sequencing depth” means in your own words.

A2: Sequencing depth describes the amount of reads that were used to
reconstruct the genome.

``` r
# mean, max and min of sample read counts
summary(sample_sum_df)
```

    ##       sum      
    ##  Min.   :6047  
    ##  1st Qu.:7153  
    ##  Median :7695  
    ##  Mean   :7695  
    ##  3rd Qu.:8453  
    ##  Max.   :9250

# Beta Diversity

Beta diversity involves calculating metrics such as distances or
dissimilarities based on pairwise comparisons of samples – they don’t
exist for a single sample, but rather only as metrics that relate
samples to each other. i.e. beta diversity = patterns in community
structure between samples

Since differences in sampling depths between samples can influence
distance/dissimilarity metrics, we first need to somehow normalize the
read depth across our samples.

## Subsample

We will rarefy (random subsample with replacement) the read depth of the
samples first (scale to the smallest library size).

[Case for not
subsampling](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

[Response blog for
subsampling](https://www.polarmicrobes.org/how-i-learned-to-stop-worrying-and-love-subsampling-rarifying/)

Read depth is an artefact of a machine made by a company in San Diego,
not anything about your samples or their biology. It is totally
artifactual, and controlling for artifacts is critical in science.
Subsampling randomly is the simplest way to control for this, and the
question is whether this is the “best” way of controlling for it. See
links above for alternative arguments about what the best way of
controlling for this artefact is.

A strong reason to subsample is to standardize effort. The bottom line
is that in all experimental design you should not be comparing things to
which you devote different effort in resolution. For instance, you don’t
sample one site once a week and another once a month if you want to
compare the dynamics between the sites. You standardize effort.

With that said, the bigger your differential in mean (or median) read
depth (reads/sample) between pre- and post-subsampling, the greater the
“effect” on beta diversity.

Examples:

-   means reads before = 40k, mean reads after = 1k, big effect.
-   mean reads before = 40k, mean reads after = 20k, small effect.
-   mean reads before = 2k, mean reads after = 1k, small effect.

We will subsample to the minimum read depth of all samples and not
subsample. We’ll then compare the mean reads pre- and post-subsampling
and also compare beta diversity patterns

``` r
ps_min <-  rarefy_even_depth(sub_ps, sample.size = min(sample_sums(sub_ps)))
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 2OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
mean(sample_sums(sub_ps)) #7686
```

    ## [1] 7695.25

``` r
mean(sample_sums(ps_min)) #6048 this is also the same as min(sample_sums(sub)ps) 
```

    ## [1] 6047

Q3: Do you think that the subsampling we did here will have a large
effect on our beta diveristy analyses? why or why not?

A3: It will likely have a small effect because the mean is 7695 reads
while the minimum is 6047. There is not a huge differential between
these two numbers, so we do not expect a large effect of beta diversity
analysis.

## NMDS

One of the best exploratory analyses for amplicon data is unconstrained
ordinations. Here we will look at non-metric multidimensional scaling
(NMDS) ordinations of our full community samples. For NMDS plots it’s
important to set a seed since the starting positions of samples in the
algorithm is random.

``` r
set.seed(1)
# Ordinate
nmds <- ordinate(sub_ps, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.04055426 
    ## Run 1 stress 0.1688387 
    ## Run 2 stress 0.04055427 
    ## ... Procrustes: rmse 1.41369e-05  max resid 2.390827e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.04055426 
    ## ... Procrustes: rmse 7.134317e-06  max resid 1.225786e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.3118493 
    ## Run 5 stress 0.04055426 
    ## ... Procrustes: rmse 3.474084e-06  max resid 5.620063e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.04055426 
    ## ... Procrustes: rmse 5.460043e-06  max resid 9.278492e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.04055427 
    ## ... Procrustes: rmse 1.594292e-05  max resid 2.654566e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.04055426 
    ## ... New best solution
    ## ... Procrustes: rmse 3.651678e-06  max resid 6.778531e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.04055426 
    ## ... Procrustes: rmse 3.221665e-06  max resid 5.495844e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.04055426 
    ## ... Procrustes: rmse 3.161775e-06  max resid 5.890661e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.04055427 
    ## ... Procrustes: rmse 1.348029e-05  max resid 2.200048e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.1745143 
    ## Run 13 stress 0.04055426 
    ## ... Procrustes: rmse 2.16648e-06  max resid 4.307359e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.04055427 
    ## ... Procrustes: rmse 1.648896e-05  max resid 2.84659e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.04055426 
    ## ... Procrustes: rmse 2.901513e-06  max resid 4.809554e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.04055426 
    ## ... Procrustes: rmse 6.237738e-06  max resid 1.029078e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.04055427 
    ## ... Procrustes: rmse 2.436672e-05  max resid 4.198842e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.04055426 
    ## ... Procrustes: rmse 2.140739e-06  max resid 3.796994e-06 
    ## ... Similar to previous best
    ## Run 19 stress 0.04055426 
    ## ... Procrustes: rmse 1.700096e-06  max resid 3.96457e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.04055426 
    ## ... Procrustes: rmse 2.698669e-06  max resid 5.048448e-06 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
set.seed(1)
# Ordinate
nmds_min <- ordinate(ps_min, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.04483483 
    ## Run 1 stress 0.174781 
    ## Run 2 stress 0.04483483 
    ## ... New best solution
    ## ... Procrustes: rmse 1.468714e-06  max resid 2.939261e-06 
    ## ... Similar to previous best
    ## Run 3 stress 0.04483483 
    ## ... Procrustes: rmse 2.036675e-06  max resid 4.282906e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.3197986 
    ## Run 5 stress 0.04483483 
    ## ... Procrustes: rmse 9.91949e-07  max resid 2.089411e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.04483483 
    ## ... Procrustes: rmse 2.569033e-06  max resid 5.352974e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.04483483 
    ## ... Procrustes: rmse 3.835202e-06  max resid 7.917893e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.04483483 
    ## ... Procrustes: rmse 1.808762e-06  max resid 3.902639e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.04483483 
    ## ... Procrustes: rmse 5.487316e-07  max resid 9.718967e-07 
    ## ... Similar to previous best
    ## Run 10 stress 0.04483484 
    ## ... Procrustes: rmse 8.571112e-06  max resid 1.730373e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.04483483 
    ## ... Procrustes: rmse 1.322656e-06  max resid 2.37637e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.174781 
    ## Run 13 stress 0.04483483 
    ## ... Procrustes: rmse 2.97407e-06  max resid 6.281312e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.2758772 
    ## Run 15 stress 0.04483483 
    ## ... Procrustes: rmse 3.619888e-06  max resid 7.496934e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.04483483 
    ## ... Procrustes: rmse 1.058356e-06  max resid 2.271466e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.04483483 
    ## ... Procrustes: rmse 1.542963e-06  max resid 3.254704e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.3658214 
    ## Run 19 stress 0.04483483 
    ## ... Procrustes: rmse 3.329236e-06  max resid 6.971648e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.04483483 
    ## ... Procrustes: rmse 1.181675e-06  max resid 2.588868e-06 
    ## ... Similar to previous best
    ## *** Solution reached

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

NMDS plots attempt to show ordinal distances between samples as
accurately as possible in two dimensions. It is important to report the
stress of these plots, because a high stress value means that the
algorithm had a hard time representing the distances between samples in
2 dimensions. The stress of this plot was good - it was .04 (generally
anything below .2 is considered acceptable).

Q4: Which datasets were used to make the two NMDS plots above? Based on
how these two plots look, which dataset are we going to move forward
with and why?

A4: The two datasets that were used were the subsampled data and the not
subsampled data. Based on the two plots, we can move forward with the
subsampled data because the patterns of diversity are very similar
between both plots.

Q5: Describe in a couple of sentences what patterns you see in the beta
diversity of the communities in the control and Ash Leachate treatments.

A5: While the ash leachate and control treatments began the experiment
with similar community composition, meaning that the beta diversity was
low. At day 4, the control and ash leachate treatment are quite
dissimilar compared to day 0, and the beta diversity accross the
treatments increases. In Santa Barbara, it looks like there is a lower
beta diversity accross treatments than San Diego on day 4 since the
points are closer together.

# Alpha Diversity

Estimating alpha diversity of microbial communities is
[problematic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC93182/) no
matter what you do.

We are going to calculate the Chao1 index for richness and the Shannon
diversity index.

**it is important to note that the alpha diversity values are not
interpretable as “real” numbers of anything (due to the nature of
amplicon data), but they can still be useful as relative metrics of
comparison. If Chao1 richness goes up, but Shannon diversity goes down,
it indicates that the sample may have more ASVs but is dominated by a
few of them.**

We will use the subsampled library, which retains estimates of the
species abundance of the real population while standardizing sampling
effort.

[subsampling and alpha diversity
paper](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)

[Chao1: nonparametric estimation of minimum community
richness](https://www.jstor.org/stable/4615964?seq=1#metadata_info_tab_contents)

``` r
richness <- estimate_richness(ps_min, measures = c("Chao1", "Shannon")) %>% 
  rownames_to_column(., var = "DNA_ID") %>% 
  mutate_at(vars(DNA_ID), str_replace_all, pattern = "171.", "171-") %>% 
   mutate_at(vars(DNA_ID), str_replace_all, pattern = "172.", "172-")
```

Let’s add the sample metadata into this dataframe

``` r
alphadiv <- left_join(richness, sample.tab %>% rownames_to_column(., var = "DNA_ID")) 
```

    ## Joining, by = "DNA_ID"

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

Boxes represent the 1.5 interquartile range, with the internal solid
line representing the median. Circles represent data points. p-values
are reported the non-parametric two sample Wilcoxon test, which tests
whether the means between two groups are equal (ns: p \> 0.05, \* : p≤
0.05, \*\* : p ≤ 0.01).

Difference in the alpha diversity indexes among conditions were tested
using pairwise Wilcoxon tests; p \< 0.05 was considered the threshold
significance for a difference between conditions.

From this plot we can see within the treatments that the richness (via
Chao index) of our samples significantly changed, while overall
diversity (via Shannon index) did not change. This suggest that while
richness decreased in both the control and ash leachate treatments, the
eveness was similar between the initial and final conditions.

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

From this plot we can see between the treatments that the richness of
the control samples were higher at the initial condition than the ash
leachate, suggesting that there may have been some quality control
issues as we would expect the initial samples to all have the same
richness. By timepoint 6, it looks like the richness was about the same
between the control and the ash leachate. Overall diversity was similar
between the treatments at the initial condition, but not by the end of
the experiment. The ash leachate samples at timepoint 6 may have been
less even.

Q6: Summarize the major takeaways from the alpha diversity plots you
generated.

A6: The major takeaway of the alpha diversity plots is that at time
point 6, the ash leachate treatment has a significantly lower alpha
diversity than the control and the ash leachate treatment also has a
lower richness than the control. The ash leachate treatment resulted in
lower richness and lower evenness.

# Who??

Which taxa were important? Which taxa were contributing to the change in
community compositon?

**Note: Recovered 16S rRNA gene copy numbers do not equal organism
abundance.**

That said, we can generate a heat map of our samples showing us how the
relative abundance of different taxonomic groups change…potentially
giving us a visual of which taxa are most important to the alpha and
beta diversity patterns we observed. First, we’re going to generate a
custom table that will be easier to work with than a phyloseq object.

## Generate relative abundances

Our data currently shows number gene copies recovered, so we’ll convert
to percentages (relative abundances)

``` r
ps_std <- transform_sample_counts(ps_min, function(x) x/sum(x))
#extract the relative abundance table and coerce into dataframe
ps_std.tab <- as(otu_table(ps_std), "matrix")
ps_std.df = as.data.frame(ps_std.tab) 
```

## Make table

``` r
#first coerce the taxa table into a data frame
tax.df <-  as.data.frame(tax.tab) 
#then combine the data frames
custom.tab <- tax.df %>% 
  rownames_to_column(., var = "asv") %>% 
  left_join(., ps_std.df %>% rownames_to_column(., var = "asv")) %>% 
  #create a new index of that combines the  class, order, family, and genus values, you can play around here!!
  mutate(#pcofg = paste(Phylum, "_", Class, "_", Order,"_", Family, "_", Genus),
         # pcof = paste(Phylum, "_", Class, "_", Order,"_", Family,),
         pco = paste(Phylum, "_", Class, "_", Order)) %>% 
  select(-c(asv:Genus)) %>% 
  # select(pcof,everything()) %>% 
  # group_by(pcof) %>% 
  select(pco,everything()) %>% 
  group_by(pco) %>% 
  #here we are combining the relative abundances based on our grouping
  summarise_at(vars(contains(c("ASH171", "ASH172"))), sum, na.rm = T) %>% 
  ungroup()
```

    ## Joining, by = "asv"

``` r
#save the row names and then make them into the column names
colnames <- custom.tab[,1] 

#transpose the dataframe so we can merge with the sample info table
t_custom.tab <-  as.data.frame(t(custom.tab[,-1]))
# colnames(t_custom.tab) <- colnames$pcof
colnames(t_custom.tab) <- colnames$pco

#merge
sweet.tab <- t_custom.tab %>% 
  rownames_to_column(., var = "sample") %>% 
  left_join(., sample.tab %>% rownames_to_column(., var = "sample") %>% select(sample, Experiment, Location, Bottle, Treatment, Timepoint, days, cells)) %>% 
  select(sample, Experiment:cells, everything())
```

    ## Joining, by = "sample"

``` r
relabund <- sweet.tab %>% 
  select(-c(sample:cells)) %>% 
  #remove groups that are completely absent
  .[ , colSums(.) > 0] %>% 
  #arrange by biggest contributors
  .[, order(colSums(-.))] %>% 
  bind_cols(sweet.tab %>% select(sample:cells), .)
```

## Heatmap

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

Q7: Who uniquely increased in relative abundance in the ash leachate
treatment that did not in the control? what about decrease?

A7: Oceanospirillales, Campylobacterales, and Alteromonadales all
increase in relative abundance in the ash leachate treatment, but not
really for the control. There were no noticable decreases in the ash
leachate that did not occur in the control as well, but the SAR11 clade
had a larger decrease in abundance in the ash leachate treatment than
the control. Synechococcales also had a larger decrease in relative
abundance in the ash leachate than the control.

Everything shown here is just a snapshot of what you can look at with
your community composition data. There are many other resources you can
use to get ideas of how to look at different aspects of your data,
including the [phyloseq tutorial](https://joey711.github.io/phyloseq/)
and [happy belly bioinformatics](https://astrobiomike.github.io). It’s
up to you and your questions!!

# Save and knit

``` r
saveRDS(sweet.tab, "~/Desktop/github/144l_students_2021/Output_Data/week7/Custom_ASV_Table.rds")
saveRDS(sub_ps, "~/Desktop/github/144l_students_2021/Output_Data/week7/phyloseq_obj.rds")
saveRDS(ps_min, "~/Desktop/github/144l_students_2021/Output_Data/week7/subsampled_phyloseq_obj.rds")
saveRDS(alphadiv, "~/Desktop/github/144l_students_2021/Output_Data/week7/alphadiv.rds")
```

# Stacked Barplots

<img src="211010_Phyloseq_submission_files/figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />
