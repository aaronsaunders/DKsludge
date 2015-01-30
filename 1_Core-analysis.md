# Wastewater treatment plants contain a conserved core community of bacteria


```r
source("R/functions.R")
```

```
## Loading required package: phyloseq
## Loading required package: ade4
## Loading required package: picante
## Loading required package: ape
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
## 
## Attaching package: 'vegan'
## 
## The following object is masked from 'package:ade4':
## 
##     cca
## 
## Loading required package: nlme
## Note: the specification for S3 class "AsIs" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "connection" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "file" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "pipe" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Note: the specification for S3 class "textConnection" in package 'BiocGenerics' seems equivalent to one from package 'RJSONIO': not turning on duplicate class definitions for this class.
## Loading required package: gridExtra
## Loading required package: grid
## 
## Attaching package: 'dplyr'
## 
## The following object is masked _by_ '.GlobalEnv':
## 
##     n
## 
## The following object is masked from 'package:nlme':
## 
##     collapse
## 
## The following objects are masked from 'package:plyr':
## 
##     arrange, desc, failwith, id, mutate, summarise
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
source("R/core_community.R")
sessionInfo()
```

```
## R version 3.0.3 (2014-03-06)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      graphics  grDevices utils     datasets  stats     methods  
## [8] base     
## 
## other attached packages:
##  [1] dplyr_0.1.1     gridExtra_0.9.1 phyloseq_1.7.3  picante_1.6-1  
##  [5] nlme_3.1-115    vegan_2.0-10    lattice_0.20-27 permute_0.8-3  
##  [9] ape_3.0-11      ade4_1.6-2      knitr_1.5       plyr_1.8       
## [13] reshape2_1.2.2  ggplot2_0.9.3.1
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.40.0         AnnotationDbi_1.24.0   
##  [3] assertthat_0.1          Biobase_2.22.0         
##  [5] BiocGenerics_0.8.0      biom_0.3.11            
##  [7] Biostrings_2.30.1       cluster_1.15.1         
##  [9] codetools_0.2-8         colorspace_1.2-4       
## [11] DBI_0.2-7               DESeq2_1.2.10          
## [13] dichromat_2.0-0         digest_0.6.3           
## [15] evaluate_0.5.1          foreach_1.4.1          
## [17] formatR_0.10            genefilter_1.44.0      
## [19] GenomicRanges_1.14.4    gtable_0.1.2           
## [21] igraph_0.7.0            IRanges_1.20.6         
## [23] iterators_1.0.6         labeling_0.2           
## [25] locfit_1.5-9.1          MASS_7.3-30            
## [27] Matrix_1.1-2-2          multtest_2.18.0        
## [29] munsell_0.4.2           parallel_3.0.3         
## [31] proto_0.3-10            RColorBrewer_1.0-5     
## [33] Rcpp_0.11.0             RcppArmadillo_0.4.000.2
## [35] RJSONIO_1.0-3           RSQLite_0.11.4         
## [37] scales_0.2.3            splines_3.0.3          
## [39] stats4_3.0.3            stringr_0.6.2          
## [41] survival_2.37-7         tools_3.0.3            
## [43] XML_3.98-1.1            xtable_1.7-1           
## [45] XVector_0.2.0
```


# Data

The source dataset contains all:

  1. Q3 samples from 2008 and 2008 from all plants. 
  2. time series from AAW
  

```r
identities <- as.character(c(94, 97, 99))
fname_template <- "data/otuXX/seqs_XX_otutable.biom"
filenames <- sapply(identities, function(id) gsub("XX", id, fname_template))

datasets <- lapply(identities, function(identity) LoadData(biompath = filenames[identity], 
    mapfpath = "data/mapfile.txt"))
names(datasets) <- paste0("otu", identities)

print(datasets_summary <- data.frame(samples = sapply(datasets, function(identity) sum(nsamples(identity))), 
    reads = sapply(datasets, function(identity) sum(sample_sums(identity))), 
    OTUs = sapply(datasets, function(identity) sum(ntaxa(identity)))))
```

```
##       samples   reads OTUs
## otu94      48 2374197 2541
## otu97      48 2374197 4586
## otu99      48 2374197 8276
```


### Datasets

 1) `coreDatasets`  Two samples from each plant from the summer 2008 and 2009 were used to calculate the core microbial community in the cross-section of Danish plants. 
 2) `tsDatasets`  All the samples from Aalborg West from 2006 and 2010 were used to calculate the core community in the time-series.


```r
subsample_depth <- 40000
# List containing phyloseq objects for core dataset at 94, 97 and 99%
# identity
coreDatasets <- lapply(datasets, function(identity) selectCoreDataset(identity, 
    depth = subsample_depth, seed = 1234))
```

```
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 187 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 389 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 863 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r
names(coreDatasets) <- paste0("otu", identities)
lapply(coreDatasets, function(identity) printDatasetStats(identity))
```

```
## [1] "number of plants: 13"
##      
##       2008 2009
##   AAE    1    1
##   AAW    1    1
##   BJE    1    1
##   EGA    1    1
##   EJB    1    1
##   HJO    1    1
##   HLV    1    1
##   ONE    1    1
##   ONW    1    1
##   RAN    1    1
##   SKI    1    1
##   SOE    1    1
##   VIB    1    1
## [1] 1040000
## [1] "total reads: 1040000"
## [1] "number of plants: 13"
##      
##       2008 2009
##   AAE    1    1
##   AAW    1    1
##   BJE    1    1
##   EGA    1    1
##   EJB    1    1
##   HJO    1    1
##   HLV    1    1
##   ONE    1    1
##   ONW    1    1
##   RAN    1    1
##   SKI    1    1
##   SOE    1    1
##   VIB    1    1
## [1] 1040000
## [1] "total reads: 1040000"
## [1] "number of plants: 13"
##      
##       2008 2009
##   AAE    1    1
##   AAW    1    1
##   BJE    1    1
##   EGA    1    1
##   EJB    1    1
##   HJO    1    1
##   HLV    1    1
##   ONE    1    1
##   ONW    1    1
##   RAN    1    1
##   SKI    1    1
##   SOE    1    1
##   VIB    1    1
## [1] 1040000
## [1] "total reads: 1040000"
```

```
## $otu94
## [1] "total reads: 1040000"
## 
## $otu97
## [1] "total reads: 1040000"
## 
## $otu99
## [1] "total reads: 1040000"
```

```r

# List containing phyloseq objects for time-series dataset at 94, 97 and 99%
# id
tseriesDatasets <- lapply(datasets, function(identity) selectAAWDataset(identity, 
    depth = subsample_depth))
```

```
## [1] "sample seed = 1234"
## [1] "subsampled at 40000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 1467 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## [1] "sample seed = 1234"
## [1] "subsampled at 40000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 2965 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## [1] "sample seed = 1234"
## [1] "subsampled at 40000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 5865 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r
names(tseriesDatasets) <- paste0("otu", identities)
lapply(tseriesDatasets, function(identity) printDatasetStats(identity))
```

```
## [1] "number of plants: 1"
##      
##       2006 2007 2008 2009 2010 2011
##   AAW    2    1    2    2    3    3
## [1] 520000
## [1] "total reads: 520000"
## [1] "number of plants: 1"
##      
##       2006 2007 2008 2009 2010 2011
##   AAW    2    1    2    2    3    3
## [1] 520000
## [1] "total reads: 520000"
## [1] "number of plants: 1"
##      
##       2006 2007 2008 2009 2010 2011
##   AAW    2    1    2    2    3    3
## [1] 520000
## [1] "total reads: 520000"
```

```
## $otu94
## [1] "total reads: 520000"
## 
## $otu97
## [1] "total reads: 520000"
## 
## $otu99
## [1] "total reads: 520000"
```


## 1. Core community

## Quantifying the core community


```r
coredataframe <- calcSummaryData(coreDatasets, identities, core_cutoff = 1)

coredata_by_id <- group_by(coredataframe, identities)
summarise(coredata_by_id, coretotalOTUs = sum(taxsum), coreOTUs = sum(taxsum[corestatus == 
    "core"]), percentcorereads = round(readprop[corestatus == "core"] * 100, 
    1))
```

```
## Source: local data frame [3 x 4]
## 
##   identities coretotalOTUs coreOTUs percentcorereads
## 1         99          7413       77             36.6
## 2         97          4197      100             54.3
## 3         94          2354       86             68.1
```

```r

coreplots <- plotCore(coredataframe)
otuplot <- coreplots[[1]]
readsplot <- coreplots[[2]]
print(Figure1 <- grid.arrange(otuplot, readsplot, nrow = 2))
```

![plot of chunk Figure 1: Core OTU conservation plot](figure/Figure_1:_Core_OTU_conservation_plot.png) 

```
## NULL
```

```r

pdf(file = "figs/Figure_1_core_otu_conservation40k.pdf")
Figure1
```

```
## NULL
```

```r
dev.off()
```

```
## pdf 
##   2
```


### Fraction of Transient OTUs (observed no more than 5 times)


```r
filter(coredataframe, numsamples < 6) %.% group_by(identities) %.% summarise(total_OTUs = sum(taxsum), 
    percent_reads = round(sum(readprop) * 100, 1))
```

```
## Source: local data frame [3 x 3]
## 
##   identities total_OTUs percent_reads
## 1         99       5796           8.9
## 2         97       3003           5.8
## 3         94       1546           2.3
```


### Fraction of Frequently occuring OTUs (observed 20 + times)


```r
filter(coredataframe, numsamples >= 20) %.% group_by(identities) %.% summarise(total_OTUs = sum(taxsum), 
    percent_reads = round(sum(readprop) * 100, 1))
```

```
## Source: local data frame [3 x 3]
## 
##   identities total_OTUs percent_reads
## 1         99        371          69.1
## 2         97        321          79.3
## 3         94        242          86.2
```


## Cumulative abundance distribution

The distribution of abundances in bacterial communities is uneven - some 
organisms are abundant some are rare.

Plotting percent cumulative abundance vs. rank ordered OTUs describes this 
distribution.


```r
cum_abun_plot <- plot_CumulRankAbundance(dataset = coreDatasets[["otu94"]], 
    fname = "figs/figure2_cum_abun_plot_94.pdf", dominant_fraction = 0.8)
```

```
## [1] "percent for count of 1: 0.0025"
## [1] "80% cutoff: 65"
## [1] "quantiles for #otus for 80% of reads"
##     0%    25%    50%    75%   100% 
##  19.00  48.50  59.50  72.25 134.00 
## [1] "abundance of otu at 80th percentile: 0.15 (+- 0.08, n=26)"
## [1] "wrote plot to: figs/figure2_cum_abun_plot_94.pdf"
```

```r
print(cum_abun_plot)
```

![plot of chunk Figure 2: Cumulative abundance curve](figure/Figure_2:_Cumulative_abundance_curve.png) 


### Distribution of abundances per OTU

There is a lot of variation in the abundances but many of the abundant organisms are consitently abundant. These are the organisms that are core. The size of the core will increase as a function of depth for these consistent organisms.



```r
plotTopN(coreDatasets[["otu94"]], topN = 50, plot_filename = "figs/Figure_3_top50_id94_boxplot.pdf", 
    core_cutoff = 1)
```

```
## sampled  top 50 OTUs
## 62.9% of the total reads
```

![plot of chunk Figure 3: boxplot the top OTUs](figure/Figure_3:_boxplot_the_top_OTUs.png) 


### Temporal stability of a single plant through time

Time series in AAW

What proportion of the core OTUs are core in: 
 - the two AAW samples used for the core
 - all the AAW samples


```r
tsdataframe <- calcSummaryData(tseriesDatasets, identities, core_cutoff = 1)
group_by(tsdataframe, identities) %.% summarise(coretotalOTUs = sum(taxsum), 
    coreOTUs = sum(taxsum[corestatus == "core"]), percentcorereads = round(readprop[corestatus == 
        "core"] * 100, 1))
```

```
## Source: local data frame [3 x 4]
## 
##   identities coretotalOTUs coreOTUs percentcorereads
## 1         99          2411      329             87.3
## 2         97          1621      254             90.6
## 3         94          1074      190             93.4
```

```r

tscoreplot <- plotCore(tsdataframe)
tsotuplot <- tscoreplot[[1]]
tsreadsplot <- tscoreplot[[2]]
print(tscoreplot <- grid.arrange(tsotuplot, tsreadsplot, nrow = 2))
```

![plot of chunk calc AAW timeseries core](figure/calc_AAW_timeseries_core.png) 

```
## NULL
```

```r

pdf(file = "figs/Figure_S2_aawtimeseries_core_conservation_40k.png")
tscoreplot
```

```
## NULL
```

```r
dev.off()
```

```
## pdf 
##   2
```



### Binning of OTUs by observation frequency and frequency of high-abundance


```r
dataset <- coreDatasets[["otu94"]]
core_prop <- 1

dominant <- sapply(sample_names(dataset), function(samplename) fill_ha_data(dataset, 
    samplename))
row.names(dominant) <- taxa_names(dataset)

observed <- as.data.frame(otu_table(transform_sample_counts(dataset, function(x) ifelse(x > 
    0, 1, 0))))

summary.df <- data.frame(OTU = taxa_names(dataset), nHA = apply(dominant, 1, 
    sum), nObs = apply(observed, 1, sum), median = apply(otu_table(dataset), 
    1, median), geomean = apply(otu_table(dataset), 1, function(x) round(exp(mean(log(x))), 
    1)), max = apply(otu_table(dataset), 1, max), min = apply(otu_table(dataset), 
    1, min), n1per = apply(otu_table(dataset), 1, function(x) sum(x > 400)))
# how does median relate to nHA and nObs?
ggplot(summary.df, aes(nObs, median)) + geom_point(alpha = 0.2) + scale_y_log10()
```

![plot of chunk Dominant OTUs](figure/Dominant_OTUs1.png) 

```r
ggplot(summary.df, aes(nHA, median)) + geom_point() + scale_y_log10()
```

![plot of chunk Dominant OTUs](figure/Dominant_OTUs2.png) 


These plots are the justification for having nHA > 10 as the cutoff for significance. 

## Figure 4

Sets up the empty Figure 4 plot. The data was filled manually using Inkscape.


```r
summary.df <- mutate(summary.df, group = factor(cut(nHA, breaks = c(-1, 0, 9, 
    25, 26), labels = c("4", "3", "2", "1")), levels = 1:4), Obsclass = cut(nObs, 
    breaks = c(0, 19, 25, 26), labels = c("ob1", "ob20", "ob26")))
summary.df <- cbind(summary.df, as.data.frame(otu_table(dataset)))
summary.df <- arrange(summary.df, desc(nHA), desc(nObs), desc(median))
# Empty plot
print(Figure4 <- plotFigure4())
```

![plot of chunk Figure 4: Dominant vs. Frequency](figure/Figure_4:_Dominant_vs__Frequency.png) 

```r
ggsave(filename = "figs/Figure4.pdf", plot = Figure4, width = 8, height = 8, 
    units = "cm")
```


## Data for plotting onto Figure 4

How many OTUs/reads are in each category


```r
r <- melt(select(summary.df, OTU, group, Obsclass, AMPA057:AMPA724), id.vars = c("OTU", 
    "group", "Obsclass"), variable.name = "sample", value.name = "count") %.% 
    group_by(Obsclass, group) %.% summarise(readpercent = round((sum(count)/(26 * 
    40000)) * 100, 1), nOTUs = n_distinct(OTU))
acast(r, group ~ Obsclass, value.var = "nOTUs", fun.aggregate = sum, margins = TRUE)
```

```
##        ob1 ob20 ob26 (all)
## 1        0    0    3     3
## 2        0   10   51    61
## 3      130   94   28   252
## 4     1982   52    4  2038
## (all) 2112  156   86  2354
```

```r
acast(r, group ~ Obsclass, value.var = "readpercent", fun.aggregate = sum, margins = TRUE)
```

```
##        ob1 ob20 ob26 (all)
## 1      0.0  0.0 25.7  25.7
## 2      0.0  6.1 36.9  43.0
## 3      7.6 10.4  5.4  23.4
## 4      6.2  1.7  0.2   8.1
## (all) 13.8 18.2 68.2 100.2
```




```r
# NB rounding error# 3. transiently highly-abundant
group3otus <- as.vector(with(summary.df, summary.df[group == 3, "OTU"]))
tHA_percent <- round(as.data.frame(otu_table(dataset)[group3otus, ]/40000 * 
    100), 1)


df <- data.frame(rank = 1:26, percentHA = sort(colSums(tHA_percent), decreasing = TRUE))
plant <- samData(dataset)[row.names(df), "plant"]
df$plant <- as.character(plant$plant)  # crazy phyloseq object!!

# TODO fix the label on this? why does the sample_data object persist after
# as.character() ggplot(df, aes(x=rank, y= percentHA)) + geom_point(stat =
# 'identity') + scale_x_discrete(labels= plant) + xlab('samples') +
# ylab('transiently HA read percent')

# num trans.abun OTUs per sample
apply(tHA_percent, 2, function(sample) sum(sample > 1))
```

```
## AMPA057 AMPA578 AMPA562 AMPA110 AMPA615 AMPA047 AMPA032 AMPA034 AMPA632 
##       0       1       4       0       2       0       1       1       2 
## AMPA566 AMPA568 AMPA581 AMPA563 AMPA564 AMPA726 AMPA580 AMPA570 AMPA056 
##       5       6       5       8       4       8       0       2       0 
## AMPA053 AMPA055 AMPA561 AMPA725 AMPA054 AMPA102 AMPA727 AMPA724 
##       3       1       0       2       1       0       2       1
```

```r
samData(dataset)[names(sort(colSums(tHA_percent), decreasing = TRUE)), "plant"]
```

```
## Sample Data:        [26 samples by 1 sample variables]:
##         plant
## AMPA563   HLV
## AMPA726   VIB
## AMPA725   ONE
## AMPA724   ONW
## AMPA566   ONW
## AMPA564   VIB
## AMPA568   SOE
## AMPA615   SOE
## AMPA562   SKI
## AMPA581   ONE
## AMPA580   HLV
## AMPA727   EJB
## AMPA632   BJE
## AMPA570   RAN
## AMPA561   EGA
## AMPA056   SKI
## AMPA047   AAE
## AMPA578   HJO
## AMPA034   AAW
## AMPA110   AAE
## AMPA032   AAW
## AMPA054   EGA
## AMPA053   BJE
## AMPA055   HJO
## AMPA102   EJB
## AMPA057   RAN
```

```r

HA.tran <- subset(summary.df, nObs <= 20 & nHA > 0 & (n1per > 0))
tax_table(dataset)[as.vector(HA.tran$OTU), 1:6]
```

```
## Taxonomy Table:     [22 taxa by 6 taxonomic ranks]:
##      Kingdom    Phylum             Class                
## 1401 "Bacteria" "Actinobacteria"   "Actinobacteria"     
## 1628 "Bacteria" "Actinobacteria"   "Actinobacteria"     
## 2389 "Bacteria" "Nitrospirae"      "Nitrospira"         
## 792  "Bacteria" "Proteobacteria"   "Deltaproteobacteria"
## 511  "Bacteria" NA                 NA                   
## 1720 "Bacteria" "Proteobacteria"   "Betaproteobacteria" 
## 2095 "Bacteria" "Chloroflexi"      "Chloroflexi"        
## 2504 "Bacteria" "Cyanobacteria"    "4C0d-2"             
## 1988 "Bacteria" "Gemmatimonadetes" "Gemmatimonadetes"   
## 785  "Bacteria" "Actinobacteria"   "Actinobacteria"     
## 637  "Bacteria" NA                 NA                   
## 1178 "Bacteria" "Chloroflexi"      "Anaerolineae"       
## 2169 "Bacteria" "Chloroflexi"      NA                   
## 1658 "Bacteria" "Proteobacteria"   "Alphaproteobacteria"
## 863  "Bacteria" "Proteobacteria"   "Alphaproteobacteria"
## 1892 "Bacteria" "Proteobacteria"   "Gammaproteobacteria"
## 365  "Bacteria" "Chloroflexi"      "Anaerolineae"       
## 2207 "Bacteria" "Actinobacteria"   "Actinobacteria"     
## 703  "Bacteria" "Planctomycetes"   "Planctomycea"       
## 1419 "Bacteria" "Proteobacteria"   "Gammaproteobacteria"
## 1127 "Bacteria" "Proteobacteria"   "Deltaproteobacteria"
## 1748 "Bacteria" "Chloroflexi"      "Chloroflexi"        
##      Order              Family               Genus              
## 1401 "Acidimicrobiales" "Microthrixaceae"    "Microthrix"       
## 1628 "Actinomycetales"  "Intrasporangiaceae" "Tetrasphaera_etal"
## 2389 "Nitrospirales"    "Nitrospiraceae"     "Nitrospira"       
## 792  "Myxococcales"     "OM27"               NA                 
## 511  NA                 NA                   NA                 
## 1720 "Methylophilales"  "Methylophilaceae"   "Methylotenera"    
## 2095 "Roseiflexales"    "Kouleothrixaceae"   "Kouleothrix"      
## 2504 "mle1-12"          NA                   NA                 
## 1988 "Gemmatimonadales" NA                   NA                 
## 785  "Acidimicrobiales" NA                   NA                 
## 637  NA                 NA                   NA                 
## 1178 "Anaerolineales"   "Anaerolinaceae"     "Anaerolinea"      
## 2169 NA                 NA                   NA                 
## 1658 "Rhodospirillales" "Rhodospirillaceae"  "Defluviicoccus"   
## 863  "Rhodobacterales"  "Rhodobacteraceae"   "Amaricoccus"      
## 1892 "Thiotrichales"    "Thiotrichaceae"     "Thiothrix"        
## 365  "WCHB1-50"         NA                   NA                 
## 2207 "Actinomycetales"  "Intrasporangiaceae" "Tetrasphaera_etal"
## 703  "Gemmatales"       "Isosphaeraceae"     "Nostocoida"       
## 1419 "Pseudomonadales"  "Moraxellaceae"      "Psychrobacter"    
## 1127 "Myxococcales"     "Haliangiaceae"      NA                 
## 1748 "Roseiflexales"    NA                   NA
```

```r
nrow(tax_table(dataset)[as.vector(HA.tran$OTU), 1:6])
```

```
## [1] 22
```



## Single sequence resolution on the abundant OTUs



# export the OTU data to a table for the Table S2


```r
# list of tax names
eco.core.taxa <- as.vector(with(summary.df, summary.df[HAclass %in% c("HA10", 
    "HA26"), "OTU"]))
```

```
## Error: object 'HAclass' not found
```

```r
trans.abun.noncore.taxa <- as.vector(with(summary.df, summary.df[HAclass == 
    "HA1" & n1per > 0, "OTU"]))
```

```
## Error: object 'HAclass' not found
```

```r

length(trans.abun.noncore.taxa)
```

```
## Error: object 'trans.abun.noncore.taxa' not found
```

```r

# lists of taxnames ecocore and abun.noncore
taxa.to.export <- c(ecocore.taxa, trans.abun.noncore.taxa)
```

```
## Error: object 'ecocore.taxa' not found
```

```r
outstats <- select(summary.df, 1:12) %.% filter(OTU %in% taxa.to.export)
```

```
## Error: binding not found: 'taxa.to.export'
```

```r
outtaxtable <- as.data.frame(tax_table(coreDatasets[["data.94"]])[taxa.to.export, 
    1:6])
```

```
## Error: tax_table slot is empty.
```

```r
outdata <- cbind(outstats, outtaxtable) %.% arrange(desc(median))
```

```
## Error: object 'outstats' not found
```

```r

table(outdata$group)
```

```
## Error: object 'outdata' not found
```

```r

write.table(outdata, file = "figs/Table_S2_newcore_otutable.txt", sep = "\t", 
    quote = FALSE, row.names = FALSE, col.names = TRUE)
```

```
## Error: object 'outdata' not found
```



### Nitrotoga and Nitrospira


```r
# Compare replicate data

data.97 <- datasets[["otu97"]]
amplibs <- sample_data(data.97)$SampleID[sample_data(datasets[["otu97"]])$sample_id == 
    "293"]
data.97 <- prune_samples(x = data.97, samples = as.character(amplibs))
table(sam_data(data.97)$sample_id, sam_data(data.97)$dna_id)
```

```
##      
##       148 149 150 151 152
##   293   1   6   1   1   1
```

```r
data.97 <- rarefy_even_depth(data.97, rngseed = 1234, sample.size = 14000, trimOTUs = TRUE)
```

```
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 2  samples removed because they contained fewer reads than `sample.size`.
## Up to first five removed samples are: 
## AMPA232	AMPA161
## ...
## 3742 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r

data.97 <- transform_sample_counts(physeq = data.97, function(x) x/sum(x) * 
    100)
NOB <- subset_taxa(data.97, (Family == "Gallionellaceae") | (Genus == "Nitrospira"))

plot_bar(NOB, "Genus", facet_grid = . ~ SampleID) + geom_bar(aes(fill = Genus), 
    stat = "identity", position = "stack")
```

![plot of chunk NOB](figure/NOB1.png) 

```r

plot_bar(NOB, "Genus", facet_grid = dna_id ~ SampleID) + geom_bar(aes(fill = Genus), 
    stat = "identity", position = "stack")
```

![plot of chunk NOB](figure/NOB2.png) 

```r

# compare core samples
data.97 <- coreDatasets[["otu97"]]
data.97 <- transform_sample_counts(physeq = data.97, function(x) x/sum(x) * 
    100)
NOB <- subset_taxa(data.97, (Genus == "Nitrotoga_etal") | (Genus == "Nitrospira"))

p <- plot_bar(NOB, "Genus", facet_grid = year ~ plant) + geom_bar(aes(fill = Genus), 
    stat = "identity", position = "stack")
levels(p$data$Genus) <- c("Nitrospira", "Nitrotoga")
p$data$plant_name <- factor(gsub(x = p$data$plant_name, pattern = "oe", replacement = "ø"))
p$data$plant_name <- factor(gsub(x = p$data$plant_name, pattern = "aa", replacement = "å"))

totals <- as.data.frame(select(p$data, OTU, Sample, Abundance, year, Genus, 
    plant_name) %.% filter(year == 2009, Genus == "Nitrotoga") %.% acast(plant_name ~ 
    Genus, value.var = "Abundance", sum))
totals$plant_name <- row.names(totals)
plant_names_by <- with(totals, totals[order(Nitrotoga, decreasing = TRUE), "plant_name"])
p$data$plant_name <- factor(p$data$plant_name, levels = plant_names_by)

# formatted in greyscale for publication
p2 <- ggplot(p$data, aes(x = plant_name, y = Abundance, fill = Genus)) + geom_bar(position = position_dodge(), 
    stat = "identity") + facet_grid(year ~ .) + theme_bw() + ylab(label = "Read abundance (%)") + 
    xlab(label = "Plant") + scale_fill_manual(values = c("grey50", "black")) + 
    annotate("text", x = "Aalborg East", y = 2.5, label = "2008", size = 2) + 
    theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8, 
        vjust = 0.2), axis.text.x = element_text(size = 6, hjust = 1, vjust = 1, 
        angle = 45), axis.text.y = element_text(size = 6), axis.ticks.x = element_line(size = 0.3), 
        axis.ticks.y = element_line(size = 0.3), panel.grid.minor = element_blank(), 
        strip.background = element_rect(linetype = "blank", fill = "white"), 
        strip.text.y = element_blank(), axis.line = element_line(size = 0.3), 
        legend.title = element_blank(), legend.text = element_text(size = 5, 
            face = "italic"), legend.key.size = unit(0.4, "lines"), legend.key = element_rect(size = 1, 
            colour = "white"), legend.justification = c(1, 1), legend.position = c(1.055, 
            1.086))

fname <- "figs/Figure_7_NOB_core.pdf"
ggsave(plot = p2, file = fname, width = 8, height = 6.37, units = "cm")

# Compare time series data
ts.97 <- tseriesDatasets[["otu97"]]
ts.97 <- transform_sample_counts(physeq = ts.97, function(x) x/sum(x) * 100)
NOB <- subset_taxa(ts.97, (Genus == "Nitrotoga_etal") | (Genus == "Nitrospira"))

p <- plot_bar(NOB, "Genus", facet_grid = month ~ year) + geom_bar(aes(fill = Genus), 
    stat = "identity", position = "stack")
levels(p$data$Genus) <- c("Nitrospira", "Nitrotoga")
p
```

![plot of chunk NOB](figure/NOB3.png) 

```r
p$data$month <- 
p2 <- ggplot(p$data, aes(x = Genus, y = Abundance, fill = Genus)) + geom_bar(position = position_dodge(), 
    stat = "identity") + facet_grid(month ~ year) + theme_bw() + ylab(label = "Abundance (%)") + 
    xlab(label = "Plant") + theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8, 
    vjust = 0.2), axis.text.x = element_text(size = 6, hjust = 1, vjust = 1, 
    angle = 45), axis.text.y = element_text(size = 6), axis.ticks.x = element_line(size = 0.3), 
    axis.ticks.y = element_line(size = 0.3), strip.text = element_text(size = 6), 
    panel.grid.minor = element_blank(), strip.background = element_rect(linetype = "blank", 
        fill = "white"), axis.line = element_line(size = 0.3), legend.position = "none")
```

```
## Error: replacement has 9 rows, data has 52
```

```r

fname <- "figs/Figure_S5.pdf"
ggsave(plot = p2, file = fname, width = 16, units = "cm")
```

```
## Saving 16 x 17.8 cm image
```

```r

d <- read.delim("data/MIDAS_combined.csv", row.names = 1)
```

```
## Warning: cannot open file 'data/MIDAS_combined.csv': No such file or
## directory
```

```
## Error: cannot open the connection
```

```r
mean(d$T_lc, na.rm = TRUE)
```

```
## Error: object 'd' not found
```

```r

# TODO change data to Table S1 yearly temp means in DK plants
d.means <- ddply(d, .(Plant, Quarter), summarize, mean = mean(T_lc, na.rm = TRUE), 
    sd = sd(T_lc, na.rm = TRUE), n = sum(!is.na(T_lc)))
```

```
## Error: object 'd' not found
```

```r
d.means <- d.means[d.means$n > 3, ]
```

```
## Error: object 'd.means' not found
```

```r
d.means <- subset(d.means, Quarter %in% c(1, 3))
```

```
## Error: object 'd.means' not found
```

```r

ggplot(d.means, aes(Plant, mean, color = factor(Quarter))) + geom_point() + 
    theme(axis.text.x = theme_text(angle = 45, hjust = 1))
```

```
## Error: object 'd.means' not found
```

```r

summarySE(d.means, measurevar = "mean", groupvars = c("Quarter"), na.rm = TRUE)
```

```
## Error: object 'd.means' not found
```




```r
# what percentage of reads have a genus level classification?
data.97 <- coreDatasets[["otu97"]]
data.97.genus <- tax_glom(data.97, taxrank = "Genus", NArm = FALSE)
ntaxa(data.97.genus)
```

```
## [1] 635
```

```r
total_reads <- sum(sample_sums(data.97.genus))
data.97.genus_na <- subset_taxa(data.97.genus, is.na(Genus))
ntaxa(data.97.genus_na)
```

```
## [1] 299
```

```r
total_genusna_reads <- sum(sample_sums(data.97.genus_na))
total_genusna_reads/total_reads
```

```
## [1] 0.6219
```

