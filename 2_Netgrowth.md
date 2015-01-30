# 2. Immigration and netgrowth


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
biomfname = "data/otu97/netgrowth_seqs_97_otutable.biom"
netgrowthdata <- LoadData(biompath = biomfname, mapfpath = "data/mapfile.txt")
```

```
## Warning: No greengenes prefixes were found. 
## Consider using parse_taxonomy_default() instead if true for all OTUs. 
## Dummy ranks may be included among taxonomic ranks now.
## Warning: No greengenes prefixes were found. 
## Consider using parse_taxonomy_default() instead if true for all OTUs. 
## Dummy ranks may be included among taxonomic ranks now.
```

```r

data.k <- SampleEvenDepth(mydata = netgrowthdata, depth = 20000, seed = 1234)
```

```
## [1] "sample seed = 1234"
## [1] "subsampled at 20000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 450 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r
print(sample_sums(data.k))
```

```
## AMPA636 AMPA645 AMPA649 AMPA642 AMPA648 AMPA609 AMPA635 AMPA644 AMPA610 
##   20000   20000   20000   20000   20000   20000   20000   20000   20000 
## AMPA640 AMPA639 AMPA647 
##   20000   20000   20000
```



Investigate/visualise the influent data. Compare to sludge. Compare to McLellan.

Calculate the net growth of each OTU

Compare net growth to abundance and find out which of the OTUs' growth is 
overestimated by their abundance.

Can we find single sequences that are abundant in the influent and sludge.


```r
ord.cca <- ordinate(physeq = data.k, method = "CCA", distance = "bray")
plot1 <- plot_ordination(physeq = data.k, ordination = ord.cca, shape = "plant", 
    type = "samples", color = "sample_type") + geom_text(aes(label = pair), 
    color = "black", size = 5, vjust = 2) + geom_point(size = 5)

ord.rda <- ordinate(physeq = data.k, method = "RDA", distance = "euclidian")
plot2 <- plot_ordination(physeq = data.k, ordination = ord.rda, type = "samples", 
    shape = "plant") + geom_text(aes(label = pair, color = sample_type, vjust = 2)) + 
    geom_point(size = 5)

ggsave(path = "figs", filename = "CCA_AS_vs_WW.png", plot = plot1, dpi = 400)
```

```
## Saving 7 x 7 in image
```


The euclidian PCA splits the samples by sample type in the first PC and roughly by plant in the second.

The CCA has the same trend but splits the plants Aalborg and Hjørring on the second axis more clearly.


```r
ids <- as.character(c(94, 97, 99))
fname_template <- "data/otuXX/netgrowth_seqs_XX_otutable.biom"
fnames <- sapply(ids, function(id) gsub("XX", id, fname_template))

datasets <- lapply(ids, function(id) LoadData(biompath = fnames[id], mapfpath = "data/mapfile.txt"))
```

```
## Warning: No greengenes prefixes were found. 
## Consider using parse_taxonomy_default() instead if true for all OTUs. 
## Dummy ranks may be included among taxonomic ranks now.
## Warning: No greengenes prefixes were found. 
## Consider using parse_taxonomy_default() instead if true for all OTUs. 
## Dummy ranks may be included among taxonomic ranks now.
## Warning: No greengenes prefixes were found. 
## Consider using parse_taxonomy_default() instead if true for all OTUs. 
## Dummy ranks may be included among taxonomic ranks now.
```

```r
ids <- sapply(ids, function(id) paste0("otu", id))
names(datasets) <- ids

sapply(datasets, function(dataset) sum(sample_sums(dataset)))
```

```
##  otu94  otu97  otu99 
## 811161 811161 811161
```

```r
sapply(datasets, function(dataset) sum(nsamples(dataset)))
```

```
## otu94 otu97 otu99 
##    12    12    12
```

```r
sapply(datasets, function(dataset) sum(ntaxa(dataset)))
```

```
## otu94 otu97 otu99 
##  2214  3609  5584
```

```r

depth <- 40000

rarefyDataset <- function(mydata, depth = FALSE, seed = FALSE) {
    mydata <- subset_samples(mydata, type == "WW")
    mydata <- rarefy_even_depth(mydata, rngseed = seed, sample.size = depth, 
        trimOTUs = TRUE)
    mydata
}

WWDatasets <- lapply(ids, function(id) rarefyDataset(datasets[[id]], depth, 
    1234))
```

```
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 1058 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 1728 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 2560 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r
names(WWDatasets) <- ids

lapply(WWDatasets, function(dataset) printDatasetStats(dataset))
```

```
## [1] "number of plants: 3"
##      
##       2012 2013
##   AAE    0    2
##   AAW    1    1
##   HJO    0    2
## [1] 240000
## [1] "total reads: 240000"
## [1] "number of plants: 3"
##      
##       2012 2013
##   AAE    0    2
##   AAW    1    1
##   HJO    0    2
## [1] 240000
## [1] "total reads: 240000"
## [1] "number of plants: 3"
##      
##       2012 2013
##   AAE    0    2
##   AAW    1    1
##   HJO    0    2
## [1] 240000
## [1] "total reads: 240000"
```

```
## $otu94
## [1] "total reads: 240000"
## 
## $otu97
## [1] "total reads: 240000"
## 
## $otu99
## [1] "total reads: 240000"
```

```r

source("R/core_community.R")

coredataframe <- calcData(WWDatasets, ids, core_cutoff = 1)
```

```
## Error: could not find function "calcData"
```

```r

OTUtotals <- ddply(coredataframe, "identities", function(id) c(totalOTUs = sum(id$taxsum), 
    coreOTUs = id$taxsum[id$corestatus == "core"], percentcoreOTUs = round(id$taxsum[id$corestatus == 
        "core"]/sum(id$taxsum) * 100, 1), percentreads = round(id$readprop[id$corestatus == 
        "core"] * 100, 1)))
```

```
## Error: object 'coredataframe' not found
```

```r
print(OTUtotals)
```

```
## Error: object 'OTUtotals' not found
```

```r

p <- plotCore(coredataframe, fig_fname = "figs/infeff_conservation40k.png", 
    cumulative = FALSE)
```

```
## Loading required package: RColorBrewer
```

```
## Error: object 'coredataframe' not found
```

```r
otuplot <- p[[1]]
```

```
## Error: object 'p' not found
```

```r
readsplot <- p[[2]]
```

```
## Error: object 'p' not found
```

```r

pdf(file = "figs/ww_core40k.pdf")
grid.arrange(otuplot, readsplot, nrow = 2)
```

```
## Error: object 'otuplot' not found
```

```r
dev.off()
```

```
## pdf 
##   2
```

```r
grid.arrange(otuplot, readsplot, nrow = 2)
```

```
## Error: object 'otuplot' not found
```


## Compare the top OTUs in sludge and influent


```r
source("R//net_growth_rate.R")
```

```
## Warning: cannot open file 'scr/functions.R': No such file or directory
```

```
## Error: cannot open the connection
```

```r
plot3 <- compareTop20AsWw(data.k)
```

```
## Error: could not find function "compareTop20AsWw"
```

```r
ggsave(path = "figs", filename = "top20OTUsWW_vs_AS.png", plot = plot3, width = 8, 
    dpi = 600)
```

```
## Error: object 'plot3' not found
```

```r

# what percentage of reads are from the top genera?
ds <- WWDatasets[["otu94"]]
ds.p <- transform_sample_counts(ds, function(x) x/sum(x))
n <- 20

topNww_otus <- names(sort(taxa_sums(ds.p), decreasing = TRUE)[1:n])
topN.ww <- subset_taxa(ds, taxa_names(ds) %in% topNww_otus)
df.rawdata <- as.data.frame(otu_table(topN.ww))
otus <- row.names(df.rawdata)
df.outdata <- data.frame(OTU = otus, as.data.frame(tax_table(topN.ww)))
row.names(df.outdata) <- 1:n
df.outdata$median <- apply(df.rawdata, 1, median)
df.outdata$geomean <- apply(df.rawdata, 1, function(x) round(exp(mean(log(x))), 
    1))
df.outdata$max <- apply(df.rawdata, 1, max)
df.outdata$min <- apply(df.rawdata, 1, min)
df.outdata <- df.outdata[order(df.outdata$median, decreasing = TRUE), ]
df.outdata
```

```
##     OTU  Kingdom         Phylum                 Class             Order
## 6   650 Bacteria Proteobacteria Epsilonproteobacteria Campylobacterales
## 18 2009 Bacteria Proteobacteria    Betaproteobacteria   Burkholderiales
## 3    56 Bacteria     Firmicutes               Bacilli   Lactobacillales
## 5   643 Bacteria Proteobacteria   Gammaproteobacteria   Pseudomonadales
## 19 2046 Bacteria Proteobacteria Epsilonproteobacteria Campylobacterales
## 15 1514 Bacteria Proteobacteria   Gammaproteobacteria     Aeromonadales
## 17 1913 Bacteria  Bacteroidetes         Flavobacteria  Flavobacteriales
## 11  954 Bacteria  Bacteroidetes           Bacteroidia     Bacteroidales
## 16 1691 Bacteria Proteobacteria   Gammaproteobacteria   Pseudomonadales
## 2    40 Bacteria Proteobacteria   Gammaproteobacteria   Pseudomonadales
## 9   796 Bacteria  Bacteroidetes           Bacteroidia     Bacteroidales
## 13 1203 Bacteria     Firmicutes            Clostridia     Clostridiales
## 12  981 Bacteria     Firmicutes            Clostridia     Clostridiales
## 4    62 Bacteria  Bacteroidetes           Bacteroidia     Bacteroidales
## 1    20 Bacteria Proteobacteria    Betaproteobacteria     Rhodocyclales
## 7   677 Bacteria     Firmicutes            Clostridia     Clostridiales
## 8   700 Bacteria   Fusobacteria          Fusobacteria   Fusobacteriales
## 14 1509 Bacteria     Firmicutes               Bacilli   Lactobacillales
## 10  830 Bacteria Proteobacteria   Gammaproteobacteria   Pseudomonadales
## 20 2139 Bacteria Proteobacteria   Gammaproteobacteria   Pseudomonadales
##                Family              Genus median geomean   max  min
## 6  Campylobacteraceae         Arcobacter 5947.5  5683.1  7036 4162
## 18     Comamonadaceae Limnohabitans_etal 3896.0  4241.8 10930 2270
## 3   Carnobacteriaceae               <NA> 3187.0  3368.6  5438 2594
## 5       Moraxellaceae      Acinetobacter 3091.0  2589.6  4565 1193
## 19 Campylobacteraceae         Arcobacter 1326.5  1204.7  1551  894
## 15     Aeromonadaceae          Aeromonas  972.0   931.0  1427  583
## 17  Flavobacteriaceae     Flavobacterium  898.0   834.2  1482  407
## 11     Bacteroidaceae        Bacteroides  671.5   601.8   875  261
## 16   Pseudomonadaceae        Pseudomonas  621.0   559.7  1904  154
## 2    Pseudomonadaceae               <NA>  589.5   582.9  5542  147
## 9      Prevotellaceae         Prevotella  544.5   466.8   703  173
## 13    Lachnospiraceae            Blautia  533.5   525.6   893  268
## 12    Ruminococcaceae   Faecalibacterium  515.0   484.3   709  219
## 4      Bacteroidaceae        Bacteroides  443.5   384.0   527  149
## 1      Rhodocyclaceae      Dechloromonas  435.0   347.3  1030  109
## 7     Lachnospiraceae          Roseburia  378.0   363.2   483  194
## 8    Fusobacteriaceae               <NA>  366.0   344.6   947  130
## 14   Streptococcaceae        Lactococcus  264.5   712.1  8160  149
## 10      Moraxellaceae      Psychrobacter  189.0   268.9  1193  148
## 20      Moraxellaceae      Enhydrobacter  109.0   108.3  1145    8
```






```r
source("R/net_growth_rate.R")
```

```
## Warning: cannot open file 'scr/functions.R': No such file or directory
```

```
## Error: cannot open the connection
```

```r
df.netgrowth <- calcNetGrowthData(data.k, k_lower_bound = -0.2)
```

```
## Error: could not find function "calcNetGrowthData"
```

```r

plist <- plotNetGrowthDistribution(df.netgrowth)
```

```
## Error: could not find function "plotNetGrowthDistribution"
```

```r

print(plist[[1]])
```

```
## Error: object 'plist' not found
```

```r
print(plist[[2]])
```

```
## Error: object 'plist' not found
```

```r
grid.arrange(plist[[1]], plist[[2]], nrow = 2)
```

```
## Error: object 'plist' not found
```

```r

dev.off()
```

```
## null device 
##           1
```

```r
ggsave(path = "figs", filename = "cumabun_by_netgrowth.pdf", plot = plist[[2]], 
    width = 8, height = 6, units = "cm")
```

```
## Error: object 'plist' not found
```

```r


ggplot(data = df.netgrowth, aes(x = km, fill = class)) + geom_histogram(binwidth = 0.05) + 
    # coord_cartesian(ylim=c(0, 500)) +
scale_fill_brewer(palette = "Set1") + labs(x = NULL, y = "OTU count") + theme(legend.position = "none")
```

```
## Error: object 'df.netgrowth' not found
```

```r

df.netgrowth$km_class <- cut(df.netgrowth$km, breaks = c(-0.22, -0.1, 0, 0.03, 
    10), labels = c("inactive", "low", "active", "max"))
```

```
## Error: object 'df.netgrowth' not found
```

```r
df.netgrowth[(df.netgrowth$prop_ww == 0), 11] = "max"
```

```
## Error: object 'df.netgrowth' not found
```

```r

as_by_class <- with(df.netgrowth, aggregate(prop_as * 100, by = list(pair, km_class), 
    FUN = length))
```

```
## Error: object 'df.netgrowth' not found
```

```r
names(as_by_class) = c("pair", "class", "n_otus")
```

```
## Error: object 'as_by_class' not found
```

```r

df.netgrowth.abun <- subset(df.netgrowth, df.netgrowth$prop_as > 0.001)
```

```
## Error: object 'df.netgrowth' not found
```

```r
res <- with(df.netgrowth.abun, aggregate(prop_as * 100, by = list(pair, km_class), 
    FUN = length))
```

```
## Error: object 'df.netgrowth.abun' not found
```

```r
as_by_class$n_abun_otus <- res$x
```

```
## Error: object 'res' not found
```

```r

res <- with(df.netgrowth, aggregate(prop_as * 100, by = list(pair, km_class), 
    FUN = sum))
```

```
## Error: object 'df.netgrowth' not found
```

```r
as_by_class$read_percent <- res$x
```

```
## Error: object 'res' not found
```

```r

res <- with(df.netgrowth.abun, aggregate(prop_as * 100, by = list(pair, km_class), 
    FUN = sum))
```

```
## Error: object 'df.netgrowth.abun' not found
```

```r
as_by_class$abun_read_percent <- res$x
```

```
## Error: object 'res' not found
```

```r

as_by_class$class_per_abun <- round(as_by_class$abun_read_percent/as_by_class$read_percent * 
    100, 0)
```

```
## Error: object 'as_by_class' not found
```

```r

ggplot(data = as_reads_by_class, aes(x = pair, y = read_percent, fill = class)) + 
    geom_bar(stat = "identity")
```

```
## Error: object 'as_reads_by_class' not found
```

```r

n_per_kmclass <- with(as_reads_by_class, aggregate(read_percent, by = list(class), 
    FUN = sum))
```

```
## Error: object 'as_reads_by_class' not found
```

```r

mean_per_kmclass <- with(as_reads_by_class, aggregate(read_percent, by = list(class), 
    FUN = mean))
```

```
## Error: object 'as_reads_by_class' not found
```

```r

sd_per_kmclass <- with(as_reads_by_class, aggregate(read_percent, by = list(class), 
    FUN = sd))
```

```
## Error: object 'as_reads_by_class' not found
```

```r

# How many reads are active (km > 0)
active <- subset(as_by_class, class %in% c("active", "max"))
```

```
## Error: object 'as_by_class' not found
```

```r
act_totals <- with(active, aggregate(read_percent, by = list(pair), FUN = sum))
```

```
## Error: object 'active' not found
```

```r
round(mean(act_totals$x), 0)
```

```
## Error: object 'act_totals' not found
```

```r
sd(act_totals$x)
```

```
## Error: object 'act_totals' not found
```

```r

act_totals <- with(active, aggregate(n_abun_otus, by = list(pair), FUN = sum))
```

```
## Error: object 'active' not found
```

```r
round(mean(act_totals$x), 0)
```

```
## Error: object 'act_totals' not found
```

```r
sd(act_totals$x)
```

```
## Error: object 'act_totals' not found
```

```r

# How many reads are inactive (km < 0)
inactive <- subset(as_by_class, class %in% c("low", "inactive"))
```

```
## Error: object 'as_by_class' not found
```

```r
inact_totals <- with(inactive, aggregate(read_percent, by = list(pair), FUN = sum))
```

```
## Error: object 'inactive' not found
```

```r
round(mean(inact_totals$x), 0)
```

```
## Error: object 'inact_totals' not found
```

```r
sd(inact_totals$x)
```

```
## Error: object 'inact_totals' not found
```

```r

inact_totals <- with(inactive, aggregate(n_otus, by = list(pair), FUN = sum))
```

```
## Error: object 'inactive' not found
```

```r
round(mean(inact_totals$x), 0)
```

```
## Error: object 'inact_totals' not found
```

```r
sd(inact_totals$x)
```

```
## Error: object 'inact_totals' not found
```

```r

sum(inactive$n_abun_otus)
```

```
## Error: object 'inactive' not found
```

```r
inact_totals <- with(inactive, aggregate(n_abun_otus, by = list(pair), FUN = sum))
```

```
## Error: object 'inactive' not found
```

```r
round(mean(inact_totals$x), 0)
```

```
## Error: object 'inact_totals' not found
```

```r
sd(inact_totals$x)
```

```
## Error: object 'inact_totals' not found
```

```r

df.netgrowth.abuninact <- subset(df.netgrowth.abun, df.netgrowth.abun$km < 0)
```

```
## Error: object 'df.netgrowth.abun' not found
```

```r
length(unique(df.netgrowth.abuninact$OTU))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r
length(unique(df.netgrowth.abuninact$tax_name))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r
barchart(table(df.netgrowth.abuninact$tax_name))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r
hist(table(df.netgrowth.abuninact$OTU))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r

abun_by_taxon <- with(df.netgrowth.abuninact, aggregate(prop_as, by = list(tax_name), 
    FUN = sum))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r
abun_by_taxon$x <- round(abun_by_taxon$x * 100/6, 1)
```

```
## Error: object 'abun_by_taxon' not found
```

```r
abun_by_taxon <- abun_by_taxon[order(abun_by_taxon$x), ]
```

```
## Error: object 'abun_by_taxon' not found
```

```r
names(abun_by_taxon) <- c("name", "mean")
```

```
## Error: object 'abun_by_taxon' not found
```

```r
per_plant_tax_counts <- with(df.netgrowth.abuninact, aggregate(prop_as * 100, 
    by = list(pair, tax_name), FUN = max))
```

```
## Error: object 'df.netgrowth.abuninact' not found
```

```r
names(per_plant_tax_counts) <- c("pair", "taxon", "max_count")
```

```
## Error: object 'per_plant_tax_counts' not found
```

```r
per_plant_tax_counts <- per_plant_tax_counts[order(per_plant_tax_counts$max_count), 
    ]
```

```
## Error: object 'per_plant_tax_counts' not found
```

```r


#### 

df.netgrowth$tax_name <- tax_table(data.k)[df.netgrowth$OTU, "Genus"]
```

```
## Error: error in evaluating the argument 'i' in selecting a method for function '[': Error: object 'df.netgrowth' not found
```

```r
df.netgrowth$taxlab <- ifelse((df.netgrowth$prop_ww > 0.01) & (df.netgrowth$km > 
    -0.2) & (df.netgrowth$km < 0.01), as.character(df.netgrowth$tax_name), "")
```

```
## Error: object 'df.netgrowth' not found
```

```r
# percent reads near k max
cumsum(df.netgrowth$prop_as/6)[sum(df.netgrowth$k > 0.02857)]
```

```
## Error: object 'df.netgrowth' not found
```

```r
high_growth_rate <- df.netgrowth[df.netgrowth$k > 0.02857, ]
```

```
## Error: object 'df.netgrowth' not found
```

```r
max(high_growth_rate$prop_ww)
```

```
## Error: object 'high_growth_rate' not found
```

```r

# percent of reads with k > -0.15
cumsum(df.netgrowth$prop_as/6)[sum(df.netgrowth$km < -0)]
```

```
## Error: object 'df.netgrowth' not found
```


The assumption is that the net growth rate at steady state will be equal to the 
max net growth rate (1/35 = 0.0286).

There are a lot of OTUs that are at or near the upper limit. 


```r
printNetGrowthStats(df.netgrowth)
```

```
## Error: could not find function "printNetGrowthStats"
```


so the are on average 64% of the reads in the sludge detected in the influent, 
but this number is highly variable probably because the sampling of the sludge
taxa in the nfluent is near the level of detection.

Can I recalculate the data resampling the wastewater at a higher cutoff.


```r
data.ww <- subset_samples(netgrowthdata, type == "WW")
sample_sums(data.ww)
```

```
## AMPA649 AMPA635 AMPA644 AMPA610 AMPA639 AMPA647 
##   74369   71700   62439   67449   75525   58297
```

```r
data.k.ww <- SampleEvenDepth(mydata = data.ww, depth = 58000, seed = 1234)
```

```
## [1] "sample seed = 1234"
## [1] "subsampled at 58000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 1651 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r

data.as <- subset_samples(netgrowthdata, type == "AS")
data.k.as <- SampleEvenDepth(mydata = data.as, depth = 40000, seed = 1234)
```

```
## [1] "sample seed = 1234"
## [1] "subsampled at 40000 reads"
## `set.seed(1234)` was used to initialize repeatable random subsampling.
## Please record this for your records so others can reproduce.
## Try `set.seed(1234); .Random.seed` for the full vector
## ...
## 1052 OTUs were removed because they are no longer 
##  present in any sample after random subsampling
## ...
```

```r

data.k.uneven <- merge_phyloseq(data.k.ww, data.k.as)
df.netgrowth.uneven <- calcNetGrowthData(data.k.uneven)
```

```
## Error: could not find function "calcNetGrowthData"
```

```r

printNetGrowthStats(df.netgrowth.uneven)
```

```
## Error: could not find function "printNetGrowthStats"
```


This did not make much difference...



```r
myBreaks <- c(1, 10)

ggplot(data = df.netgrowth, aes(x = prop_as * 100, y = km, color = prop_ww * 
    100, label = taxlab)) + geom_point() + geom_text(size = 3, color = "black", 
    vjust = 1, hjust = 1) + scale_x_log10(breaks = c(0.01, 0.1, 1, 10), limits = c(0.01, 
    15)) + scale_colour_gradient(trans = "log", low = "white", high = "red", 
    breaks = myBreaks, labels = myBreaks) + facet_wrap(pair ~ plant) + labs(y = "bounded net growth rate", 
    x = "biomass abundance (%)", colour = "ww abundance (%)")
```

```
## Error: object 'df.netgrowth' not found
```

```r
ggsave(filename = "figs/sludge_vs_k_all.png")
```

```
## Saving 7 x 7 in image
```

```r


myBreaks <- c(1, 10)
df.netgrowth$tax_name <- makeTaxLabels(df.netgrowth$OTU, data.k)
```

```
## Error: object 'df.netgrowth' not found
```

```r

df.netgrowth$taxlab <- NA
```

```
## Error: object 'df.netgrowth' not found
```

```r
df.netgrowth$taxlab <- with(df.netgrowth, ifelse((prop_ww >= 0.01) & (prop_as >= 
    0.0012), tax_name, taxlab))
```

```
## Error: object 'df.netgrowth' not found
```

```r

p <- plot_netgrowth(mydata = df.netgrowth, value = "AAE-2")
```

```
## Error: could not find function "plot_netgrowth"
```

```r
p
```

```
## Error: object 'p' not found
```

```r

info <- c(as.character(samData(data.k)$pair), "all")
per_plant_plots <- lapply(info, FUN = plot_netgrowth)
```

```
## Error: object 'plot_netgrowth' not found
```

```r
names(per_plant_plots) <- info
```

```
## Error: object 'per_plant_plots' not found
```

```r

```



```r
class_by_pair <- group_by(df.netgrowth, pair, class)
```

```
## Error: object 'df.netgrowth' not found
```

```r
d <- summarise(class_by_pair, total = n()) %.% arrange(class, pair)
```

```
## Error: object 'class_by_pair' not found
```

```r
d1 <- d[d$class == "Both", "total"]/d[d$class == "AS", "total"] * 100
```

```
## Error: object 'd' not found
```

```r
names(d1) <- unique(d$pair)
```

```
## Error: object 'd' not found
```

```r
# result per pair
d1
```

```
## Error: object 'd1' not found
```

```r

by_class <- group_by(df.netgrowth, class)
```

```
## Error: object 'df.netgrowth' not found
```

```r
d2 <- summarise(by_class, total = n()) %.% arrange(class)
```

```
## Error: object 'by_class' not found
```

```r
d2
```

```
## Error: object 'd2' not found
```

```r
# overall fraction
d2[d2$class == "Both", "total"]/d2[d2$class == "AS", "total"] * 100
```

```
## Error: object 'd2' not found
```






```r
source("R/net_growth_rate.R")
```

```
## Warning: cannot open file 'scr/functions.R': No such file or directory
```

```
## Error: cannot open the connection
```

```r

data.genus <- tax_glom(physeq = data.k, taxrank = "Genus", NArm = TRUE)
df.genus.k <- calcNetGrowthData(data.genus, k_lower_bound = -0.2)
```

```
## Error: could not find function "calcNetGrowthData"
```

```r

plotlist <- plotNetGrowthDistribution(df.genus.k)
```

```
## Error: could not find function "plotNetGrowthDistribution"
```

```r

plotlist[[1]]
```

```
## Error: object 'plotlist' not found
```

```r
plotlist[[2]]
```

```
## Error: object 'plotlist' not found
```

```r

dev.off()
```

```
## null device 
##           1
```

```r

myBreaks <- c(1, 10)
df.genus.k$tax_name <- makeTaxLabels(df.genus.k$OTU, data.k)
```

```
## Error: object 'df.genus.k' not found
```

```r

df.genus.k$taxlab <- NA
```

```
## Error: object 'df.genus.k' not found
```

```r
df.genus.k$taxlab <- with(df.genus.k, ifelse((prop_ww >= 0.01) & (prop_as >= 
    0.0012), tax_name, taxlab))
```

```
## Error: object 'df.genus.k' not found
```

```r
df.genus.k$taxlab <- with(df.genus.k, ifelse((prop_as >= 0.02), tax_name, taxlab))
```

```
## Error: object 'df.genus.k' not found
```

```r

p <- plot_netgrowth(mydata = df.genus.k, value = "AAE-2")
```

```
## Error: could not find function "plot_netgrowth"
```

```r
p
```

```
## Error: object 'p' not found
```

```r

```




# Which taxa have a strong negative growth rate?


```r

df.netgrowth.neg <- subset(df.netgrowth, (k <= -0) & (prop_ww > 0.001)  )
```

```
## Error: object 'df.netgrowth' not found
```

```r
neg_otus         <- unique(df.netgrowth.neg$OTU)
```

```
## Error: object 'df.netgrowth.neg' not found
```

```r

neg_taxnames     <- as.data.frame(tax_table(data.k)[neg_otus])[1:6]
```

```
## Error: error in evaluating the argument 'i' in selecting a method for function '[': Error: object 'neg_otus' not found
```

```r
neg_taxnames$OTU <- row.names(neg_taxnames)
```

```
## Error: object 'neg_taxnames' not found
```

```r
stats            <- ddply(df.netgrowth.neg, "OTU", function (otu) 
                                        c("max" = max(otu$k),
                                          "min" = min(otu$k),
                                         "range"= max(otu$k) -  min(otu$k),
                                    "total_as" = sum(otu$prop_as),
                                     "total_w" = sum(otu$prop_ww)))
```

```
## Error: object 'df.netgrowth.neg' not found
```

```r

neg_taxa <- merge(neg_taxnames, stats)
```

```
## Error: object 'neg_taxnames' not found
```

```r

neg_taxa.min <- neg_taxa[ order(neg_taxa$min), ]
```

```
## Error: object 'neg_taxa' not found
```

```r
write.table(neg_taxa.min, file="exp/neg_taxa.tsv", sep="\t", 
              quote=FALSE, row.names=FALSE, col.names=TRUE)
```

```
## Error: object 'neg_taxa.min' not found
```

```r


pl <- ggplot(data= df.netgrowth.neg, 
             aes(y=abs(k), 
             x =reorder(OTU, abs(k), order= TRUE, function(x) median(x))) ) +
    scale_y_reverse(name = "k") + 
    scale_x_discrete(name   = "highly negative OTUs", 
                  labels = function(x) 
                     df.netgrowth.neg[ df.netgrowth.neg$OTU == x, "taxname"]) +
    coord_cartesian(ylim = c(0, 0.5)) +
    geom_boxplot(  ) + # aes(fill = tax_name)
    theme(axis.text.x  = element_text(size = 8,  hjust = 1, vjust = 1, 
                                      angle = 45)) +
    theme(axis.title.x = element_text(size = 12)) + 
    theme(axis.title.y = element_text(size = 12, vjust = 0.2)) +  
    theme(axis.text.y  = element_text(size = 12)) +
    theme(legend.text  = element_text(size = 10)) + 
    theme(legend.title = element_text(size = 12)) +
    theme(legend.justification = c(1, 0))
```

```
## Error: object 'df.netgrowth.neg' not found
```

```r
#+    theme(legend.position      = "none")
  
pl
```

```
## Error: object 'pl' not found
```

```r



df.netgrowth.asneg <- subset(df.netgrowth, k <= -0.1)
```

```
## Error: object 'df.netgrowth' not found
```

```r
stats <- ddply(df.netgrowth.asneg, .(pair), summarize, 
               sumas=sum(prop_as),
               sumww=sum(prop_ww))
```

```
## Error: object 'df.netgrowth.asneg' not found
```

```r
mean(stats$sumas) * 100
```

```
## Error: object 'stats' not found
```

```r
sd  (stats$sumas) * 100
```

```
## Error: object 'stats' not found
```

```r

df.netgrowth.asneg.abun <- subset(df.netgrowth, (k <= -0.1) & (prop_as > 0.005))
```

```
## Error: object 'df.netgrowth' not found
```

```r

stats <- ddply(df.netgrowth.asneg.abun, .(tax_name, pair), summarize, 
               sumas=round(sum(prop_as * 100), 1),
               sumww=round(sum(prop_ww * 100), 1))
```

```
## Error: object 'df.netgrowth.asneg.abun' not found
```

```r
stats
```

```
## Error: object 'stats' not found
```


The two big players in the neg k taxa are Clostridia and Bacteroidales. These
are typically anaerobes and therefore these are likely inactive.

The Bacteroides had a k of -0.2 which is considerably less than the the -0.04
we have used as a lower bound.

THere are some OTUs of known sludge genera that actually have a low k.

Thauera, Defluvicoccus, Thiothrix

TODO: Move the lower bound to say -0.2 and look at the 


## What is the range of the OTUs? 

Are there some in the higher ranges with a larger range?

## not sure what I am doing here...



```r



calcKSummary <- function(ps, df) {
    result <- as.data.frame(tax_table(ps))[1:6]
    result$OTU <- row.names(result)
    stats <- ddply(df, "OTU", function(otu) c(max = max(otu$k), min = min(otu$k), 
        range = max(otu$k) - min(otu$k), total_as = sum(otu$prop_as), total_w = sum(otu$prop_ww)))
    result <- merge(result, stats)
}


plotKSummary <- function(df) {
    a <- ggplot(data = df, aes(x = k)) + geom_histogram(binwidth = 0.005) + 
        coord_cartesian(ylim = c(0, 50), xlim = c(-0.5, 0.03)) + scale_fill_brewer(palette = "Set1") + 
        labs(x = NULL, y = "OTU count") + theme(legend.position = "none")
    a
    
    df <- df[order(df$km, decreasing = TRUE), ]
    prop_at_zerok <- cumsum(df$prop_as/6)[sum(df$km > 0)]
    df <- df[order(df$km, decreasing = FALSE), ]
    
    b <- ggplot(data = df, aes(x = km, y = cumsum(prop_as/6) * 100)) + scale_x_continuous(limit = c(-0.041, 
        0.03), breaks = seq(-0.04, 0.03, 0.01)) + coord_cartesian(ylim = c(0, 
        10)) + geom_step() + geom_hline(yintercept = (1 - prop_at_zerok) * 100, 
        color = "red") + labs(x = "bounded Net growth rate (k)", y = "% cumulative abundance")
    b
    list(a, b)
}


Lachnospir <- subset_taxa(data.k, Family == "Lachnospiraceae")

Lachnospir.otus <- taxa_names(Lachnospir)
df.Lachnospir <- df.netgrowth[df.netgrowth$OTU %in% Lachnospir.otus, ]
```

```
## Error: object 'df.netgrowth' not found
```

```r

stats <- ddply(df.Lachnospir, .(tax_name, pair), summarize, sumas = round(sum(prop_as * 
    100), 1), sumww = round(sum(prop_ww * 100), 1), kmean = round(mean(k), 3))
```

```
## Error: object 'df.Lachnospir' not found
```

```r
stats
```

```
## Error: object 'stats' not found
```

```r

Lach.summary <- calcKSummary(data.k, df.Lachnospir)
```

```
## Error: object 'df.Lachnospir' not found
```

```r
Lach.plots <- plotKSummary(df.Lachnospir)
```

```
## Error: object 'df.Lachnospir' not found
```

```r

a <- subset(df.Lachnospir, (df$prop_as > 0.001) & (df$prop_ww > 0.001))
```

```
## Error: object 'df.Lachnospir' not found
```

```r

ggplot(data = a, aes(x = k)) + geom_histogram(binwidth = 0.005) + coord_cartesian(ylim = c(0, 
    50), xlim = c(-0.15, 0.03)) + scale_fill_brewer(palette = "Set1") + labs(x = NULL, 
    y = "OTU count") + theme(legend.position = "none")
```

```
## Error: object 'a' not found
```

```r

Bacteroides <- subset_taxa(data.k, Family == "Bacteroidaceae")
Bacteroides.otus <- taxa_names(Bacteroides)
df.Bacteroides <- df.netgrowth[df.netgrowth$OTU %in% Bacteroides.otus, ]
```

```
## Error: object 'df.netgrowth' not found
```

```r

stats <- ddply(df.Bacteroides, .(tax_name, pair), summarize, sumas = round(sum(prop_as * 
    100), 1), sumww = round(sum(prop_ww * 100), 1), kmean = round(mean(k), 3))
```

```
## Error: object 'df.Bacteroides' not found
```

```r
stats
```

```
## Error: object 'stats' not found
```

```r

Bacteroides.summary <- calcKSummary(data.k, df.Bacteroides)
```

```
## Error: object 'df.Bacteroides' not found
```

```r
Bacteroides.plots <- plotKSummary(df.Bacteroides)
```

```
## Error: object 'df.Bacteroides' not found
```

```r

a <- subset(df.Bacteroides, (df$prop_as > 0.001) & (df$prop_ww > 0.001))
```

```
## Error: object 'df.Bacteroides' not found
```

```r

ggplot(data = a, aes(x = k)) + geom_histogram(binwidth = 0.005) + coord_cartesian(ylim = c(0, 
    50), xlim = c(-0.15, 0.03)) + scale_fill_brewer(palette = "Set1") + labs(x = NULL, 
    y = "OTU count") + theme(legend.position = "none")
```

```
## Error: object 'a' not found
```

```r

df.Bacteroides <- calcNetGrowthData(Bacteroides)
```

```
## Error: could not find function "calcNetGrowthData"
```

```r
Bact.summary <- calcKSummary(data.k, df.Bacteroides)
```

```
## Error: object 'df.Bacteroides' not found
```

```r
Bact.plots <- plotKSummary(df.Bacteroides)
```

```
## Error: object 'df.Bacteroides' not found
```

```r
Bact.plots
```

```
## Error: object 'Bact.plots' not found
```

```r



Arcobacter <- subset_taxa(data.k, Genus == "Arcobacter")
Arco.otus <- taxa_names(Arcobacter)
df.Arcobacter <- df.netgrowth[df.netgrowth$OTU %in% Arco.otus, ]
```

```
## Error: object 'df.netgrowth' not found
```

```r

Arco.summary <- calcKSummary(data.k, df.Arcobacter)
```

```
## Error: object 'df.Arcobacter' not found
```

```r
Arco.plots <- plotKSummary(df.Arcobacter)
```

```
## Error: object 'df.Arcobacter' not found
```

```r
Arco.plots
```

```
## Error: object 'Arco.plots' not found
```

```r


Accumulibacter <- subset_taxa(data.k, Genus == "CandidatusAccumulibacter")
Accumulibacter.otus <- taxa_names(Accumulibacter)
df.Accumulibacter <- df.netgrowth[df.netgrowth$OTU %in% Accumulibacter.otus, 
    ]
```

```
## Error: object 'df.netgrowth' not found
```

```r

ggplot(df.Accumulibacter, aes(y = prop_as, x = prop_ww, color = OTU)) + geom_point()
```

```
## Error: object 'df.Accumulibacter' not found
```

```r

Accumulibacter.summary <- calcKSummary(data.k, df.Accumulibacter)
```

```
## Error: object 'df.Accumulibacter' not found
```

```r
Accumulibacter.plots <- plotKSummary(df.Accumulibacter)
```

```
## Error: object 'df.Accumulibacter' not found
```

```r
Accumulibacter.plots
```

```
## Error: object 'Accumulibacter.plots' not found
```

```r

data.k.p <- transformSampleCounts(data.k, function(x) x/sum(x) * 100)
Tetrasphaera <- subset_taxa(data.k.p, Genus == "Tetrasphaera_etal")
Tetrasphaera.otus <- taxa_names(Tetrasphaera)
df.Tetrasphaera <- df.netgrowth[df.netgrowth$OTU %in% Tetrasphaera.otus, ]
```

```
## Error: object 'df.netgrowth' not found
```

```r
p <- plot_bar(Tetrasphaera, "Genus", fill = "OTU", facet_grid = type ~ pair)
```

```
## Error: error in evaluating the argument 'j' in selecting a method for function '[': Error in apply(!apply(TT, 2, is.na), 2, any) : 
##   dim(X) must have a positive length
## Calls: which -> apply
```

```r
p$data$OTU <- factor(p$data$OTU)
```

```
## Error: object 'p' not found
```

```r
p
```

```
## Error: object 'p' not found
```

```r

Tetrasphaera.summary <- calcKSummary(data.k, df.Tetrasphaera)
```

```
## Error: object 'df.Tetrasphaera' not found
```

```r
Tetrasphaera.plots <- plotKSummary(df.Tetrasphaera)
```

```
## Error: object 'df.Tetrasphaera' not found
```

```r
Tetrasphaera.plots
```

```
## Error: object 'Tetrasphaera.plots' not found
```

```r

AOB <- subset_taxa(data.k, Family == "Nitrosomonadaceae")
AOB <- transform_sample_counts(physeq = AOB, function(x) x/sum(x))

p <- plot_bar(AOB, "Genus", fill = "OTU", facet_grid = type ~ pair)
p$data$OTU <- factor(p$data$OTU)
p
```

```
## Warning: Removed 7 rows containing missing values (position_stack).
## Warning: no non-missing arguments to min; returning Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: position_stack requires constant width: output may be incorrect
## Warning: Removed 7 rows containing missing values (position_stack).
## Warning: no non-missing arguments to min; returning Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: position_stack requires constant width: output may be incorrect
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


## NOB in influent and effluent


```r
data.p <- transform_sample_counts(physeq = data.k, function(x) x/sum(x) * 100)
NOB <- subset_taxa(data.p, (Genus == "Nitrotoga_etal") | (Genus == "Nitrospira"))
NOB.prop <- transform_sample_counts(physeq = NOB, function(x) x/sum(x) * 100)

p <- plot_bar(NOB, "Genus", fill = "OTU", facet_grid = type ~ pair)
p$data$OTU <- factor(p$data$OTU)
p$data$Genus <- gsub(p$data$Genus, pattern = "_etal", replacement = "")
p <- p + ylab("Absolute abundance") + xlab(label = "") + theme(axis.text.x = element_blank(), 
    legend.position = "none")

p2 <- plot_bar(NOB.prop, "Genus", fill = "OTU", facet_grid = type ~ pair)
p2$data$OTU <- factor(p2$data$OTU)
p2$data$Genus <- gsub(p2$data$Genus, pattern = "_etal", replacement = "")
p2 <- p2 + ylab("Relative abundance") + xlab(label = "") + theme(legend.position = "none")

pdf(file = "figs/NOB_influent_sludge.pdf")
grid.arrange(p, p2, nrow = 2)
```

```
## Warning: Removed 7 rows containing missing values (position_stack).
## Warning: no non-missing arguments to min; returning Inf
## Warning: no non-missing arguments to max; returning -Inf
## Warning: position_stack requires constant width: output may be incorrect
```

```r
dev.off()
```

```
## pdf 
##   2
```

```r

d <- select(p$data, OTU, pair, type, Abundance) %.% filter(Abundance > 0.1)
mutate(d, k = df.netgrowth$k[OTU])
```

```
## Error: binding not found: 'df.netgrowth'
```

```r
d_group <- group_by(d, OTU)
dcast(d.group, OTU ~ type + pair, value.var = "Abundance", fun.aggregate = sum)
```

```
## Error: object 'd.group' not found
```

