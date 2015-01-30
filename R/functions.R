require(phyloseq)
packageVersion("phyloseq")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(dplyr)


# Data functions ----------------------------------------------------------

LoadData <- function(biompath, mapfpath){
  # Helper functions that wraps the steps of compiling a phyloseq object.
  #
  # Args:
  #   biompath: Path to a BIOM file
  #   mapfpath: Path to the qiime mapping OR sample_data file
  #
  # Returns:
  #   a phyloseq object 
  mapdata  <- import_qiime_sample_data(mapfpath)
  biomdata <- import_biom(BIOMfilename=biompath, parseFunction = parse_taxonomy_greengenes )
  mydata   <- merge_phyloseq(biomdata, sample_data(mapdata))
  mydata
}  

selectCoreDataset <- function(mydata, depth= FALSE, seed= FALSE){
  # Selects the samples that are compared for the core calculation. 
  #  - 2 samples from each plant
  #  - 2008, 2009 in quarter 3
  #
  # Args:
  #   mydata: A phyloseq object
  #   depth: the number of reads to resample
  #
  # Returns:
  #   a phyloseq object with samples at even depth 
  
  mydata.two <- subset_samples(mydata, 
                               (year %in% c(2008, 2009)) &
                                 (quarter == 3 ) &
                                 (!(plant %in% c("FOR", "FAA", "AVE", "MAR"))) &
                                 (!(ext_protocol == "Spin Kit")) )
  if (depth){
    mydata.two <-  rarefy_even_depth(mydata.two, rngseed = seed,
                                     sample.size = depth, trimOTUs= TRUE)  
  }
  mydata.two
}

selectAAWDataset <- function(myData, depth = 24000){
  # Selects the samples that are compared for the core calculation. 
  #  - 2 samples from each plant
  #  - 2008, 2009 in quarter 3
  #
  # Args:
  #   mydata: A phyloseq object
  #   depth: the number of reads to resample
  #
  # Returns:
  #   a phyloseq object with samples at even depth 
  mydata.subset <- subset_samples(physeq=myData, 
                                  (plant == "AAW") & 
                                    (ext_protocol == "phenol") &
                                    (year %in% 2006:2011) )
  nondups <- !duplicated(sample_data(mydata.subset)$date)
  mydata.subset <- prune_samples(x=mydata.subset, samples=nondups)
  mydata.subset <- SampleEvenDepth(mydata.subset, depth=depth)  
}

printDatasetStats <- function(data){
  print(paste("number of plants:", length(unique(sample_data(data)$plant))))
  print(table(sample_data(data)$plant, sample_data(data)$year))
  total_reads <- sum(sample_sums(data))
  print(total_reads)
  print(paste("total reads:", total_reads))
}


makeTaxLabels <- function(OTU, mydata){
  # Vector wrapper for "makeLabel" 
  # Makes a string label using the lowest informative tax level
  #
  # Args:
  #   OTUs: vector of OTU numbers
  #   mydata: the phyloseq object with the tax table
  #
  # Returns:
  #   a vector of a tax names
  sapply(OTU, FUN=makeTaxLabel, mydata)
}

makeTaxLabel <- function(OTU, mydata){
  # Makes a string label using the lowest informative tax level
  #
  # Args:
  #   OTU: OTU number
  #   mydata: the phyloseq object with the tax table
  #
  # Returns:
  #   a tax name
  OTU <- as.character(OTU)  # the OTU numbers are stored as character
  taxstrings <- as.character(tax_table(mydata)[OTU])
  empty_strings <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax_name <- NA
  tax_level <- length(taxstrings)  # start at lowest tax level
  
  while(is.na(tax_name) | 
          (tax_name %in% empty_strings)){
    tax_name  <- taxstrings[tax_level]
    tax_level <- tax_level -1
  }
  tax_name <- as.character(tax_name)
  tax_name
}

findBestTaxLevel <- function(OTU, mydata){
  # Makes a string label using the lowest informative tax level
  #
  # Args:
  #   OTU: OTU number
  #   mydata: the phyloseq object with the tax table
  #
  # Returns:
  #   a tax name
  OTU <- as.character(OTU)  # the OTU numbers are stored as character
  taxstrings <- as.character(tax_table(mydata)[OTU])
  empty_strings <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax_name <- NA
  tax_level <- length(taxstrings)  # start at lowest tax level
  
  while(is.na(tax_name) | 
          (tax_name %in% empty_strings)){
    tax_name  <- taxstrings[tax_level]
    tax_level <- tax_level -1
  }
  tax_level + 1 
}

SampleEvenDepth <- function(mydata, depth, seed=1234) {
  # A helper function around phyloseq rarefy_even_depth to remove samples that 
  # are less then the cutoff and to remove taxa that are zero after resampling
  #
  # Args:
  #   mydata: A phyloseq object containing samples to be rarefied.
  #   depth: the number of reads to resample
  #   seed: the random number seed for reproducible resampling
  #
  # Returns:
  #   a phyloseq object with samples at even depth  
  mydata <- prune_samples(sample_sums(mydata) > depth, mydata)
  print(paste(c("sample seed = ", seed), collapse=""))
  print(paste(c("subsampled at ", depth, " reads"), collapse=""))
  mydata.even <- rarefy_even_depth(mydata, sample.size = depth, 
                                   rngseed = seed, trimOTUs= TRUE)
  return(mydata.even)
}

resample <- function(depth, seqlib, nrare, seed=1) {
  # Multiple rarefaction of a sample 
  #
  # Args:
  #   depth: the number of reads to resample
  #   seqlib: A phyloseq object of a single library
  #   nrare: the number of resamplings
  #   seed: random number seed for reproducible resampling
  #
  # Returns:
  #   a phyloseq object with samples at even depth 
  res.matrix <- matrix(rep(0, times=ntax*nrare), nrow=ntax)
  colnames(res.matrix)  <- paste("rep", 1:nrare, sep="_")
  row.names(res.matrix) <- row.names(otu_table(seqlib))
  
  for (n in 1:nrare){
    set.seed(seed+(n*2))
    r = RarefyEvenDepth(seqlib, sample.size = depth)
    res.matrix[, n] = otu_table(r) / depth
  }
  res.matrix 
}

writeOTUList <- function(mydata, filename){
  # write a file of OTU names one per line
  #
  # mydata = phyloseq object
  # filename
  # 
  n <- taxa_names(mydata)
  write.table(n, file=filename, sep="\t", 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Analysis Functions ------------------------------------------------------



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

CalcRankOrderMatrix <- function(mydata, depth){
  # Calculates a matrix from the otu_table with each sample sorted by abundance.
  # Then summed to cumulative abundance.
  #
  # Args:
  #   mydata: a phyloseq object
  #   depth: depth of sampling (must be the same for each sample)
  #
  # Returns:
  #    a matrix
  otu.m = otu_table(mydata)
  otu.s = apply(otu.m, 2, function(sample) sort(sample, decreasing=T))
  otu.s = otu.s / depth
  otu.c = apply(otu.s, MARGIN=2, function(sample) cumsum(sample) )
  otu.c
}

CompareClassifications <- function(bioms, prefixes){
  # Calculates the number of OTUs classified at each tax level in a list of 
  # phyloseq objects. Used to compare the results of classification of the same
  # dataset using different reference taxonomies.
  # 
  # Calculates: taxa (number unique taxa)
  #             OTUs (total number OTUs, each taxa can have >1 OTUs)
  #             reads (number of total reads from OTUs classified at that level)
  #
  # Args:
  #   bioms: ta list of phyloseq objects
  #   prefixes: a vector of prefixes to append to the result col headings
  #
  # Returns:
  #    a data frame of results (colnames: prefix_taxa, prefix_OTUs, prefix_reads )
  tax.levels = c('kingdom', 'phylum', 'class', 'order', 
                 'family', 'genus', 'total')
  dataHeaders = c( outer(prefixes, c("taxa", "OTUs", "reads"),  
                         FUN = paste, sep = "_") )
  results <- data.frame(matrix(rep(0, 7*(3*length(bioms))), 
                               nrow=7, ncol=(3 * length(bioms)) ))
  rownames(results) = tax.levels
  colnames(results) = dataHeaders
  empty_strings <- c("k__", "p__", "c__", "o__", "f__", "g__") 
  
  for (n in 1:length(bioms)){
    biom = bioms[[n]]
    cols = dataHeaders[ seq(1, by=length(bioms), length.out=3) + (n - 1) ]
    
    # count OTUs with tax assignments
    v = apply(tax_table(biom), 2, function(tax.level) {
      x = ifelse(tax.level %in% empty_strings, NA, tax.level)
      sum(!(is.na(x)))
    })
    results[1:6, cols[2] ] = v
    
    for (level in 2:6) {
      tax.level = tax.levels[level]
      if (tax.level == "kingdom") {
        biom.m <- subset_taxa(biom, (kingdom != "k__"))
      }
      if (tax.level == "phylum") {
        biom.m <- subset_taxa(biom, (phylum != "p__"))
      }
      if (tax.level == "class") {
        biom.m <- subset_taxa(biom, (class != "c__"))
      }
      if (tax.level == "order") {
        biom.m <- subset_taxa(biom, (order != "o__"))
      }
      if (tax.level == "family") {
        biom.m <- subset_taxa(biom, (family != "f__"))
      }
      if (tax.level == "genus") {
        biom.m <- subset_taxa(biom, (genus != "g__"))
      }
      
      m = tax_glom(biom.m, tax.level)
      # count unique assignments
      results[tax.level, cols[1]] = ntaxa(m)
      # count reads with assignments 
      results[tax.level, cols[3]] = sum(sample_sums(m))
    }
    results[7, cols[1] ] = ntaxa(biom)
    results[1, cols[2] ] = ntaxa(biom)
    results[1, cols[3] ] = sum(sample_sums(biom))
  }
  results
  
}



topNTable <- function (mydata, topN, func, taxlevels=3:6) {
  if (topN > ntaxa(mydata)){
    topN <- ntaxa(mydata)
  }
  mydata.ordered <- names(sort(taxa_sums(mydata), decreasing=TRUE))
  cbind(tax_table(mydata)[ mydata.ordered[1:topN], taxlevels], 
        round(apply(otu_table(mydata)[mydata.ordered[1:topN], ], 1, func), 2))
}

plotTopN <- function(psdata, topN, plot_filename, core_cutoff=1){
  # boxplot of all data for the top 100 OTUs
  nsamp        <- nsamples(psdata)
  req_for_core <- ceiling(nsamp * core_cutoff)
  depth        <- as.numeric(sample_sums(psdata)[1])
  topnames     <- names(sort(apply(otu_table(psdata), 1, median), 
                             decreasing = TRUE )[1:topN])
  keep         <- taxa_names(psdata) %in% topnames
  topNtaxa     <- prune_taxa(x=psdata, taxa=keep )
  message(paste("sampled  top", topN, "OTUs", sep= " "))
  topN_fraction <- sum(sample_sums(topNtaxa)) / sum(sample_sums(psdata)) *100
  message(paste(round(topN_fraction, 1), "% of the total reads", sep= "") )
  
  #convert count to percent
  df.topN        <- data.frame(as.matrix(otu_table(topNtaxa) / depth * 100))
  num_pos_plants <- apply(df.topN, 1, function(x) sum(x > 0))
  df.topN[ df.topN == 0] <- 0.00025
  
  df.topN$median_abun  <- apply(df.topN[ , 1:nsamp], 1, median)
  df.topN$OTU          <- factor(row.names(df.topN))
  df.topN$n_pos_plants <- num_pos_plants 
  df.topN$core     <- factor(ifelse(df.topN$n_pos_plants >= req_for_core, 1, 2),
                             levels = c(1, 2), labels = c("core", "non-core" ))
  df.topN$taxname  <- makeTaxLabels(df.topN$OTU, topNtaxa)
  df.topN   <- df.topN[ order(df.topN$median_abun, decreasing=TRUE) , ]
  
  df.topN[ df.topN$OTU     == "681", "taxname"]                      <- "Tetrasphaera"
  df.topN[ df.topN$taxname == "CandidatusAccumulibacter", "taxname"] <- "Accumulibacter"
  df.topN[ df.topN$taxname == "CandidatusXenovorus", "taxname"]      <- "Xanthamonadaceae"
  df.topN[ df.topN$taxname == "CandidatusEpiflobacter", "taxname"]   <- "Epiflobacter"
  df.topN[ df.topN$taxname == "Unk01", "taxname"]                    <- "Ellin_Unk01"  
  df.topN[ df.topN$taxname == "Methyloversatilis", "taxname"]        <- "Sulfuritalia"
  df.topN[ df.topN$taxname == "Unk04", "taxname"]                    <- "Ellin_Unk04"
  
  df.topN.l <- melt(data = df.topN, 
                    id.vars = c("OTU", "taxname", "median_abun", "n_pos_plants", "core"), 
                    variable.name = "sample", 
                    value.name = "per_abun")
  
  pl <- ggplot(data= df.topN.l, aes(x = reorder(OTU, per_abun, order= TRUE,
                                                function(x) -median(x)), 
                                    y = per_abun) ) +
    scale_y_log10(   name = "Abundance (%)",
                     breaks = c(0.001, 0.01, 0.1, 1, 10), expand=c(0,0) ) + 
    scale_x_discrete(name   = "genus-level OTUs", 
                     labels = function(x) df.topN.l[ df.topN.l$OTU == x,
                                                     "taxname"]) +
    coord_cartesian(ylim = c(0.0025, 50)) +
    scale_fill_manual(values=c("white", "grey80")) +
    geom_boxplot( aes(fill = core) ) +
    theme_bw() +
    theme(legend.position="none", 
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8, vjust = 0.2), 
          axis.text.x  = element_text(size = 6, hjust = 1, vjust = 1, 
                                      angle = 45), 
          axis.text.y  = element_text(size = 8))
  
  if (plot_filename != FALSE){
    ggsave(filename= plot_filename, plot= pl, 
           width= 18, height= 12.8, units= "cm")    
  }
  return(pl)
} 

formatAverage <- function(mydata, decimals=1){
  # takes a vector and calculate and format an output string as: mean (+-SD, n=)
  results <- rep(0, 3)
  results[1] <- round(mean  (mydata), decimals)
  results[2] <- round(sd    (mydata), decimals)
  results[3] <- round(length(mydata), decimals)
  s <- (paste(
    results[1], " (+- ", results[2], ",", " n=", results[3], ")", sep=""))
  return(s)
}

makeOTUTable <- function (mydata, outfilename) {
  colheaders <- c(colnames(otu_table(mydata)), rank_names(mydata)[1:6])
  df.OTUtable <- data.frame(matrix(nrow=ntaxa(mydata), 
                                   ncol=(nsamples(mydata) + 6), 
                                   dimnames=list(taxa_names(mydata), 
                                                 colheaders)))
  df.OTUtable[ taxa_names(mydata), sample_names(mydata)] <- 
    otu_table(mydata)
  df.OTUtable[ taxa_names(mydata), rank_names(mydata)[1:6]] <- 
    tax_table(mydata)
  write.table(df.OTUtable, file=outfilename, sep="\t", na="NA", 
              row.names=TRUE, col.names=TRUE)
  df.OTUtable
}

calculate_cell_OM <- function(cell_volume, cell_density, OM_per_cell_WW){
  # calculation of cell fraction of organic matter
  # cell mass = 95 fg cell⁻¹     [3]
  # cell mass = 9.5e-13 g cell⁻¹
  cell_volume     <- cell_volume  * (1 / um3_per_cm3)
  # cm⁻³ cell⁻¹
  cell_WW         <- cell_density * cell_volume    
  # g-cell cell⁻¹    g-cell cm³-cell⁻¹ *  cm³-cell cell⁻¹  
  cell_OM         <- cell_WW * OM_per_cell_WW 
  # g-OM g-cell⁻¹    g-cell cell⁻¹ * g-OM g-cell⁻¹ 
  cell_OM  
}