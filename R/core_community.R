plot_CumulRankAbundance <- function(dataset, core_prop=1, fname = NULL, 
                                    dominant_fraction, average = TRUE){
  # Plot Cumulative Rank abundance
  depth = as.numeric(sample_sums(dataset)[1])
  print(paste0("percent for count of 1: ", signif((1 / depth) * 100, 3)))
  
  # plot rank cumulative abundance --------------
  otu.m <- otu_table(dataset)
  otu.s <- apply(otu.m, 2, function(sample) sort(sample, decreasing=T))
  otu.s <- otu.s / depth
  otu.c <- apply(otu.s, MARGIN=2, function(sample) cumsum(sample) )
  
  res.mean <- apply(otu.c, MARGIN=1, mean )
  res.sd   <- apply(otu.c, MARGIN=1, sd )
  df.otu   <- data.frame("mean"=res.mean * 100, "sd"=res.sd * 100) 
  
  # where is the dominant cutoff?
  gm <- sum(res.mean < dominant_fraction)
  print(paste0( dominant_fraction * 100, "% cutoff: ", round(gm, 1)))
  
  otus.dom <- apply(otu.c, 2, function(x) sum(x < dominant_fraction))
  print(paste0("quantiles for #otus for ", dominant_fraction * 100, "% of reads"))
  print(quantile(otus.dom))
  otus_dom_abun <- otu.s[ otus.dom ]
  print(paste0("abundance of otu at ", dominant_fraction * 100, 
               "th percentile: ", formatAverage(otus_dom_abun * 100, 2 )))
  
  #cheat! #TODO
  cutoffs <- data.frame("otus" = gm, "reads"= c(dominant_fraction * 100))
  #  print(cutoffs)
  
  p1 <- ggplot(df.otu, aes(x = 1:length(mean), y = mean)) + 
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70") +
    geom_line() +
    geom_point(size = 1) +
    xlab(label= "Species-level OTUs in rank order") +
    ylab(label= "Cumulative read abundance (%)") +
    scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100),  expand=c(0,0)) +
    scale_x_log10(breaks=c(1, 10, 100, 1000),  limits=c(1, 1000), expand=c(0,0)) +
    geom_segment(data= cutoffs, aes(x = 0, y = reads, xend = otus, 
                                    yend = reads), size=0.3) +
    geom_segment(data= cutoffs, aes(x = otus, y = 0, xend = otus, 
                                    yend = reads), size=0.3) +
    geom_text(data= cutoffs, aes(x=otus, y= 5, label=round(otus, 0), size = 0.5)) +
    annotation_logticks(sides = "b", size = 0.3) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8, vjust = 0.2),
          axis.text.x  = element_text(size = 6),
          axis.text.y  = element_text(size = 6),
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color= "black"),
          axis.line        = element_line(size = 0.3),
          legend.position  = "none" )
  
  
  if (!(is.null(fname))){
    print(paste0("wrote plot to: ", fname))
    ggsave(plot=p1, file = fname, 
           width=8, height=6.37, units="cm", dpi=200)
  }
  
  #   print("%reads of top10 OTUs:")
  #   print(paste(signif(mean(otu.c[10,])*100, 2), 
  #               "(+-", signif(sd(otu.c[10,])*100, digits=2), ")"))
  #   
  #   print("%reads of top100 OTUs")
  #   print(paste(signif(mean(otu.c[100,])*100, digits=2), 
  #               "(+-", signif(sd(otu.c[100,])*100, digits=2), ")"))
  #   
  #   
  #   print("%reads of >0.1%")
  #   print(paste(signif(mean(otu.c[100,])*100, digits=2), 
  #               "(+-", signif(sd(otu.c[100,])*100, digits=2), ")"))
  
  return(p1)
}

calcDataPerID <- function(id, identities, datasets, frame, 
                          abun_cutoff=0, core_cutoff=1){

  data  <- datasets[[grep(pattern=id, identities)]]
  depth <- sample_sums(data)[1]
  cutoff_counts <- abun_cutoff * sample_sums(data)[1]
  if (cutoff_counts == 0) cutoff_counts <- 1
  data.binary <- transform_sample_counts(
    data, function(x) ifelse(x >= cutoff_counts, 1, 0 ) )
  
  # OTUs per core conservation
  sampletotals    <- sapply(datasets, function(dataset) nsamples(dataset))
  numsamples      <- sampletotals[1]
  core_low_limit      <- ceiling((numsamples * core_cutoff))
  core_groups         <- core_low_limit:numsamples
  num_core_samples    <- length(core_groups)
  num_noncore_samples <- numsamples - num_core_samples
  
  x <- as.character(frame[ (frame$identities == id), "numsamples"])
  num_core_taxa   <-   sum(taxa_sums(data.binary))
  taxsum_per_core <- table(taxa_sums(data.binary))
  frame[ frame$identities == id , "taxsum" ]  <- taxsum_per_core[x]
  
  # read abundance per core conservation
  reads_per_core    <- rowsum(otu_table(data), group=taxa_sums(data.binary))
  readprop_per_core <- reads_per_core / depth
  meanprop_per_core <- apply(readprop_per_core, 1, mean)
  meansdpercore     <- apply(readprop_per_core, 1, sd) 

  frame[ (frame$identities == id) , "readprop" ]    <- meanprop_per_core[x]
  frame[ (frame$identities == id) , "readpropSD" ]  <- meansdpercore[x]
  
  return(frame)
}

calcSummaryData <- function(datasets, identities, abun_cutoff=0, core_cutoff=1){
  # set up data frame
  sampletotals <- sapply(datasets, function(dataset) nsamples(dataset))
  if ((sum(sampletotals) %% sampletotals[1]) != 0 ){
    stop("num samples must be the same in seach dataset")
  }
  numsamples <- as.numeric(sampletotals[1])
  totalrows  <- as.numeric(numsamples * length(identities))
  
  frame <- data.frame("identities"= factor(rep(identities, each=numsamples ), 
                                           levels=sort(identities, 
                                                       decreasing=TRUE)),
                    "numsamples" = rep(1:numsamples, length(identities)),
                    "taxsum"     = rep(NA, totalrows),
                    "readprop"   = rep(NA, totalrows),
                    "readpropSD" = rep(NA, totalrows),
                    "corestatus" = gl(n=2, k= totalrows, length=totalrows, 
                                     labels=c("non-core", "core")) )
  
  core_low_limit <- ceiling((numsamples * core_cutoff))
  core_groups    <- core_low_limit:numsamples
  frame[ frame$numsamples %in% core_groups, "corestatus"] <- "core"
  
  for (id in identities){
    frame <- calcDataPerID(id, identities, datasets, frame, 
                           abun_cutoff, core_cutoff)
  }
  return(frame)
}

plotCore <- function(dframe, cumulative = TRUE){
  require(RColorBrewer)
  require(gridExtra)
  
  # TODO core cut off line
  
  otuplot <- ggplot(data=dframe, aes(numsamples, taxsum, fill=identities)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_colour_brewer(palette="Set1") +
    geom_bar(stat="identity", position= "dodge") +
    ylab("OTU count") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, vjust = 0.2),
          axis.text.x  = element_blank(),
          axis.text.y  = element_text(size = 10),
          legend.justification=c(1, 1), legend.position=c(1, 1))
  
  identities <- unique(dframe$identities)
  for (id in identities){
    rp <- rev(dframe[ dframe$identities == id, "readprop" ] * 100 )
    dframe[ dframe$identities == id, "cumsum" ] <- rev(cumsum(rp))
  }
  
  readsplot <- ggplot(data=dframe, aes(numsamples, readprop * 100, 
                                      fill= identities)) +
    geom_bar(stat="identity", position= "dodge") +
    ylab("Read abundance (% of total)") +
    xlab("Frequency (observed in n samples)") +
    scale_y_continuous(limits = c(0, 100), 
                       breaks = seq(0, 100, 20), 
                       expand=c(0,0)) +
    theme_bw() +
    theme(legend.position="none", 
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12, vjust = 0.2), 
          axis.text.x  = element_text(size = 10), 
          axis.text.y  = element_text(size = 10),
          legend.text  = element_text(size = 10), 
          legend.title = element_text(size = 12))

  if (cumulative == TRUE){
    readsplot <- readsplot + 
      geom_line(aes(y=cumsum, group= identities,color= identities), 
                size=0.5)
  }
  return( list(otuplot, readsplot))
}



calcCoreCommunity <- function(data, fig_filename, core_cutoff=0.9){
  library(RColorBrewer)
  
  data.binary <- transform_sample_counts(
    data, function(x) ifelse(x >= 1, 1, 0 ) )
  
  # OTUs per core conservation
  samplenum <-  nsamples(data)
  core_low_limit <- ceiling((samplenum * core_cutoff))
  core_groups    <- core_low_limit:samplenum
  num_core_samples    <- length(core_groups)
  num_noncore_samples <- samplenum - num_core_samples
  print(paste(samplenum, "samples" ) )
  print(paste(core_low_limit, "samples minimum to be core" ) )
  
  num_core_taxa <- sum(taxa_sums(data.binary) > core_low_limit )
  taxsum_per_core <- table(taxa_sums(data.binary))
  
  # read abundance per core conservation
  reads_per_core <- rowsum(otu_table(data), group=taxa_sums(data.binary))
  readsum_per_core <- apply(reads_per_core, 1, sum)
  readprop_per_core <- readsum_per_core / sum(readsum_per_core)
  print("percent readsum per core:")
  print(round(readprop_per_core * 100, 1))
  print(paste("core OTUs:", num_core_taxa ) )
  print(paste("percent reads in the core:", 
         round(sum(readprop_per_core[core_low_limit:samplenum]) * 100, 2) ) )  
  
  df.core <- data.frame(numsamples= 1:nsamples(data), 
                        taxsum= as.vector(taxsum_per_core), 
                        read_prop= readprop_per_core,
                        core= as.factor(c(rep("non-core", num_noncore_samples), 
                                          rep("core", num_core_samples)))  )
  
  png(file = fig_filename, width=700, height=500)
  
  pal  <- brewer.pal(4,"Set1")
  cols <- c(rep(pal[2], num_noncore_samples), rep(pal[1], num_core_samples))
  par(fig= c(0, 1, 0.5, 1))
  par(mar= c(0.5, 5, 2, 2))
  
  barplot(taxsum_per_core, xlab="", xaxt="n", cex.names=0.8,
          ylab="#OTUs", main="", col = cols)
  
  par(fig=c(0, 1, 0, 0.5), new= TRUE)
  par(mar=c(5, 5, 0.5, 2))
  plot(cumsum(readprop_per_core) * 100, type= "l",
       ylab= "cumul. % reads", main= "", xlab= "#samples" )
  abline(h= 100 - (sum(readprop_per_core[core_low_limit:samplenum]) * 100), 
         v= (core_low_limit - 1), col= pal[1])
  dev.off()

  df.core
}

calcCoreRep <- function (datasets, depth, ids, rep, stats=FALSE) {
  # wrapper around calcCoreRep to repeat the calculation a number of times.
  #
  # args: datasets: list of datasets
  #       depth:    sequence depth dataset at which data has been resampled
  #       ids:      vector of identities (names of the datasets)
  #       rep:      rep iterator
  # returns: result of 1 rep; data.frame in long format
  #
  coreDatasets        <- lapply(datasets, 
                                function(dataset) selectCoreDataset(dataset, depth))
  names(coreDatasets) <- paste("data.", ids, sep="")
  
  if (stats) {
    lapply(coreDatasets, function(dataset) printDatasetStats(dataset) )
  }
  
  coredataframe <- calcData(coreDatasets, ids, core_cutoff=1)
  
  reptotals <- ddply(coredataframe, 'identities', function(id) 
    c("depth"      = depth,
      "rep"        = rep,
      "totalOTUs"  = sum(id$taxsum),
      "coreOTUs"   = id$taxsum[ id$corestatus == "core"] , 
      "percentcoreOTUs"= round(id$taxsum[id$corestatus == "core"] / 
                                 sum(id$taxsum) *100, 1),
      "percentreads"   = round(id$readprop[id$corestatus == "core"] * 
                                 100, 1)))
  return(reptotals)
}

calcCoreReps <- function (datasets, ids, depth, reps=1, stats=FALSE) {
  # wrapper to resample at differnet depths and compare observed core.
  #
  # args: datasets: list of datasets
  #       depth:    vector of depths at which to resample
  #       ids:      vector of identities (names of the datasets)
  #       reps:     num reps (default=1)
  #       stats:    print stats? (default = FALSE)
  # returns: result of all depths/reps; data.frame in long format
  #
  OTUtotals <- lapply(seq_len(reps), function(rep) 
                                calcCoreRep(datasets, depth, ids, rep, stats ))  
  OTUtotals <- rbind.fill(OTUtotals)
  
  return(OTUtotals)
} 

get_ha_names <- function(samplecounts){
  counts         <- as.vector(samplecounts)
  sample         <- sampleNames(samplecounts)
  names(counts)  <- rownames(samplecounts)
  counts         <- sort(counts, decreasing=TRUE)
  percents       <- counts / sum(counts)
  scumsum        <- cumsum(percents)
  hanames        <- names(scumsum[scumsum < 0.8])
}
fill_ha_data <- function(dataset, samplename){
  ha_names <- get_ha_names(otu_table(dataset)[ , samplename])
  ifelse(taxa_names(dataset) %in% ha_names, 1, 0)
}

#Figure 4
plotFigure4 <- function(){
  plotlimits <- data.frame(nObs_range = c(1, 26), nHA_range = c(0, 26))
  
  nObs_breaks = c(1, 20, 25, 26)
  nHA_breaks  = c(0, 1, 10, 25, 26) 
  
  ggplot(plotlimits, aes(nObs_range, nHA_range) ) +
    scale_y_continuous(breaks= nHA_breaks,  labels= c(0, 1, 10, 25, 26), expand=c(0,0)) +
    scale_x_continuous(breaks= nObs_breaks, labels= c(1, 20, 25, 26), expand=c(0,0)) +
    geom_rect(aes(xmin=25, xmax=26, ymin=25, ymax=26), fill="red", color = "black") +
    xlab("Observed in n samples") +
    ylab("Highly-abundant in n samples") +
    geom_rect(aes(xmin=25, xmax=26, ymin=10, ymax=25), fill="orange", color = "black") +
    geom_rect(aes(xmin=20, xmax=25, ymin=10, ymax=25), fill="orange", color = "black") +
    geom_rect(aes(xmin=25, xmax=26, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
    geom_rect(aes(xmin=20, xmax=25, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
    geom_rect(aes(xmin= 1, xmax=20, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
    geom_rect(aes(xmin= 1, xmax=20, ymin= 0, ymax= 1), fill="grey", color = "black") +
    geom_rect(aes(xmin=20, xmax=25, ymin= 0, ymax= 1), fill="grey", color = "black") +
    geom_rect(aes(xmin=25, xmax=26, ymin= 0, ymax= 1), fill="grey", color = "black") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8, vjust = 0.2),
          axis.text.x  = element_text(size = 6),
          axis.text.y  = element_text(size = 6),
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color= "black"),
          axis.line        = element_line(size = 0.3),
          legend.position  = "none" )
}
#Figure 4
plotNewFigure4 <- function(mydata, alphafactor = 0.5){
  nObs_breaks = c(1, 20, 25, 26)
  nHA_breaks  = c(0, 1, 10, 25, 26) 
  library(scales)
  
  ggplot(mydata, aes(y = nHA, x = nObs, color = group) ) +
    geom_point(stat= "identity", size = I(3), alpha = alphafactor) +
    scale_y_continuous(breaks= nHA_breaks,  labels= c(0, 1, 10, 25, 26)) +
    scale_x_continuous(breaks= nObs_breaks, labels= c(1, 20, 25, 26)) +
#    geom_rect(aes(xmin=25.5, xmax=26.5, ymin=25.5, ymax=26.5), fill="red", color = "black", alpha = 0.2) +
    xlab("Observed in n samples") +
    ylab("Highly-abundant in n samples") +
#     geom_rect(aes(xmin=25, xmax=26, ymin=10, ymax=25), fill="orange", color = "black") +
#     geom_rect(aes(xmin=20, xmax=25, ymin=10, ymax=25), fill="orange", color = "black") +
#     geom_rect(aes(xmin=25, xmax=26, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
#     geom_rect(aes(xmin=20, xmax=25, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
#     geom_rect(aes(xmin= 1, xmax=20, ymin= 1, ymax=10), fill="darkgreen", color = "black") +
#     geom_rect(aes(xmin= 1, xmax=20, ymin= 0, ymax= 1), fill="grey", color = "black") +
#     geom_rect(aes(xmin=20, xmax=25, ymin= 0, ymax= 1), fill="grey", color = "black") +
#     geom_rect(aes(xmin=25, xmax=26, ymin= 0, ymax= 1), fill="grey", color = "black") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8, vjust = 0.2),
          axis.text.x  = element_text(size = 6),
          axis.text.y  = element_text(size = 6),
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color= "black"),
          axis.line        = element_line(size = 0.3),
          legend.position  = "none" )
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}