# Select the influent/effluent data

CalcAbsoluteCount <- function(prop, totalcount){ 
  # calcualte the absolute populations
  #
  # Args:
  #   prop: the proportional abundance of the OTU
  #   totalcount: the total bacterial count
  #
  # Returns:
  #   absolute abundance of the OTU   
  prop * totalcount 
}

CalcNetGrowthRate <- function(abscount_as, abscount_ww, yield, srt) {
  # calculate the net growth rate
  #
  # Args:
  #   abscount_as: count of cells inside system
  #   abscount_ww: count of cells into system
  #   yield: growth yield 
  #   srt: sludge retention time
  #
  # Returns:
  #   net growth rate
  (abscount_as - (abscount_ww * 1/yield )) / (abscount_as * srt)
}

plot_netgrowth <- function(value, mydata=df.netgrowth){
  sub.data <- mydata
  if (value != "all"){
    sub.data <-  mydata[ mydata$pair == value, ]    
  }
  
  myBreaks <- c(1,10)
  p <- ggplot(data = sub.data, aes(y = prop_as * 100, x = km, 
                                 color = WW_abun, label = taxlab)) +
    geom_point(size = 3) +
    geom_text(size = 3, color = "black", vjust = 1.7, hjust = 0.5 ) +
    xlab(expression(paste("Net growth rate (", d^-1, ")", sep=""))) +
    scale_y_log10(breaks = c(0.1, 1, 10),
                  limits = c(0.1, 15), expand = c(0,0),
                  name = 'Activated sludge abundance (%)') +
    annotation_logticks(sides = "l", size = 0.2) +
    scale_x_continuous(breaks = c(seq(-0.2, 0, 0.05), 0.03), 
                       limits = c(-0.22, 0.05), expand=c(0,0) ) +
    scale_colour_brewer(type = "seq", palette = "YlOrRd",
                        name = "Wastewater abund. (%)") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10, vjust = 0.2),
          axis.text.x  = element_text(size = 9),
          axis.text.y  = element_text(size = 9),
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.justification = c(0, 1), 
          legend.position      = c(0, 1)  )
  
  pu <- ggplot(data=sub.data, aes(y=prop_as * 100, x= km, 
                                  color = WW_abun, label = taxlab)) +
    geom_point(size=3) +
    #    geom_text(size = 3, color = "black", vjust= 1, hjust= 0 ) +
    xlab(expression(paste("Bounded net growth rate (", d^-1, ")", sep=""))) +
    scale_y_log10(breaks= c(0.1, 1, 10), 
                  limits= c(0.1, 15),
                  name= 'Activated sludge abundance (%)') +
    scale_colour_brewer(type="seq", palette="YlOrRd",
                        name= "WW % abund.") +
    theme(axis.title.x = element_text(size = 10)) + 
    theme(axis.title.y = element_text(size = 10, vjust = 0.2)) + 
    theme(axis.text.x  = element_text(size = 9)) + 
    theme(axis.text.y  = element_text(size = 9)) +
    theme(legend.text  = element_text(size = 9)) + 
    theme(legend.title = element_text(size = 10)) +
    theme(legend.justification=c(0, 1), legend.position=c(0, 1))
  
  ggsave(path="figs", 
         filename=paste("abun_vs_netgrowth_", value, ".pdf", sep=""),
         plot=p, width=18, height=10, units="cm")
  
  ggsave(path="figs", 
         filename=paste("abun_vs_netgrowth", value, "unlabelled.pdf", sep="_"),
         plot=pu, width=18, height=10, units="cm")  
  p
}

reshapeNetgrowthData <- function(mydata.k.p) {
  # reshape data into long format for ggplot
  # plant | pair | OTU# | prop_ww | prop_as
  #
  # Args:
  #   mydata.k.p: phyloseq object
  #
  # Returns:
  #   a data.frame 
  
  df_map  <- data.frame(plant  = sample_data(mydata.k.p)$plant, 
                        pair   = as.factor(sample_data(mydata.k.p)$pair),
                        type   = sample_data(mydata.k.p)$type,
                        sample = sample_names(mydata.k.p))
  
  reformData <- function(sampletype, data_bytype, df_map){
    df_type      <- data_bytype[[sampletype]]
    myData       <- as.data.frame(otu_table(df_type))
    myData$OTU   <- row.names(myData)
    myData       <- melt(myData, id.vars="OTU", 
                       variable.name="sample", value.name="prop")
    myData$prop[ is.na(myData$prop) ] <- 0
    myData$type  <- sampletype
    myData$plant <- with(df_map, plant[match(myData$sample, sample)])
    myData$pair  <- with(df_map, pair [match(myData$sample, sample)])
    return(myData)
  }
  
  sampletypes <-  c("WW", "AS")
  
  mydata.k.ww <- subset_samples(mydata.k.p, type == "WW")
  mydata.k.as <- subset_samples(mydata.k.p, type == "AS")
  data_bytype <- list(mydata.k.ww, mydata.k.as)
  names(data_bytype) <- sampletypes
  
  df_bytype    <- lapply(sampletypes, function(sampletype) 
                               reformData(sampletype, data_bytype, df_map))
  print(str(df_bytype))
  df_netgrowth <- rbind_all(df_bytype)
  
  df_netgrowth <- dcast(df_netgrowth, OTU + plant + pair ~ type, value.var="prop")
  names(df_netgrowth)[names(df_netgrowth)=="WW"] <- "prop_ww"
  names(df_netgrowth)[names(df_netgrowth)=="AS"] <- "prop_as"
  
  return(df_netgrowth)
}

calcNetGrowthData <- function (myData, seed = FALSE, k_lower_bound=FALSE) {
  # Function to calcualte the netgrowth rate (k) for each OTU by comparing the
  # the absolute values in the influent and sludge. Uses the mass balance.
  # Args:
  #   myData: phyloseq object
  #   seed:   random number seed
  #   k_lower_bound: if desired a lower bound can be set for k.
  #                  All k < k_lower_bound set to k_lower_bound
  #
  # Returns:
  #   a data.frame
  
  # counts as proportion
  data.k.p <- transform_sample_counts(physeq = myData, function(x) x / sum(x))
  print(table(sample_data(myData)$plant, sample_data(myData)$type))
  
  #constants
  totalcount_as <- 2.0E+15  # cells kg-biomass^-1 from Vollertsen et al. (2001) 
  totalcount_ww <- 3.4E+14  # cells kg-biomass^-1 
  yield         <- 0.28
  srt           <- 35
  
  # reshape data into long format for ggplot
  df.netgrowth <- reshapeNetgrowthData(data.k.p)
  df.netgrowth <- mutate(df.netgrowth, OTU = factor(OTU))
  
  # remove empty rows
  df.netgrowth <- df.netgrowth[ !((df.netgrowth$prop_as == 0) & 
                                  (df.netgrowth$prop_ww == 0)), ]
  
  # Calculate absolute counts from flow and conc.
  df.netgrowth$abscount_ww <- 
    CalcAbsoluteCount(df.netgrowth$prop_ww, totalcount_ww)
  df.netgrowth$abscount_as <- 
    CalcAbsoluteCount(df.netgrowth$prop_as, totalcount_as)
  
  df.netgrowth <- mutate(df.netgrowth,
          # Calculate apparent net growth rate
          k = (abscount_as - (abscount_ww * 1/yield )) / (abscount_as * srt),
          # Calculate ratio
          ratio = abscount_as * 1/srt / abscount_ww )
  
  # Set lower limit on negative growth rate 
  df.netgrowth$km <- df.netgrowth$k
  if (k_lower_bound != FALSE){
    df.netgrowth$km[df.netgrowth$km < k_lower_bound ] <- k_lower_bound    
  }
  
  #####################################################################
  # Classify OTUs based on whether they are detected in the influent/sludge only
  # or both. 
  df.netgrowth$class <- 2
  df.netgrowth$class[ df.netgrowth$prop_as == 0 ] <- 1
  df.netgrowth$class[ df.netgrowth$prop_ww == 0 ] <- 3
  df.netgrowth$class <- factor(df.netgrowth$class, levels=c(1, 2, 3), 
                               labels=c("WW", "Both", "AS"))
  
  df.netgrowth$km_class <- cut(df.netgrowth$km, 
                               breaks=c(-1, -0.04, -0.02, 0, 0.02, 0.03), 
                               labels=c("min", "<-0.02", "<0", "<0.02", "max"))
  df.netgrowth$AS_abun <- cut(df.netgrowth$prop_as, 
                              breaks=c(-1, 1e-05, 0.001, 0.005, 0.01, 0.1, 0.2), 
                              labels=c("WW only", "<0.1%", "0.1-0.5%", "0.5-1%.", "1-10%", "10-20%"))
  df.netgrowth$WW_abun <- cut(df.netgrowth$prop_ww, 
                              breaks=c(-1, 1e-05, 0.001, 0.005, 0.01, 0.1, 0.2), 
                              labels=c("AS only", "<0.1%", "0.1-0.5%", "0.5-1%.", "1-10%", "10-20%"))
  return(df.netgrowth)
}

printNetGrowthStats <- function (df.netgrowth) {
  print("Wastewater totals")
  
  print("What proportion of the WW OTUs were found in the sludge?")
  class_totals  <- (table(df.netgrowth$class))
  print("Number of OTUs in Influent/sludge:")
  print(class_totals)
  
  ww_reads <- as.vector(c(class_totals[1:2], sum(class_totals[1:2])))
  names(ww_reads) <- c(names(class_totals[1:2]), "Total")
  print("Number of OTUs in Influent:")
  print(ww_reads)
  print("Number of OTUs in Influent (%):")
  print(round(ww_reads[1:2] / sum(class_totals[1:2]) * 100), 1)
  
  ww_reads_by_class <- with(df.netgrowth,
                            aggregate(prop_ww * 100, 
                                      by=list(pair, class), 
                                      FUN=sum))
  names(ww_reads_by_class) <- c("pair", "class", "percent")
  ww_reads_by_class$percent <- round(ww_reads_by_class$percent, 1)
  print("WW reads per class:")
  print(ww_reads_by_class[ww_reads_by_class$class %in% c("WW", "Both"),])
  
  ww_both <- ww_reads_by_class[ ww_reads_by_class$class == "Both", "percent"]
  print("%reads from WW OTUs detected in sludge:")
  print(formatAverage(ww_both))
  ww_only <- ww_reads_by_class[ ww_reads_by_class$class == "WW", "percent"]
  print("%reads from WW OTUs not detected in sludge:")
  print(formatAverage(ww_only))
  print(summary(ww_only))
  
  print("########################")
  print("Activated sludge totals")
  as_reads <- as.vector(c(class_totals[2:3], sum(class_totals[2:3])))
  names(as_reads) <- c(names(class_totals[2:3]), "Total")
  print("Number of OTUs in sludge:")
  print(as_reads[c(2, 1, 3)])
  print("OTUs in sludge (%):")
  print(round(as_reads[c(2, 1)] / sum(class_totals[2:3]) * 100), 1)
  
  as_reads_by_class <- with(df.netgrowth,
                      aggregate(prop_as * 100, by=list(pair, class), FUN=sum))
  names(as_reads_by_class) <- c("pair", "class", "percent")
  as_reads_by_class$percent <- round(as_reads_by_class$percent, 1)
  print("AS reads per class:")
  print(as_reads_by_class[as_reads_by_class$class %in% c("AS", "Both"),])
  
  as_both <- as_reads_by_class[ as_reads_by_class$class == "Both", "percent"]
  print("%reads from AS OTUs detected in influent:")
  print(formatAverage(as_both))
}

plotNetGrowthDistribution <- function(df.netgrowth, filename = FALSE){
  # sludge cum abundance vs. net growth rate plot
  a <- ggplot(data=df.netgrowth, aes(x = km, fill = class)) +
    geom_histogram(binwidth = 0.005) +
    coord_cartesian(ylim=c(0, 500)) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = NULL, y = 'OTU count') +
    theme(legend.position="none")
  
  df.netgrowth  <- df.netgrowth[ order(df.netgrowth$km, decreasing = TRUE),]
  prop_at_zerok <- cumsum(df.netgrowth$prop_as / 6)[sum(df.netgrowth$km > 0)]
  df.netgrowth  <- df.netgrowth[ order(df.netgrowth$km, decreasing = FALSE),]
  
  b <- ggplot(data=df.netgrowth, 
              aes(x = km, y = cumsum(prop_as / 6) * 100)) +
    scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0,0)) +
    scale_x_continuous(breaks = c(seq(-0.2, 0, 0.05), 0.029), expand = c(0,0)) +
#     geom_rect(aes(xmin =  0,   xmax =  0.029, ymin = 0, ymax = 100), fill= "#d9d9d9", alpha = 0.1) +
#     geom_rect(aes(xmin = -0.1, xmax =  0,     ymin = 0, ymax = 100), fill= "#f0f0f0", alpha = 0.1) +
#     geom_rect(aes(xmin = -0.2, xmax = -0.1,   ymin = 0, ymax = 100), fill= "#ffffff", alpha = 0.1) +
    geom_step() +
    xlab(expression(paste("Net growth rate (", d^-1, ")", sep=""))) +
    ylab('Cumulative read abundance (%)') +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8, vjust = 0.2),
          axis.text.x  = element_text(size = 6),
          axis.text.y  = element_text(size = 6),
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          axis.line        = element_line(size = 0.3),
          legend.position  = "none" )
  
  
  #    geom_hline(yintercept = (1 - prop_at_zerok) * 100, color = "red") +
  #     scale_x_continuous(limit  = c(-0.041, 0.03), 
  #                        breaks = seq(-0.04, 0.03, 0.01)) +
  
  print("Percent of reads at k < zero:")
  print(round(100 - (prop_at_zerok * 100), 2) )
  
  return(list(a, b))
}

compareTop20AsWw <- function(myData){
  data.k.p      <- transform_sample_counts(myData, function(x) x / sum(x))
  data.k.ww     <- subset_samples(data.k.p, type == "WW")
  top20ww_otus  <- names(sort(taxa_sums(data.k.ww), decreasing=TRUE)[1:20])
  top20ww_names <- sapply(top20ww_otus, makeTaxLabel, data.k)
  
  data.k.as     <- subset_samples(data.k.p, type == "AS")
  top20as_otus  <- names(sort(taxa_sums(data.k.as), decreasing=TRUE)[1:20])
  top20as_names <- sapply(top20as_otus, makeTaxLabel, data.k)
  
  top_otus <- prune_taxa(x = data.k.p, 
                         taxa = unique(c(top20ww_otus, top20as_otus)))
  # transforms to percentage and divides by 6 to give average percent over the 
  # samplesfor using the "identity" stat in the next graph
  top_otus_6 <- transformSampleCounts(physeq=top_otus, function(x) x *100) 
  
  p3 <- plot_bar(top_otus_6, "Family", fill = "Genus", 
                 facet_grid = "pair~sample_type" ) +
    geom_bar(aes(fill = Genus, color = Genus), 
                stat = "identity", position = "stack") + 
    ylab("percent abundance") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10)) + 
    theme(axis.title.y = element_text(size = 10, vjust = 0.2)) + 
    theme(axis.text.x  = element_text(size = 7, hjust = 1, vjust = 1, 
                                      angle = 45, color= "black")) + 
    theme(axis.text.y  = element_text(size = 7, color= "black")) +
    theme(legend.text  = element_text(size = 8)) + 
    theme(legend.key   = element_rect(size = 0.3)) +
    theme(legend.title = element_text(size = 10),
          strip.text   = element_text(size = 7),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(linetype= "blank", fill = "white"),
          axis.line        = element_line(size = 0.3)) 
  
  return(p3)
}
