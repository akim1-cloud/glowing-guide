### install a new package; you only need to do this once.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")

### libraries
library(ggplot2)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(2)
library(Rsamtools)

### this
### specify the bam file
fl <- system("ls -d /standard/BerglandTeach/mapping_output/*PRJNA315172*/*.original.bam", intern=T)

### How many reads map to each chromosome?
  rd <- foreach(bamFile=fl, .combine="rbind") %dopar% {
    ### tell me what file we are working on
      message(bamFile)
      
### get the information about number of reads for each chromosome
      if(!file.exists(paste(bamFile, ".bai", sep=""))) indexBam(bamFile)
      stats <- as.data.table(idxstatsBam(bamFile))
      
### we need to subset to the main autosomal arms of melanogaster and simulans
      setkey(stats, seqnames)
      stats_small <- stats[J(c("2L", "2R", "3L", "3R", "X"))]
      stats_small[seqnames == "X", chrType := "X"]
      stats_small[seqnames != "X", chrType := "Auto"]
      
### this command aggregates the data. For each unique value of species, we calculate the total number of mapped reads.
  stats_small.ag <- stats_small[,list(mapped=sum(mapped), seqlength=sum(seqlength)), list(chrType)]
 
 ###depth
 auto_depth <- stats_small.ag[chrType == "Auto"]$mapped / stats_small.ag[chrType == "Auto"]$seqlength
 x_depth <- stats_small.ag[chrType == "X"]$mapped / stats_small.ag[chrType == "X"]$seqlength
 
 ###ratio?
 x_auto_ratio <- x_depth / auto_depth
 propFemale <- pmax(0, pmin(1, 2 * (x_auto_ratio -0.5)))
 
  ### format output
  out <- data.table(
    x_depth = stats_small.ag[chrType=="X"]$mapped / stats_small.ag[chrType =="X"]$seqlength,
    auto_depth = stats_small.ag[chrType == "Auto"]$mapped / stats_small.ag[chrType == "Auto"]$seqlength
  )
  
  out[, x_auto_ratio := x_depth / auto_depth]
  out[, propFemale := pmax(0, pmin(1, 2 * (x_auto_ratio - 0.5)))]
  out[, samp := last(tstrsplit(bamFile, "/"))]
  
  ### return output
  return(out)
}

print(rd)

### estimated vs real
p <- ggplot(data=rd, aes(x=samp, y=propFemale, fill=propFemale)) +
    geom_col(color="black", width=0.6) +
    geom_hline(yintercept=1.0, linetype="dashed", color="red", linewidth=0.8) +
    geom_hline(yintercept=0.0, linetype="dashed", color="blue", linewidth=0.8) +
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.2)) +
    labs(
      x = "Sample",
      y = "Estimated Proportion Female",
      title = "Sex Ratio Estimation via X:Autosome Coverage Ratio",
      subtitle = "Red line = 100% Female (1.0), Blue line = 100% Male (0.0)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
      
ggsave("sex_ratio_figure.jpeg", plot=p, width=8, height=6, dpi=300, device="jpeg", type="cairo")