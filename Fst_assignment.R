### ============================================================================
### Optimized Candidate SNP Allele Frequency Trajectory Pipeline
### ============================================================================

library(SeqArray)
library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)
library(doMC)
library(poolfstat)

registerDoMC(10)

### Load utility function for data extraction
source("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/utils/misc/getData_function.R")

### Open GDS file
genofile <- seqOpen("/scratch/uzm5ue/dest.all.PoolSNP.001.50.22Sept2026_ExpEvo.norep.eff.reheader.gds")

### Load meta-data file
samps <- fread("https://raw.githubusercontent.com/biol4559-uva/CompEvoBio_modules/refs/heads/main/data/expEvo_2026_samplemetadata.csv")

### Common SNP.dt filtering
seqResetFilter(genofile)
snp.dt <- data.table(chr = seqGetData(genofile, "chromosome"),
                     pos = seqGetData(genofile, "position"),
                     nAlleles = seqGetData(genofile, "$num_allele"),
                     id = seqGetData(genofile, "variant.id"),
                     refAllele = seqGetData(genofile, "$ref"),
                     altAllele = seqGetData(genofile, "$alt"))

snp.dt <- snp.dt[nAlleles == 2] 
seqSetFilter(genofile, snp.dt$id)
snp.dt[, af := seqGetData(genofile, "annotation/info/AF")$data]

### Load Big Fst object
load(file = "/scratch/uzm5ue/big_fstOutput.Rdata")

### Define candidates from target window on chromosome 3L
candidates <- data.table(chr = "3L", pos = seq(7219979, 7319979))
setkey(candidates, chr, pos)
candidates <- merge(snp.dt, candidates, allow.cartesian = TRUE)

### Extract allele frequencies for all candidate SNPs simultaneously (Vectorized Batch)
freqs <- getData(snps = candidates[, c("chr", "pos", "id", "refAllele", "altAllele", "nAlleles"), with = FALSE], samples = samps)

### Plot experimental evolution trajectories for project ID PRJNA315172
ggplot(data = na.omit(freqs[grepl("PRJEB5713", sampleId), .(exp_rep, af_nEff, variant.id)]), 
       aes(x = exp_rep, y = af_nEff, group = variant.id)) +
  geom_line(alpha = 0.6, color = "steelblue") +
  theme_minimal() +
  labs(title = "Allele Frequency Trajectories across Experimental Replicates", x = "Experimental Replication", y = "Effective Allele Frequency")

### Aggregate mean frequencies across replicates
freqs.ag <- freqs[grepl("PRJNA315172", sampleId)][, list(p = mean(af_nEff, na.rm = TRUE)), list(exp_rep, variant.id)]
ggplot(data = na.omit(freqs.ag), aes(x = as.factor(exp_rep), y = p, group = variant.id)) +
  +     geom_line(alpha = 0.6, color = "darkgreen") +
  +     theme_minimal() +
  +     labs(title = "Mean Allele Frequency Trajectory", x = "Experimental Replication", y = "Mean Frequency")

### Merge complete sample metadata including continent and recommendation status
setkey(freqs, sampleId)
freqs <- merge(freqs, samps[, c("sampleId", "Recommendation", "continent", "lat")], by = "sampleId", all.x = TRUE)

### Plot latitudinal cline for passing European samples
ggplot(data = freqs[Recommendation == "Pass" & continent == "Europe"], 
       aes(x = as.numeric(as.character(lat)), y = af_nEff, group = variant.id)) + 
  geom_line(alpha = 0.4, color = "purple") +
  theme_minimal() +
  labs(title = "European Latitudinal Cline for Target Window", x = "Latitude", y = "Effective Allele Frequency")
