### ============================================================================
### Drosophila Genomic Differentiation & Outlier Enrichment Pipeline
### ============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)
library(doMC)
library(dplyr)

registerDoMC(4)

# 1. Collect and combine all chunk files from scratch directory
fl <- list.files("/scratch/uzm5ue/Fst_out", pattern = "\\.Rdata$", full.names = TRUE)

out_list <- lapply(fl, function(f) {
  env <- new.env()
  load(f, envir = env)
  if (exists("out", envir = env)) {
    dt <- get("out", envir = env)
    if (is.data.frame(dt) && nrow(dt) > 0) return(as.data.table(dt))
  }
  return(NULL)
})
out <- rbindlist(out_list, fill = TRUE)
setDT(out)
names(out) <- tolower(names(out))

# Save master combined object
save(out, file = "/scratch/uzm5ue/big_fstOutput.Rdata")

# 2. Standardize annotations into simplified functional classes
out[, annotation_simple := case_when(
  annotation == "5_prime_UTR_premature_start_codon_gain_variant" ~ "start_stop",
  annotation == "downstream_gene_variant" ~ "intergenic_region",
  annotation == "initiator_codon_variant" ~ "splice",
  annotation == "initiator_codon_variant&splice_region_variant" ~ "splice",
  annotation == "missense_variant&splice_region_variant" ~ "missense_variant",
  annotation == "splice_acceptor_variant&intron_variant" ~ "splice",
  annotation == "splice_donor_variant&intron_variant" ~ "splice",
  annotation == "splice_region_variant" ~ "splice",
  annotation == "splice_region_variant&intron_variant" ~ "splice",
  annotation == "splice_region_variant&non_coding_transcript_exon_variant" ~ "splice",
  annotation == "splice_region_variant&stop_retained_variant" ~ "splice",
  annotation == "splice_region_variant&synonymous_variant" ~ "splice",
  annotation == "start_lost&splice_region_variant" ~ "start_lost",
  annotation == "stop_gained&splice_region_variant" ~ "stop_gained",
  annotation == "stop_lost&splice_region_variant" ~ "stop_gained",
  annotation == "upstream_gene_variant" ~ "intergenic_region",
  annotation %in% c("stop_gained", "stop_lost", "stop_retained_variant", "start", "start_lost") ~ "start_stop",
  annotation == "start_lost" ~ "start_stop",
  annotation == "stop_gained" ~ "start_stop",
  annotation == "initiator_codon_variant&non_canonical_start_codon" ~ "start_stop",
  TRUE ~ annotation
)]
out[is.na(annotation_simple), annotation_simple := annotation]

# 3. Build sliding window dictionary (100kb window, 50kb step)
win.bp <- 1e5
step.bp <- 5e4
setkey(out, "chr")

wins <- foreach(chr.i = c("2L", "2R", "3L", "3R"), .combine = "rbind", .errorhandling = "remove") %dopar% {
  tmp <- out[J(chr.i)]
  data.table(chr = chr.i,
             start = seq(from = min(tmp$pos), to = max(tmp$pos) - win.bp, by = step.bp),
             end = seq(from = min(tmp$pos), to = max(tmp$pos) - win.bp, by = step.bp) + win.bp)
}
wins[, window := 1:nrow(wins)]

# 4. Compute Normalized Fst and Windowed Outlier Enrichment (Binomial P-values)
# 4. Compute Normalized Fst and Windowed Outlier Enrichment (Binomial P-values)
out.norm <- out[, list(normFst = 1 - rank(avefst) / (1 + length(avefst)), chr, pos, avefst), list(trt1, trt2)]
setkey(out.norm, chr)

windowedFst.pval_list <- foreach(i = 1:nrow(wins), .errorhandling = "remove") %dopar% {
  tmp <- out.norm[J(wins[i]$chr)][pos >= wins[i]$start & pos <= wins[i]$end]
  if (nrow(tmp) == 0) return(NULL)
  
  foreach(j = c(0.001, 0.005, 0.01, 0.05), .combine = "rbind") %do% {
    tmp.ag <- tmp[, list(
      nSNPs = .N, 
      TR = sum(normFst <= j, na.rm = TRUE),
      thr = j,
      avefst = mean(avefst, na.rm = TRUE),
      meanPos = mean(pos, na.rm = TRUE), 
      window = i
    ), list(trt1, trt2, chr)]
    tmp.ag[, pvalue := pbinom(TR, nSNPs, thr, lower.tail = FALSE)]
    tmp.ag
  }
}

# Remove any NULL elements from empty windows before combining
windowedFst.pval_list <- windowedFst.pval_list[!sapply(windowedFst.pval_list, is.null)]

# Now combine cleanly into a single data.table
windowedFst.pval <- rbindlist(windowedFst.pval_list)

# 5. Generate Individual Publication Plots
p2 <- ggplot(data = out[avefst > 0.15], aes(x = pos, y = avefst, color = trt1)) +
  geom_point(alpha = 0.7, size = 1.2) +
  facet_grid(trt1 + trt2 ~ chr, scales = "free_x") +
  labs(title = "Manhattan Plot of Extreme Fst Values (avefst > 0.15)", x = "Genomic Position", y = "Average Fst", color = "Treatment") +
  theme_minimal(base_size = 11) + theme(strip.text = element_text(face = "bold"), axis.text.x = element_blank(), legend.position = "bottom")

p3 <- ggplot(windowedFst.pval, aes(x = meanPos, y = -log10(pvalue), color = factor(thr), group = thr)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.6) +
  facet_grid(trt1 + trt2 ~ chr, scales = "free_x") +
  scale_color_viridis_d(name = "Quantile Threshold") +
  labs(title = "Windowed Outlier Enrichment Across Genome", x = "Genomic Position", y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))

p4 <- ggplot(data = out[annotation_simple %like% "missense" | annotation_simple %like% "synon"], 
             aes(x = annotation_simple, y = avefst, fill = annotation_simple)) + 
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  facet_grid(trt1 ~ trt2) + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Fst Comparison: Missense vs. Synonymous Variants", x = "Variant Annotation", y = "Average Fst", fill = "Annotation") +
  theme_minimal(base_size = 11) + theme(strip.text.x = element_blank(), legend.position = "bottom")

# 6. Combine Panels and Save Composite Figure
combined_plot <- (p2 / p3 / p4) + 
  plot_annotation(
    title = "Comprehensive Genomic Differentiation Analysis",
    subtitle = "A: Extreme Fst Manhattan Scan | B: Windowed Outlier Enrichment | C: Functional Variant Comparison",
    tag_levels = 'A',
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"))
  )

ggsave(filename = "/scratch/uzm5ue/Fst_out/combined_genomic_analysis.pdf", plot = combined_plot, width = 12, height = 16, device = "pdf", dpi = 300)
print("Pipeline execution and composite figure export complete!")