###############################################
## Supplementary/ Extended data figure 1     ##
###############################################

# Requires running of Figure_1E_GWAS_Contrast.R script prior to this!!!!

# Filter out specific chromosomes

filtered_data_CHR10 <- subset(AllSNP_GWAS, CHR == 10)
filtered_data_CHR15 <- subset(AllSNP_GWAS, CHR == 15)
filtered_data_CHR17 <- subset(AllSNP_GWAS, CHR == 17)
filtered_data_CHR19 <- subset(AllSNP_GWAS, CHR == 19)


###########################
# Plot Chr 10            ##
###########################

filtered_data_CHR10


library(rtracklayer)

# Read GTF file
cluhar_v2.0.2_gtf <- import("/pathtodata/Clupea_harengus.Ch_v2.0.2.108.gtf")

# Filter annotations for chromosome 10
chr10_annotations <- subset(cluhar_v2.0.2_gtf, seqnames(cluhar_v2.0.2_gtf) == "10")

# Convert filtered annotations to GRanges object
cluhar_v2.0.2_gtf_gr <- as(chr10_annotations, "GRanges")

#########
#########
peak <- 22428922
start_pos <- peak - 1000000
end_pos <- peak + 1000000

filtered_data_Chr_10_SIGN <- subset(filtered_data_CHR10, BONF <= 0.05)
Start_Sign <- 22021875
End_Sign <- 22480344
Size_Sign <- (End_Sign - Start_Sign)/ 1000000
print(Size_Sign)

filtered_data_Chr_10_peak <- subset(filtered_data_CHR10, POS >= Start_Sign & POS <= End_Sign)


#########
# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_10_core <- filtered_data_CHR10[, c("CHR", "POS")]

# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_10_peak <- filtered_data_Chr_10_peak[, c("CHR", "POS")]

#########
# Convert snp_chr_pos and gtf_data into GRanges objects
Chr_10_core_GR <- GRanges(seqnames = 10,
                          ranges = IRanges(start = Start_Sign, end = End_Sign),
                          strand = "*")

Chr_10_peak_GR <- GRanges(seqnames = filtered_data_Chr_10_peak$CHR,
                          ranges = IRanges(start = filtered_data_Chr_10_peak$POS, end = filtered_data_Chr_10_peak$POS),
                          strand = "*")

# Find overlaps between SNPs and genomic features
Chr_10_peak_gtf_hits <- findOverlaps(Chr_10_core_GR, cluhar_v2.0.2_gtf_gr, maxgap = 5e4)
Chr_10_peak_core_hits <- findOverlaps(Chr_10_core_GR, cluhar_v2.0.2_gtf_gr)

#########
Chr_10_peak_gtf <- cluhar_v2.0.2_gtf_gr[Chr_10_peak_gtf_hits@to]
Chr_10_core_gtf <- cluhar_v2.0.2_gtf_gr[Chr_10_peak_core_hits@to]
Chr_10_peak_gtf$peak <- Chr_10_core_GR[Chr_10_peak_gtf_hits@from]
Chr_10_peak_gtf$core <- F
Chr_10_peak_gtf$core[Chr_10_peak_gtf$gene_id %in% Chr_10_core_gtf$gene_id] <- T
Chr_10_peak_genes <- Chr_10_peak_gtf[Chr_10_peak_gtf$type == "gene"]
#########

#########
library(biomaRt)
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
mart = useEnsembl(biomart = "ensembl", dataset = "charengus_gene_ensembl")
listAttributes(mart = mart)

#Chr_10_peak_3_BM <- getBM(attributes = c("description", "go_id", "name_1006", "external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_10_peak_genes$gene_id)
Chr_10_peak_3_BM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_10_peak_genes$gene_id)
Chr_10_peak_genes$gene_name <- Chr_10_peak_3_BM$external_gene_name[match(Chr_10_peak_genes$gene_id, Chr_10_peak_3_BM$ensembl_gene_id)]
Chr_10_peak_gtf$gene_name <- Chr_10_peak_3_BM$external_gene_name[match(Chr_10_peak_gtf$gene_id, Chr_10_peak_3_BM$ensembl_gene_id)]

write.table(Chr_10_peak_genes, "/pathtoexport/Chr10_Genomic_Ranges_Annotations.txt")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 5)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(0, 60), main = "")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 35, y1 = 35, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      #if(genes_gtf$core[gene_idx]){
      #  text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      #}
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

# Usage:
DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR10,
                        #site_filter = site_filter,
                        peak_GR = Chr_10_peak_GR, 
                        full_gtf = Chr_10_peak_gtf, 
                        plot_reg = c(min(filtered_data_CHR10$POS) - 1e5, 
                                     max(filtered_data_CHR10$POS) + 1e5), 
                        pdf_file = "Chr_10_SNPChip_plot.pdf")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

#Zoom in - GENE TRACKS

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(29, 70), xlim = c(Start_Sign, End_Sign),  main = "",  yaxt = "n")
  
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(30, 70), xlim = c(Start_Sign, End_Sign), main = "", yaxt = "n")
  #axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  
  
  # Plot SNPs
  #points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
  #       col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  #points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 25, y1 = 25, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 48.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      if(genes_gtf$core[gene_idx]){
        text(x = mid(genes_gtf[target_gene_entry]), srt=45, cex=0.7, y = c(40, 61)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      }
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR10,
                        #site_filter = site_filter,
                        peak_GR = Chr_10_peak_GR, 
                        full_gtf = Chr_10_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_10_SNPChip_zoomed_plot_GENE_TRACK.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(0, 20), xlim = c(Start_Sign, End_Sign), main = "", xaxt = "n")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 42, y1 = 42, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR10,
                        #site_filter = site_filter,
                        peak_GR = Chr_10_peak_GR, 
                        full_gtf = Chr_10_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_10_SNPChip_zoomed_plot_CONTRAST.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 20), , main = "", xaxt = "n")
  axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  mtext("Position (Mb)", side = 3, line = 2)
  
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 18, y1 = 18, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR10,
                        #site_filter = site_filter,
                        peak_GR = Chr_10_peak_GR, 
                        full_gtf = Chr_10_peak_gtf, 
                        plot_reg = c(start_pos, 
                                     max(filtered_data_CHR10$POS) + 1e3), 
                        pdf_file = "Chr_10_SNPChip_zoomed_plot_CONTRAST_WIDE.pdf")




###########################
# Plot Chr 15            ##
###########################

filtered_data_CHR15

library(rtracklayer)

# Read GTF file
cluhar_v2.0.2_gtf <- import("/pathtodata/Clupea_harengus.Ch_v2.0.2.108.gtf")

# Filter annotations for chromosome 15
Chr_15_annotations <- subset(cluhar_v2.0.2_gtf, seqnames(cluhar_v2.0.2_gtf) == "15")

# Convert filtered annotations to GRanges object
cluhar_v2.0.2_gtf_gr <- as(Chr_15_annotations, "GRanges")

#########
peak <- 9023505
start_pos <- peak - 1000000
end_pos <- peak + 1000000

filtered_data_Chr_15_SIGN <- subset(filtered_data_CHR15, BONF <= 0.05)
Start_Sign <- 8813384
End_Sign <- 9323382
Size_Sign <- (End_Sign - Start_Sign)/ 1000000
print(Size_Sign)

filtered_data_Chr_15_peak <- subset(filtered_data_CHR15, POS >= Start_Sign & POS <= End_Sign)

#########
# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_15_core <- filtered_data_CHR15[, c("CHR", "POS")]

# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_15_peak <- filtered_data_Chr_15_peak[, c("CHR", "POS")]

#########
# Convert snp_chr_pos and gtf_data into GRanges objects
Chr_15_core_GR <- GRanges(seqnames = 15,
                          ranges = IRanges(start = Start_Sign, end = End_Sign),
                          strand = "*")

Chr_15_peak_GR <- GRanges(seqnames = filtered_data_Chr_15_peak$CHR,
                          ranges = IRanges(start = filtered_data_Chr_15_peak$POS, end = filtered_data_Chr_15_peak$POS),
                          strand = "*")

# Find overlaps between SNPs and genomic features
Chr_15_peak_gtf_hits <- findOverlaps(Chr_15_core_GR, cluhar_v2.0.2_gtf_gr, maxgap = 5e4)
Chr_15_peak_core_hits <- findOverlaps(Chr_15_core_GR, cluhar_v2.0.2_gtf_gr)

#########
Chr_15_peak_gtf <- cluhar_v2.0.2_gtf_gr[Chr_15_peak_gtf_hits@to]
Chr_15_core_gtf <- cluhar_v2.0.2_gtf_gr[Chr_15_peak_core_hits@to]
Chr_15_peak_gtf$peak <- Chr_15_core_GR[Chr_15_peak_gtf_hits@from]
Chr_15_peak_gtf$core <- F
Chr_15_peak_gtf$core[Chr_15_peak_gtf$gene_id %in% Chr_15_core_gtf$gene_id] <- T
Chr_15_peak_genes <- Chr_15_peak_gtf[Chr_15_peak_gtf$type == "gene"]
#########

#########
library(biomaRt)
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
mart = useEnsembl(biomart = "ensembl", dataset = "charengus_gene_ensembl")
listAttributes(mart = mart)

#Chr_15_peak_3_BM <- getBM(attributes = c("description", "go_id", "name_1006", "external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_15_peak_genes$gene_id)
Chr_15_peak_3_BM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_15_peak_genes$gene_id)
Chr_15_peak_genes$gene_name <- Chr_15_peak_3_BM$external_gene_name[match(Chr_15_peak_genes$gene_id, Chr_15_peak_3_BM$ensembl_gene_id)]
Chr_15_peak_gtf$gene_name <- Chr_15_peak_3_BM$external_gene_name[match(Chr_15_peak_gtf$gene_id, Chr_15_peak_3_BM$ensembl_gene_id)]

write.table(Chr_15_peak_genes, "/pathtoexport/Chr15_Genomic_Ranges_Annotations.txt")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 5)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(0, 60), main = "")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 35, y1 = 35, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      #if(genes_gtf$core[gene_idx]){
      #  text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      #}
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

# Usage:
DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR15,
                        #site_filter = site_filter,
                        peak_GR = Chr_15_peak_GR, 
                        full_gtf = Chr_15_peak_gtf, 
                        plot_reg = c(min(filtered_data_CHR15$POS) - 1e5, 
                                     max(filtered_data_CHR15$POS) + 1e5), 
                        pdf_file = "Chr_15_SNPChip_plot.pdf")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

#Zoom in - GENE TRACKS

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(29, 70), xlim = c(Start_Sign, End_Sign),  main = "",  yaxt = "n")
  
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(30, 70), xlim = c(Start_Sign, End_Sign), main = "", yaxt = "n")
  #axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  
  # Plot SNPs
  #points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
  #       col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  #points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 25, y1 = 25, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 48.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      if(genes_gtf$core[gene_idx]){
        text(x = mid(genes_gtf[target_gene_entry]), srt=45, cex=0.7, y = c(40, 61)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      }
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR15,
                        #site_filter = site_filter,
                        peak_GR = Chr_15_peak_GR, 
                        full_gtf = Chr_15_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_15_SNPChip_zoomed_plot_GENE_TRACK.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(0, 43), xlim = c(Start_Sign, End_Sign), main = "", xaxt = "n")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 42, y1 = 42, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR15,
                        #site_filter = site_filter,
                        peak_GR = Chr_15_peak_GR, 
                        full_gtf = Chr_15_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_15_SNPChip_zoomed_plot_CONTRAST.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), , main = "", xaxt = "n")
  axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  mtext("Position (Mb)", side = 3, line = 2)
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 42, y1 = 42, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR15,
                        #site_filter = site_filter,
                        peak_GR = Chr_15_peak_GR, 
                        full_gtf = Chr_15_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_15_SNPChip_zoomed_plot_CONTRAST_WIDE.pdf")





###########################
# Plot Chr 17            ##
###########################


filtered_data_CHR17


library(rtracklayer)

# Read GTF file
cluhar_v2.0.2_gtf <- import("/pathtodata/Clupea_harengus.Ch_v2.0.2.108.gtf")

# Filter annotations for chromosome 17
chr17_annotations <- subset(cluhar_v2.0.2_gtf, seqnames(cluhar_v2.0.2_gtf) == "17")

# Convert filtered annotations to GRanges object
cluhar_v2.0.2_gtf_gr <- as(chr17_annotations, "GRanges")

#########
peak <- 27080324
start_pos <- peak - 1500000
end_pos <- peak + 1000000

filtered_data_Chr_17_SIGN <- subset(filtered_data_CHR17, BONF <= 0.05)
Start_Sign <- 25814643
End_Sign <- 27549527
Size_Sign <- (End_Sign - Start_Sign)/ 1000000
print(Size_Sign)

filtered_data_CHR17_peak <- subset(filtered_data_CHR17, POS >= Start_Sign & POS <= End_Sign)

#########
# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_17_core <- filtered_data_CHR17[, c("CHR", "POS")]

# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_17_peak <- filtered_data_CHR17_peak[, c("CHR", "POS")]

#########
# Convert snp_chr_pos and gtf_data into GRanges objects
Chr_17_core_GR <- GRanges(seqnames = 17,
                          ranges = IRanges(start = Start_Sign, end = End_Sign),
                          strand = "*")

Chr_17_peak_GR <- GRanges(seqnames = filtered_data_Chr_17_peak$CHR,
                          ranges = IRanges(start = filtered_data_Chr_17_peak$POS, end = filtered_data_Chr_17_peak$POS),
                          strand = "*")

# Find overlaps between SNPs and genomic features
Chr_17_peak_gtf_hits <- findOverlaps(Chr_17_core_GR, cluhar_v2.0.2_gtf_gr, maxgap = 5e4)
Chr_17_peak_core_hits <- findOverlaps(Chr_17_core_GR, cluhar_v2.0.2_gtf_gr)

#########
Chr_17_peak_gtf <- cluhar_v2.0.2_gtf_gr[Chr_17_peak_gtf_hits@to]
Chr_17_core_gtf <- cluhar_v2.0.2_gtf_gr[Chr_17_peak_core_hits@to]
Chr_17_peak_gtf$peak <- Chr_17_core_GR[Chr_17_peak_gtf_hits@from]
Chr_17_peak_gtf$core <- F
Chr_17_peak_gtf$core[Chr_17_peak_gtf$gene_id %in% Chr_17_core_gtf$gene_id] <- T
Chr_17_peak_genes <- Chr_17_peak_gtf[Chr_17_peak_gtf$type == "gene"]
#########

#########
library(biomaRt)
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
mart = useEnsembl(biomart = "ensembl", dataset = "charengus_gene_ensembl")
listAttributes(mart = mart)

#Chr_17_peak_3_BM <- getBM(attributes = c("description", "go_id", "name_1006", "external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_17_peak_genes$gene_id)
Chr_17_peak_3_BM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_17_peak_genes$gene_id)
Chr_17_peak_genes$gene_name <- Chr_17_peak_3_BM$external_gene_name[match(Chr_17_peak_genes$gene_id, Chr_17_peak_3_BM$ensembl_gene_id)]
Chr_17_peak_gtf$gene_name <- Chr_17_peak_3_BM$external_gene_name[match(Chr_17_peak_gtf$gene_id, Chr_17_peak_3_BM$ensembl_gene_id)]

write.table(Chr_17_peak_genes, "/pathtoexport/Chr17_Genomic_Ranges_Annotations.txt")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$
#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 5)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(0, 60), main = "")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 35, y1 = 35, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      #if(genes_gtf$core[gene_idx]){
      #  text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      #}
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

# Usage:
DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR17,
                        #site_filter = site_filter,
                        peak_GR = Chr_17_peak_GR, 
                        full_gtf = Chr_17_peak_gtf, 
                        plot_reg = c(min(filtered_data_CHR17$POS) - 1e5, 
                                     max(filtered_data_CHR17$POS) + 1e5), 
                        pdf_file = "Chr_17_SNPChip_plot.pdf")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

#Zoom in - GENE TRACKS

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(29, 70), xlim = c(Start_Sign, End_Sign),  main = "",  yaxt = "n")
  
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(30, 70), xlim = c(Start_Sign, End_Sign), main = "", yaxt = "n")
  #axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  
  # Plot SNPs
  #points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
  #       col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  #points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 25, y1 = 25, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 48.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      if(genes_gtf$core[gene_idx]){
        text(x = mid(genes_gtf[target_gene_entry]), srt=45, cex=0.7, y = c(40, 61)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      }
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR17,
                        #site_filter = site_filter,
                        peak_GR = Chr_17_peak_GR, 
                        full_gtf = Chr_17_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_17_SNPChip_zoomed_plot_GENE_TRACK.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(0, 25), xlim = c(Start_Sign, End_Sign), main = "", xaxt = "n")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 42, y1 = 42, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR17,
                        #site_filter = site_filter,
                        peak_GR = Chr_17_peak_GR, 
                        full_gtf = Chr_17_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_17_SNPChip_zoomed_plot_CONTRAST.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 25), , main = "", xaxt = "n")
  axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  mtext("Position (Mb)", side = 3, line = 2)
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 21, y1 = 21, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR17,
                        #site_filter = site_filter,
                        peak_GR = Chr_17_peak_GR, 
                        full_gtf = Chr_17_peak_gtf, 
                        plot_reg = c(start_pos, 
                                     max(filtered_data_CHR17$POS) + 1e3), 
                        pdf_file = "Chr_17_SNPChip_zoomed_plot_CONTRAST_WIDE.pdf")




###########################
# Plot Chr 19            ##
###########################

filtered_data_CHR19

library(rtracklayer)

# Read GTF file
cluhar_v2.0.2_gtf <- import("/pathtodata/Clupea_harengus.Ch_v2.0.2.108.gtf")

# Filter annotations for chromosome 19
chr19_annotations <- subset(cluhar_v2.0.2_gtf, seqnames(cluhar_v2.0.2_gtf) == "19")

# Convert filtered annotations to GRanges object
cluhar_v2.0.2_gtf_gr <- as(chr19_annotations, "GRanges")

#########
peak <- 20543808
start_pos <- peak - 1000000
end_pos <- peak + 1000000


filtered_data_Chr_19_SIGN <- subset(filtered_data_CHR19, BONF <= 0.05)
Start_Sign <- 20422070
End_Sign <- 20622999
Size_Sign <- (End_Sign - Start_Sign)/ 1000000
print(Size_Sign)

filtered_data_Chr_19_peak <- subset(filtered_data_CHR19, POS >= Start_Sign & POS <= End_Sign)

#########
# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_19_core <- filtered_data_CHR19[, c("CHR", "POS")]

# Extract chromosome and position information from the filtered SNP data
filtered_data_Chr_19_peak <- filtered_data_CHR19_peak[, c("CHR", "POS")]

#########
# Convert snp_chr_pos and gtf_data into GRanges objects
Chr_19_core_GR <- GRanges(seqnames = 19,
                          ranges = IRanges(start = Start_Sign, end = End_Sign),
                          strand = "*")

Chr_19_peak_GR <- GRanges(seqnames = filtered_data_Chr_19_peak$CHR,
                          ranges = IRanges(start = filtered_data_Chr_19_peak$POS, end = filtered_data_Chr_19_peak$POS),
                          strand = "*")

# Find overlaps between SNPs and genomic features
Chr_19_peak_gtf_hits <- findOverlaps(Chr_19_core_GR, cluhar_v2.0.2_gtf_gr, maxgap = 5e4)
Chr_19_peak_core_hits <- findOverlaps(Chr_19_core_GR, cluhar_v2.0.2_gtf_gr)

#########
Chr_19_peak_gtf <- cluhar_v2.0.2_gtf_gr[Chr_19_peak_gtf_hits@to]
Chr_19_core_gtf <- cluhar_v2.0.2_gtf_gr[Chr_19_peak_core_hits@to]
Chr_19_peak_gtf$peak <- Chr_19_core_GR[Chr_19_peak_gtf_hits@from]
Chr_19_peak_gtf$core <- F
Chr_19_peak_gtf$core[Chr_19_peak_gtf$gene_id %in% Chr_19_core_gtf$gene_id] <- T
Chr_19_peak_genes <- Chr_19_peak_gtf[Chr_19_peak_gtf$type == "gene"]
#########

#########
library(biomaRt)
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
mart = useEnsembl(biomart = "ensembl", dataset = "charengus_gene_ensembl")
listAttributes(mart = mart)

#Chr_19_peak_3_BM <- getBM(attributes = c("description", "go_id", "name_1006", "external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_19_peak_genes$gene_id)
Chr_19_peak_3_BM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = Chr_19_peak_genes$gene_id)
Chr_19_peak_genes$gene_name <- Chr_19_peak_3_BM$external_gene_name[match(Chr_19_peak_genes$gene_id, Chr_19_peak_3_BM$ensembl_gene_id)]
Chr_19_peak_gtf$gene_name <- Chr_19_peak_3_BM$external_gene_name[match(Chr_19_peak_gtf$gene_id, Chr_19_peak_3_BM$ensembl_gene_id)]

write.table(Chr_19_peak_genes, "/pathtoexport/Chr19_Genomic_Ranges_Annotations.txt")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$
#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 5)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(0, 60), main = "")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 35, y1 = 35, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      #if(genes_gtf$core[gene_idx]){
      #  text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      #}
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

# Usage:
DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR19,
                        #site_filter = site_filter,
                        peak_GR = Chr_19_peak_GR, 
                        full_gtf = Chr_19_peak_gtf, 
                        plot_reg = c(min(filtered_data_CHR19$POS) - 1e5, 
                                     max(filtered_data_CHR19$POS) + 1e5), 
                        pdf_file = "Chr_19_SNPChip_plot.pdf")

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

#Zoom in - GENE TRACKS

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(BONF)", ylim = c(29, 70), xlim = c(Start_Sign, End_Sign),  main = "",  yaxt = "n")
  
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(30, 70), xlim = c(Start_Sign, End_Sign), main = "", yaxt = "n")
  #axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  
  # Plot SNPs
  #points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
  #       col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  #points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 25, y1 = 25, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 47.5 + y_offset, ytop = 49.5 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 48.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      if(genes_gtf$core[gene_idx]){
        text(x = mid(genes_gtf[target_gene_entry]), srt=45, cex=0.7, y = c(40, 61)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
      }
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR19,
                        #site_filter = site_filter,
                        peak_GR = Chr_19_peak_GR, 
                        full_gtf = Chr_19_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_19_SNPChip_zoomed_plot_GENE_TRACK.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 3, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "", ylim = c(0, 40), xlim = c(Start_Sign, End_Sign), main = "", xaxt = "n")
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  #segments(x0 = start_pos, x1 = end_pos, y0 = 42, y1 = 42, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR19,
                        #site_filter = site_filter,
                        peak_GR = Chr_19_peak_GR, 
                        full_gtf = Chr_19_peak_gtf, 
                        plot_reg = c(start_pos - 1e5, 
                                     end_pos + 1e5), 
                        pdf_file = "Chr_19_SNPChip_zoomed_plot_CONTRAST.pdf")


#Zoom in - CONTRAST

DAF_and_gene_plot_baseR <- function(freq_df, peak_GR, full_gtf = NULL, plot_reg, pdf_file = NULL) {
  
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 5, width = 10)
  
  # Calculate start and end positions for the peak segment
  peak <- 9023505
  start_pos <- Start_Sign
  end_pos <- End_Sign
  
  # Set up plot
  #plot(plot_reg, c(0, 1), type = "n", xlab = "Position", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 43), main = "")
  plot(plot_reg, c(0, 1), type = "n", xlab = "", ylab = "-log10(Bonferroni corrected p)", ylim = c(0, 40), , main = "", xaxt = "n")
  axis(3, tcl = -0.5)  # Placing x-axis ticks at the top
  mtext("Position (Mb)", side = 3, line = 2)
  
  # Plot SNPs
  points(freq_df$POS, -log10(freq_df$BONF), pch = 20, 
         col = ifelse(abs(log10(freq_df$BONF)) > 1.30102999566, "firebrick", "black"))
  
  # Plot missense SNPs
  points(freq_df$POS[freq_df$missense_pos], -log10(freq_df$BONF[freq_df$missense_pos]), pch = 20, cex = 1, col = "lightgreen")
  
  # Plot peak segment
  segments(x0 = start_pos, x1 = end_pos, y0 = 21, y1 = 21, col = "darkorchid", lwd = 2)
  
  # Plot gene features
  #if(!is.null(full_gtf)){
  #  gene_idx <-1
  #  for(target_gene in genes_gtf$gene_id){
  #    y_offset <-  c(0, 3)[1 + (gene_idx %% 2)]
  #    target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
  #    if(sum(target_cds) > 0)rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
  #    target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
  #    if(sum(target_utr) > 0) rect(ybottom = 47 + y_offset, ytop = 48 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
  #    target_gene_entry <- genes_gtf$gene_id == target_gene
  #    segments(y0 = 47.5 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
  #    if(genes_gtf$core[gene_idx]){
  #      text(x = mid(genes_gtf[target_gene_entry]), srt=60, cex=0.5, y = c(42, 57)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry]) 
  #    }
  #    gene_idx <- gene_idx + 1
  #  }
  #}
  if(pdf_file != "") dev.off()
}

#$#$#$#$#$#$#$#$#$#$#$#$#$##$

DAF_and_gene_plot_baseR(freq_df = filtered_data_CHR19,
                        #site_filter = site_filter,
                        peak_GR = Chr_19_peak_GR, 
                        full_gtf = Chr_19_peak_gtf, 
                        plot_reg = c(start_pos, 
                                     max(filtered_data_CHR19$POS) + 1e3), 
                        pdf_file = "Chr_19_SNPChip_zoomed_plot_CONTRAST_WIDE.pdf")


