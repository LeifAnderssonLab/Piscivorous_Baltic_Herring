#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se



#Slåttersill etc frequency contrasts
#Han et al pool frequencies
load("~/Projects/Herring/data/v2.0.2_genotypes/pool_freq_60Neff.RData")


require(tidyverse)

all_large_frq <- read.table(header = F, sep = "\t", skip = 1, col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_1", "FREQ_2"), file = "~/Projects/Herring/data/slattersill_etc/All_large_herring.frq")
all_large_frq$join_check_SNP_ID <- paste(all_large_frq$CHR, all_large_frq$POS, sep = "_")
pool_freq_SNP_ID <- paste(pool_freq$CHROM, pool_freq$POS, sep = "_")

merged_freq <- inner_join(x = pool_freq, y = all_large_frq, by = c("CHROM", "POS"))
save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_v2.RData")

Baltic_Spr_pops <- grep("_Baltic_Spring", names(merged_freq), value = T)
Baltic_Spr_pops <- Baltic_Spr_pops[-grep("LandvikS17", Baltic_Spr_pops)]

merged_freq$Baltic_Spr <- rowMeans(merged_freq[,Baltic_Spr_pops], na.rm = T)
#save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_v1.RData")

Baltic_Aut_pops <- grep("_Baltic_Autumn", names(merged_freq), value = T)
Baltic_Aut_pops <- Baltic_Aut_pops[-grep("HGS5_Schlei", Baltic_Aut_pops)]

merged_freq$Baltic_Aut <- rowMeans(merged_freq[,Baltic_Aut_pops], na.rm = T)

British_pops <- c("HGS16_Orkney_NorthSea_Autumn", "HGS17_IsleOfMan_IrishSea_Autumn", "HGS18_CelticSea_Atlantic_AutumnWinter", "HGS19_TeelinBay_Atlantic_Winter", "HGS20_CapeWrath_Atlantic_Spring", "HGS21_Hebrides_Atlantic_Mixed", "HGS22_CapeWrath_Atlantic_Autumn", "HGS23_Clyde_Atlantic_Spring")
#Baltic_Aut_pops <- Baltic_Aut_pops[-grep("HGS5_Schlei", Baltic_Aut_pops)]

merged_freq$British <- rowMeans(merged_freq[,British_pops], na.rm = T)

merged_freq$Baltic_Spr_v_All_large <- abs(merged_freq$Baltic_Spr - merged_freq$FREQ_1)
names(merged_freq)[grep("FREQ_1", names(merged_freq))] <- "All_large"

#Adding remaining groups
freq_files <- dir("~/Projects/Herring/data/slattersill_etc", full.names = T)
freq_files <- freq_files[grep("frq.gz", freq_files)]
freq_files <- freq_files[-grep("All_large", freq_files)]

#Adding the spawn-switchers
freq_files_switch <- dir("~/Projects/Herring/data/slattersill_etc", full.names = T)
freq_files_switch <- freq_files_switch[grep("frq.gz", freq_files_switch)]
freq_files_switch <- freq_files_switch[grep("GA_SpringSpawner|GS_AutumnSpawner|AutumnSpawner_regular", freq_files_switch)]
#Regular spring-spawners
freq_files_cleanup <- "~/Projects/Herring/data/slattersill_etc/SpringSpawner_regular.frq.gz"

#for(freq_file in freq_files){
for(freq_file in freq_files_switch){
  tmp_large_frq <- read.table(header = F, sep = "\t", skip = 1, col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_1", "FREQ_2"), file = freq_file)
  #tmp_large_frq <- read.table(header = F, sep = "\t", skip = 1, col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_1", "FREQ_2"), file = freq_files_cleanup)
  tmp_large_frq$join_check_SNP_ID <- paste(tmp_large_frq$CHR, tmp_large_frq$POS, sep = "_")
  tmp_group_name <- sub(".+[/]([A-Za-z]+)[.]frq[.]gz", "\\1", freq_file)
  names(tmp_large_frq)[grep("FREQ_1", names(tmp_large_frq))] <- tmp_group_name
  tmp_large_frq <- tmp_large_frq[tmp_large_frq$N_CHR >= 0.8*max(tmp_large_frq$N_CHR),]
  merged_freq <- left_join(x = merged_freq, y = tmp_large_frq[, c("CHROM", "POS", tmp_group_name)], by = c("CHROM", "POS"))
  rm(tmp_large_frq)
}
#save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_v3.RData")
#save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_v4.RData")
save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_v5.RData")
names(merged_freq)[grep("Gavle_large", names(merged_freq))] <- "Gavle_large"
names(merged_freq)[grep("Gavle_regular", names(merged_freq))] <- "Gavle_regular"

names(merged_freq)[grep("GA_SpringSpawner", names(merged_freq))] <- "GA_SpringSpawner"
names(merged_freq)[grep("GS_AutumnSpawner", names(merged_freq))] <- "GS_AutumnSpawner"
names(merged_freq)[grep("AutumnSpawner_regular", names(merged_freq))] <- "AutumnSpawner_regular"


merged_freq$Baltic_Spr_v_GS_AutumnSpawner <- abs(merged_freq$Baltic_Spr - merged_freq$GS_AutumnSpawner)
merged_freq$Baltic_Aut_v_GS_AutumnSpawner <- abs(merged_freq$Baltic_Aut - merged_freq$GS_AutumnSpawner)

merged_freq$Baltic_Spr_v_GA_SpringSpawner <- abs(merged_freq$Baltic_Spr - merged_freq$GA_SpringSpawner)
merged_freq$Baltic_Aut_v_GA_SpringSpawner <- abs(merged_freq$Baltic_Aut - merged_freq$GA_SpringSpawner)

merged_freq$Baltic_Spr_v_AutumnSpawner_regular <- abs(merged_freq$Baltic_Spr - merged_freq$AutumnSpawner_regular)
merged_freq$Baltic_Aut_v_AutumnSpawner_regular <- abs(merged_freq$Baltic_Aut - merged_freq$AutumnSpawner_regular)



merged_freq$Baltic_Spr_v_Slatter <- abs(merged_freq$Baltic_Spr - merged_freq$Slattersill)

slatter_like_pops <- grep("Gavle_large|Oster|Slattersill", names(merged_freq), value = T)

merged_freq$slatter_like <- rowMeans(merged_freq[,slatter_like_pops], na.rm = T)
merged_freq$Baltic_Spr_v_slatter_like <- abs(merged_freq$Baltic_Spr - merged_freq$slatter_like)

slatterish_pops <- grep("Gavle_large|Oster|Slattersill|Sthlm|Kalmar$", names(merged_freq), value = T)
merged_freq$slatterish <- rowMeans(merged_freq[,slatterish_pops], na.rm = T)
merged_freq$Baltic_Spr_v_slatterish <- abs(merged_freq$Baltic_Spr - merged_freq$slatterish)
merged_freq$Baltic_Aut_v_slatterish <- abs(merged_freq$Baltic_Aut - merged_freq$slatterish)

odd_large_pops <- grep("Finland|Blekinge", names(merged_freq), value = T)
merged_freq$odd_large_pops <- rowMeans(merged_freq[,odd_large_pops], na.rm = T)
merged_freq$Baltic_Spr_v_odd_large <- abs(merged_freq$Baltic_Spr - merged_freq$odd_large_pops)
merged_freq$Baltic_Aut_v_odd_large <- abs(merged_freq$Baltic_Aut - merged_freq$odd_large_pops)

SOK_pops <- grep("Sthlm|Oster|Kalmar$", names(merged_freq), value = T)
merged_freq$SOK <- rowMeans(merged_freq[,SOK_pops], na.rm = T)
merged_freq$Baltic_Spr_v_SOK <- abs(merged_freq$Baltic_Spr - merged_freq$SOK)
merged_freq$Baltic_Aut_v_SOK <- abs(merged_freq$Baltic_Aut - merged_freq$SOK)
merged_freq$British_v_SOK <- abs(merged_freq$British - merged_freq$SOK)



#Adding annotation information
merged_freq$NonSyn <- Baltic_EastAtl_ann_df$NonSyn[match(merged_freq$join_check_SNP_ID, Baltic_EastAtl_ann_df$SNP_id)]
merged_freq$missense_pos <- 1:length(merged_freq$British_v_SOK) %in% which(merged_freq$NonSyn) #NA free version

save(merged_freq, file = "~/Projects/Herring/data/slattersill_etc/merged_freq_NonSyn.RData")
merged_freq$Spr_v_Aut_Balt <- abs(merged_freq$Baltic_Aut - merged_freq$Baltic_Spr)
plot(y = merged_freq$Spr_v_Aut_Balt[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)

Atlantic_Aut_pops <- grep("_Atlantic_(Autumn|Winter)", names(merged_freq), value = T)
merged_freq$Atlantic_Aut <- rowMeans(merged_freq[,Atlantic_Aut_pops], na.rm = T)

Atlantic_Spr_pops <- grep("Atlantic_Spring|LandvikS17|RingkobingFjord", names(merged_freq), value = T)
merged_freq$Atlantic_Spr <- rowMeans(merged_freq[,Atlantic_Spr_pops], na.rm = T)
merged_freq$Spr_v_Aut_Atl <- abs(merged_freq$Atlantic_Aut - merged_freq$Atlantic_Spr)

merged_freq$Atl_Spr_v_All_large <- abs(merged_freq$All_large - merged_freq$Atlantic_Spr)
merged_freq$Balt_Aut_v_All_large <- abs(merged_freq$All_large - merged_freq$Baltic_Aut)



merged_freq_GR <- GRanges(seqnames = merged_freq$CHROM, IRanges(start = merged_freq$POS, end = merged_freq$POS))
merged_freq_GR$Spr_v_Aut_Balt <- merged_freq$Spr_v_Aut_Balt
merged_freq_GR$Spr_v_Aut_Atl <- merged_freq$Spr_v_Aut_Atl
merged_freq_GR$Atl_v_Balt_Spr <- abs(merged_freq$Atlantic_Spr - merged_freq$Baltic_Spr)
merged_freq_GR$Atl_v_Balt_Aut <- abs(merged_freq$Atlantic_Aut - merged_freq$Baltic_Aut)


#Cumulative position
merged_freq$cumulative_pos <- NA
merged_freq$cumulative_pos <- merged_freq$POS + Ch_v2.0.2_sizes$offset[match(merged_freq$CHROM, Ch_v2.0.2_sizes$name)] 
merged_freq$manhattan_col <- Ch_v2.0.2_sizes$col[match(merged_freq$CHROM, Ch_v2.0.2_sizes$name)] 

manhattan_site_filter <- merged_freq$N_CHR >= 152 & grepl("chr",  merged_freq$CHROM)

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Baltic_Spr_vs_all_large.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Baltic_Spr_v_All_large[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Baltic_Spr_v_All_large[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Baltic_Aut_vs_all_large.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Balt_Aut_v_All_large[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Balt_Aut_v_All_large[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Atlantic_Spr_vs_all_large.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Atl_Spr_v_All_large[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Atl_Spr_v_All_large[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Baltic_Spr_vs_Slatter_like.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Baltic_Spr_v_slatter_like[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Baltic_Spr_v_slatter_like[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Baltic_Spr_vs_SOK.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Baltic_Spr_v_SOK[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Baltic_Spr_v_SOK[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()

png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/WG_Baltic_Spr_vs_odd_large.png"), width = 2000, height = 500 )
tmp_manhattan_roll_DAF <- stats::filter(merged_freq$Baltic_Spr_v_odd_large[manhattan_site_filter], rep(1/50, 50))
plot(y = merged_freq$Baltic_Spr_v_odd_large[manhattan_site_filter], 
     x = merged_freq$cumulative_pos[manhattan_site_filter],
     col = merged_freq$manhattan_col[manhattan_site_filter], pch = 20, cex = 0.3)
lines(y = tmp_manhattan_roll_DAF,
      x = merged_freq$cumulative_pos[manhattan_site_filter],
      col = "darkorchid")
rm(tmp_manhattan_roll_DAF)
dev.off()


plot(y = merged_freq_GR$Spr_v_Aut_Balt[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], x = merged_freq_GR@ranges@start[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
plot(y = merged_freq_GR$Spr_v_Aut_Atl[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], x = merged_freq_GR@ranges@start[merged_freq$CHR == "chr15" & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)



for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_all_large.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_All_large[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}

for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_slatter.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_Slatter[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}

for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_slatterlike.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_slatter_like[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}

for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_slatterish.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_slatterish[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}

for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Baltic_Aut/", chr, "_Baltic_Aut_vs_odd_large.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Aut_v_odd_large[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_odd_large.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_odd_large[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}

for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Baltic_Aut/", chr, "_Baltic_Aut_vs_SOK.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Aut_v_SOK[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", chr, "_Baltic_Spr_vs_SOK.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_SOK[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}
png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/", "chr12", "_British_vs_SOK.png"), width = 1000, height = 500 )
plot(y = merged_freq$British_v_SOK[merged_freq$CHR == "chr12" & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == "chr12" & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
dev.off()

#Switchers
for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Switchers/", chr, "_Baltic_Spr_vs_GS_aut.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_GS_AutumnSpawner[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Switchers/", chr, "_Baltic_Aut_vs_GS_aut.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Aut_v_GS_AutumnSpawner[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}


for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Switchers/", chr, "_Baltic_Spr_vs_GA_spr.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_GA_SpringSpawner[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Switchers/", chr, "_Baltic_Aut_vs_GA_spr.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Aut_v_GA_SpringSpawner[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}


#Controls
for(chr in paste0("chr", 1:26)){
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Controls/", chr, "_Baltic_Spr_vs_regular_aut.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Spr_v_AutumnSpawner_regular[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
  png(filename = paste0("~/Projects/Herring/doc/Large_herring_WGS/Controls/", chr, "_Baltic_Aut_vs_regular_aut.png"), width = 1000, height = 500 )
  plot(y = merged_freq$Baltic_Aut_v_AutumnSpawner_regular[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], x = merged_freq$POS[merged_freq$CHR == chr & merged_freq$N_CHR >= 140], pch = 20, cex = 0.5)
  dev.off()
}


#Heatmaps
#hm_col_filter <- grep("join|N_ALL|N_CHR|FREQ_2|_v_", names(merged_freq))
hm_col_filter_ext <- grep("join|N_ALL|N_CHR|FREQ_2|_v_|slatter_like|slatterish|odd_large_pops|SOK|GA_SpringSpawner|GS_AutumnSpawner|SpringSpawner_regular|All_large|NonSyn|missense_pos|cumulative_pos|manhattan_col", names(merged_freq))
hm_col_filter <- unique(c(hm_col_filter_ext, 3:62)) #more restriced version, emphasising the new pools
#pool_order_vec_large_ext <- (1:length(names(merged_freq)[-hm_col_filter_ext]))[-grep("POS|CHROM|NonSyn|missense_pos", names(merged_freq)[-hm_col_filter_ext])]
ext_pool_names <- names(merged_freq)[-c(1,2,hm_col_filter_ext)]
#pool_order_vec_large_ext <- (1:length(ext_pool_names))[-grep("POS|CHROM|NonSyn|missense_pos", ext_pool_names)]
pool_order_vec_large <- c(13,3,7,12,14,4,5,8,11,9,6,10)

balt_spr_block <- grep("Baltic_(Spring|Summer)|TysklandS18", ext_pool_names)
balt_aut_block <- grep("Baltic_Autumn", ext_pool_names)
atl_spr_block <- grep("Atlantic_Spring|LandvikS17|RingkobingFjord", ext_pool_names)
atl_aut_block <- grep("Atlantic_Autumn", ext_pool_names)
brit_block <- grep("Hebrides|Orkney|CapeWrath|IsleOfMan|Downs|Clyde|TeelinBay", ext_pool_names)
pac_block <- grep("HWS|Pacific", ext_pool_names)
large_block <- grep("Gavle_large|Slattersill|Sthlm|Ostergotland|Kalmar$|Finland|Blekinge", ext_pool_names)
new_block <- grep("Gavle_regular|AutumnSpawner_regular", ext_pool_names)
ref_block <- grep("Baltic_Spr$|Baltic_Aut$|British", ext_pool_names)
balt_spr_block <- balt_spr_block[!balt_spr_block %in% atl_spr_block]
atl_spr_block <- atl_spr_block[!atl_spr_block %in% pac_block]
atl_spr_block <- atl_spr_block[!atl_spr_block %in% brit_block]
atl_aut_block <- atl_aut_block[!atl_aut_block %in% brit_block]

pool_order_vec_large_ext <- c(pac_block, balt_spr_block, balt_aut_block, atl_spr_block, atl_aut_block, brit_block, new_block, large_block[c(1,2,3,6,7,5,4)], ref_block)


chr5_plot_filter <- which(merged_freq$CHROM == "chr5" & merged_freq$POS > 3.5e6 & merged_freq$POS < 4.5e6 & merged_freq$N_CHR >= 140 & (merged_freq$Baltic_Spr_v_slatterish >= 0.35 | merged_freq$Baltic_Aut_v_slatterish >= 0.35))
tmp_freq <- merged_freq[chr5_plot_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
#pool_freq_hm_v3(tmp_freq, pool_order_vec = 3:dim(tmp_freq)[2], margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_5_BS_v_Large_hm.pdf", col_lab_vec = rownames(tmp_freq))
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_5_4Mb_hm.pdf", col_lab_vec = rownames(tmp_freq))


chr8_plot_filter <- which(merged_freq$CHROM == "chr8" & merged_freq$POS > 11e6 & merged_freq$POS < 13e6 & merged_freq$N_CHR >= 140 & (merged_freq$Baltic_Spr_v_slatterish >= 0.35 | merged_freq$Baltic_Aut_v_slatterish >= 0.35))
tmp_freq <- merged_freq[chr8_plot_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
#pool_freq_hm_v3(tmp_freq, pool_order_vec = 3:dim(tmp_freq)[2], margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_8_BS_v_Large_hm.pdf", col_lab_vec = rownames(tmp_freq))
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_8_12Mb_hm.pdf", col_lab_vec = rownames(tmp_freq))


chr12_plot_filter <- which(merged_freq$CHROM == "chr12" & merged_freq$POS > 17e6 & merged_freq$POS < 27e6 & merged_freq$N_CHR >= 140 & (merged_freq$Baltic_Spr_v_slatterish >= 0.35 | merged_freq$Baltic_Aut_v_slatterish >= 0.35))
tmp_freq <- merged_freq[chr12_plot_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
#pool_freq_hm_v3(tmp_freq, pool_order_vec = 3:dim(tmp_freq)[2], margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_BS_v_Large_hm.pdf", col_lab_vec = rownames(tmp_freq))
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_22Mb_hm.pdf", col_lab_vec = rownames(tmp_freq))


chr20_plot_filter <- which(merged_freq$CHROM == "chr20" & merged_freq$POS > 12.5e6 & merged_freq$POS < 17e6 & merged_freq$N_CHR >= 140 & merged_freq$Baltic_Spr_v_All_large >= 0.25)
tmp_freq <- merged_freq[chr20_plot_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
#pool_freq_hm_v3(tmp_freq, pool_order_vec = 3:dim(tmp_freq)[2], margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_20_BS_v_Large_hm.pdf", col_lab_vec = rownames(tmp_freq))
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_20_15Mb_hm.pdf", col_lab_vec = rownames(tmp_freq))

chr20_scatter_filter <- which(merged_freq$CHROM == "chr20" & merged_freq$POS > 11.5e6 & merged_freq$POS < 17.5e6 & merged_freq$N_CHR >= 140)
chr20_freq_est_filter <- which(merged_freq$CHROM == "chr20" & merged_freq$POS > 12.5e6 & merged_freq$POS < 17e6 & merged_freq$N_CHR >= 140 & merged_freq$Baltic_Spr_v_All_large >= 0.25 & (merged_freq$British < 0.1|merged_freq$British > 0.9))


tmp_BS_ref_freq <- merged_freq$Baltic_Spr[chr20_freq_est_filter]
loci_to_invert <- tmp_BS_ref_freq > 0.5

pdf("~/Projects/Herring/doc/Large_herring_WGS/Chr_20_15_Mb_DAF.pdf")
plot(x = merged_freq$POS[chr20_scatter_filter], y = merged_freq$British_v_SOK[chr20_scatter_filter], pch = 20, ylim = c(0,1), cex = 0.6)
dev.off()

tmp_plot_freq_BS <- merged_freq$Baltic_Spr[chr20_freq_est_filter]
tmp_plot_freq_BS[loci_to_invert] <- 1 - tmp_plot_freq_BS[loci_to_invert]

pdf("~/Projects/Herring/doc/Large_herring_WGS/Chr_20_15_Mb_freq.pdf")
plot(x = merged_freq$POS[chr20_freq_est_filter], y = tmp_plot_freq_BS, col = "darkorchid", pch = 20, ylim = c(0,1))

tmp_plot_freq_Sthlm <- merged_freq$Sthlm[chr20_freq_est_filter]
tmp_plot_freq_Sthlm[loci_to_invert] <- 1 - tmp_plot_freq_Sthlm[loci_to_invert]
points(x = merged_freq$POS[chr20_freq_est_filter], y = tmp_plot_freq_Sthlm, col = "firebrick1", pch = 20)
#tmp_plot_freq_Ostergotland <- merged_freq$Ostergotland[chr20_freq_est_filter]
#tmp_plot_freq_Ostergotland[loci_to_invert] <- 1 - tmp_plot_freq_Ostergotland[loci_to_invert]
#points(x = merged_freq$POS[chr20_freq_est_filter], y = tmp_plot_freq_Ostergotland, col = "firebrick4", pch = 20)

tmp_plot_freq_Slattersill <- merged_freq$Slattersill[chr20_freq_est_filter]
tmp_plot_freq_Slattersill[loci_to_invert] <- 1 - tmp_plot_freq_Slattersill[loci_to_invert]
points(x = merged_freq$POS[chr20_freq_est_filter], y = tmp_plot_freq_Slattersill, col = "olivedrab", pch = 20)

tmp_plot_freq_British <- merged_freq$British[chr20_freq_est_filter]
tmp_plot_freq_British[loci_to_invert] <- 1 - tmp_plot_freq_British[loci_to_invert]
points(x = merged_freq$POS[chr20_freq_est_filter], y = tmp_plot_freq_British, col = "steelblue", pch = 20)

dev.off()

plot(x = tmp_plot_freq_Sthlm, y = tmp_plot_freq_Ostergotland, pch = 20, ylim = c(0,1), cex = 0.6, xlim = c(0,1))
lm1 <- summary(lm(tmp_plot_freq_Ostergotland~tmp_plot_freq_Sthlm-1))
abline(a= 0, b = lm1$coefficients[1]); text(x = 0.1, y = 0.9, labels = paste0("k = ", round(lm1$coefficients[1], digits = 2), "\nr^2 = ", round(lm1$r.squared, digits = 2)), pos = 4)
abline(a= 0, b = 1, col = "red")

plot(x = tmp_plot_freq_Sthlm, y = tmp_plot_freq_Slattersill, pch = 20, ylim = c(0,1), cex = 0.6, xlim = c(0,1))
lm1 <- summary(lm(tmp_plot_freq_Slattersill~tmp_plot_freq_Sthlm-1))
abline(a= 0, b = lm1$coefficients[1]); text(x = 0.1, y = 0.9, labels = paste0("k = ", round(lm1$coefficients[1], digits = 2), "\nr^2 = ", round(lm1$r.squared, digits = 2)), pos = 4)
abline(a= 0, b = 1, col = "red")

plot(x = tmp_plot_freq_Sthlm, y = tmp_plot_freq_BS, pch = 20, ylim = c(0,1), cex = 0.6, xlim = c(0,1))
lm1 <- summary(lm(tmp_plot_freq_BS~tmp_plot_freq_Sthlm-1))
abline(a= 0, b = lm1$coefficients[1]); text(x = 0.1, y = 0.9, labels = paste0("k = ", round(lm1$coefficients[1], digits = 2), "\nr^2 = ", round(lm1$r.squared, digits = 2)), pos = 4)
abline(a= 0, b = 1, col = "red")

plot(x = tmp_plot_freq_Sthlm, y = tmp_plot_freq_British, pch = 20, ylim = c(0,1), cex = 0.6, xlim = c(0,1))
lm1 <- summary(lm(tmp_plot_freq_British~tmp_plot_freq_Sthlm-1))
abline(a= 0, b = lm1$coefficients[1]); text(x = 0.1, y = 0.9, labels = paste0("k = ", round(lm1$coefficients[1], digits = 2), "\nr^2 = ", round(lm1$r.squared, digits = 2)), pos = 4)
abline(a= 0, b = 1, col = "red")

pdf("~/Projects/Herring/doc/Large_herring_WGS/Chr_20_15_Mb_hist.pdf")

hist(tmp_plot_freq_British, breaks = 10, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), xlim = c(0,1))
hist(tmp_plot_freq_Sthlm, breaks = 10, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5), add = T)
hist(tmp_plot_freq_BS, breaks = 10, add = T, col = rgb(red = 1, green = 0, blue = 1, alpha = 0.5))
dev.off()

#points(x = merged_freq$POS[chr20_plot_filter], y = 1-merged_freq$Slattersill[chr20_plot_filter], col = "olivedrab1", pch = 20)
#points(x = merged_freq$POS[chr20_plot_filter], y = 1-merged_freq$Finland[chr20_plot_filter], col = "olivedrab4", pch = 20)





chr19_plot_filter <- which(merged_freq$CHROM == "chr19" & merged_freq$POS > 19e6 & merged_freq$POS < 20.8e6 & merged_freq$N_CHR >= 140 & merged_freq$Baltic_Spr_v_All_large >= 0.45)
tmp_freq <- merged_freq[chr19_plot_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
#pool_freq_hm_v3(tmp_freq, pool_order_vec = 3:dim(tmp_freq)[2], margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_19_BS_v_Large_hm.pdf", col_lab_vec = rownames(tmp_freq))
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_19_20Mb_hm.pdf", col_lab_vec = rownames(tmp_freq))

#THRB (not large herring related)
chr19_THRB_filter <- which(merged_freq$CHROM == "chr19" & merged_freq$POS > 6.3e6 & merged_freq$POS < 6.5e6 & merged_freq$N_CHR >= 140)
tmp_freq <- merged_freq[chr19_THRB_filter,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_19_6.3Mb_THRB_hm.pdf", col_lab_vec = rownames(tmp_freq))

chr19_THRB_filter_NonSyn <- which(merged_freq$CHROM == "chr19" & merged_freq$POS > 6.3e6 & merged_freq$POS < 6.5e6 & merged_freq$N_CHR >= 140 & merged_freq$missense_pos)
tmp_freq <- merged_freq[chr19_THRB_filter_NonSyn,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_19_6.3Mb_THRB_NonSyn_hm.pdf", col_lab_vec = rownames(tmp_freq))

#SEC16B (not large herring related)
chr10_SEC16B_filter <- which(merged_freq$CHROM == "chr10" & merged_freq$POS > 25.1e6 & merged_freq$POS < 25.2e6 & merged_freq$N_CHR >= 140)
tmp_freq <- merged_freq[chr10_SEC16B_filter,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr10_25.15_Mb_SEC16B_hm.pdf", col_lab_vec = rownames(tmp_freq))

chr10_SEC16B_filter_NonSyn <- which(merged_freq$CHROM == "chr10" & merged_freq$POS > 25.128e6 & merged_freq$POS < 25.146e6 & merged_freq$N_CHR >= 140 & merged_freq$missense_pos)
tmp_freq <- merged_freq[chr10_SEC16B_filter_NonSyn,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr10_25.15_Mb_SEC16B_NonSyn_hm.pdf", col_lab_vec = rownames(tmp_freq))

#Red-opsin cluster (not large herring related)
chr5_red_opsin_filter <- which(merged_freq$CHROM == "chr5" & merged_freq$POS > 4.15e6 & merged_freq$POS < 4.21e6 & merged_freq$N_CHR >= 140 & abs(merged_freq$Atlantic_Spr-merged_freq$Baltic_Spr) > 0.25)
tmp_freq <- merged_freq[chr5_red_opsin_filter,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr5_4.2_Mb_red_opsin_hm.pdf", col_lab_vec = rownames(tmp_freq))
plot(x = merged_freq$POS[chr5_red_opsin_filter], y = abs(merged_freq$Atlantic_Spr[chr5_red_opsin_filter]-merged_freq$Baltic_Spr[chr5_red_opsin_filter]))
tmp_5_filter <- which(merged_freq$CHROM == "chr5" & merged_freq$N_CHR >= 140)
plot(x = merged_freq$POS[tmp_5_filter], y = abs(merged_freq$Atlantic_Spr[tmp_5_filter]-merged_freq$Baltic_Spr[tmp_5_filter]))
abline(v = 4.186e6)


#Cyt27c1 cluster (not large herring related)
chr21_Cyt27c1_filter <- which(merged_freq$CHROM == "chr21" & merged_freq$POS > 25.7e6 & merged_freq$POS < 25.9e6 & merged_freq$N_CHR >= 140) # & abs(merged_freq$Atlantic_Spr-merged_freq$Baltic_Spr) > 0.25)
tmp_freq <- merged_freq[chr5_red_opsin_filter,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr5_4.2_Mb_red_opsin_hm.pdf", col_lab_vec = rownames(tmp_freq))
plot(x = merged_freq$POS[chr21_Cyt27c1_filter], y = abs(merged_freq$Atlantic_Aut[chr21_Cyt27c1_filter]-merged_freq$Baltic_Aut[chr21_Cyt27c1_filter]))

tmp_21_filter <- which(merged_freq$CHROM == "chr21" & merged_freq$N_CHR >= 140)

#SOK-specipci regions on Chr 12 inversion
SOK_12_site_filter <- merged_freq$CHROM == "chr12" & merged_freq$N_CHR >= 140 & merged_freq$POS > 17.5e6 & merged_freq$POS < 25.7e6
SOK_12_site_GR <- GRanges(seqnames = merged_freq$CHROM[SOK_12_site_filter], ranges = IRanges(start = merged_freq$POS[SOK_12_site_filter], end = merged_freq$POS[SOK_12_site_filter]))
SOK_12_site_GR$British_v_SOK <- merged_freq$British_v_SOK[SOK_12_site_filter]
SOK_12_site_GR$missense_pos <- merged_freq$missense_pos[SOK_12_site_filter]
SOK_12_peak_GR <- reduce(SOK_12_site_GR[which(SOK_12_site_GR$British_v_SOK > 0.55)], min.gapwidth = 5e4)
SOK_12_peak_GR <- SOK_12_peak_GR[width(SOK_12_peak_GR) > 1e2]
SOK_12_peak_SNP_hits <- findOverlaps(SOK_12_peak_GR, SOK_12_site_GR)
SOK_12_peak_SNP_GR <- SOK_12_site_GR[SOK_12_peak_SNP_hits@to]
SOK_12_peak_SNP_GR$peak <- SOK_12_peak_GR[SOK_12_peak_SNP_hits@from]
SOK_12_peak_3_NonSyn <- SOK_12_peak_SNP_GR[SOK_12_peak_SNP_GR$British_v_SOK > 0.6 & SOK_12_peak_SNP_GR$missense_pos]



plot(y = merged_freq$British_v_SOK[SOK_12_site_filter], x = merged_freq$POS[SOK_12_site_filter], pch = 20, cex = 0.5, ylim = c(0,1) )
points(y = merged_freq$British_v_SOK[SOK_12_site_filter & merged_freq$missense_pos], x = merged_freq$POS[SOK_12_site_filter & merged_freq$missense_pos], pch = 20, cex = 0.6, col = "firebrick")
segments(y0 = 1, x0 = start(SOK_12_peak_GR), x1 = end(SOK_12_peak_GR), lwd = 2, col = "darkorchid")
abline(h = 0.55)


SOK_12_peak_gtf_hits <- findOverlaps(SOK_12_peak_GR, cluhar_v2.0.2_gtf, maxgap = 5e4)
SOK_12_peak_core_hits <- findOverlaps(SOK_12_peak_GR, cluhar_v2.0.2_gtf)

SOK_12_peak_gtf <- cluhar_v2.0.2_gtf[SOK_12_peak_gtf_hits@to]
SOK_12_core_gtf <- cluhar_v2.0.2_gtf[SOK_12_peak_core_hits@to]
SOK_12_peak_gtf$peak <- SOK_12_peak_GR[SOK_12_peak_gtf_hits@from]
SOK_12_peak_gtf$core <- F
SOK_12_peak_gtf$core[SOK_12_peak_gtf$gene_id %in% SOK_12_core_gtf$gene_id] <- T
SOK_12_peak_gtf$gene_name <- SOK_12_peak_3_BM$external_gene_name[match(SOK_12_peak_gtf$gene_id, SOK_12_peak_3_BM$ensembl_gene_id)]
SOK_12_peak_genes <- SOK_12_peak_gtf[SOK_12_peak_gtf$type == "gene"]

mart = useEnsembl(biomart = "ensembl", dataset = "charengus_gene_ensembl")
listAttributes(mart = mart)
#SOK_12_peak_3_BM <- getBM(attributes = c("description", "go_id", "name_1006", "external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = SOK_12_peak_genes$gene_id)
SOK_12_peak_3_BM <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = mart, filters = "ensembl_gene_id", values = SOK_12_peak_genes$gene_id)
SOK_12_peak_genes$gene_name <- SOK_12_peak_3_BM$external_gene_name[match(SOK_12_peak_genes$gene_id, SOK_12_peak_3_BM$ensembl_gene_id)]

SOK_12_peak_3_NonSyn_hits <- findOverlaps(SOK_12_peak_3_NonSyn, SOK_12_peak_genes)
SOK_12_peak_genes[SOK_12_peak_3_NonSyn_hits@to]
amacr_aln <- readAAMultipleAlignment(filepath = "~/Projects/Herring/data/slattersill_etc/amacr_gene_tree.aln", format = "clustal")
amacr_AAbin <- as.AAbin(amacr_aln)

amacr_aln_large <- readAAMultipleAlignment(filepath = "~/Projects/Herring/data/slattersill_etc/amacr_gene_tree_large.aln", format = "clustal")
amacr_large_AAbin <- as.AAbin(amacr_aln_large)


 #zgrep -E "chr12\t25104738" Ch_v2.0.2_67_pools_snpEff.vcf.gz

DAF_and_gene_plot(freq_df = merged_freq, site_filter = SOK_12_site_filter, plot_reg = c(start(SOK_12_peak_GR)[3] - 1e5, end(SOK_12_peak_GR)[3] + 1e5), full_gtf = SOK_12_peak_gtf, peak_GR = SOK_12_peak_GR,  pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/chr12_SOK_peak_3_annot.pdf")
amacr_site_filter <- merged_freq$CHROM == "chr12" & merged_freq$N_CHR >= 140 & merged_freq$POS > 25.0e6 & merged_freq$POS < 25.2e6 & merged_freq$British_v_SOK > 0.55
tmp_freq <- merged_freq[amacr_site_filter,-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_amacr_reg_hm.pdf", col_lab_vec = rownames(tmp_freq))
tmp_freq <- merged_freq[amacr_site_filter,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_amacr_ext_hm.pdf", col_lab_vec = rownames(tmp_freq))

tmp_freq <- merged_freq[amacr_site_filter & merged_freq$missense_pos,-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_amacr_NonSyn_hm.pdf", col_lab_vec = rownames(tmp_freq))



DAF_and_gene_plot(freq_df = merged_freq, site_filter = SOK_12_site_filter, plot_reg = c(start(SOK_12_peak_GR)[1] - 2e5, end(SOK_12_peak_GR)[1] + 2e5), full_gtf = SOK_12_peak_gtf, peak_GR = SOK_12_peak_GR, pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/chr12_SOK_peak_1_annot.pdf")
peak1_site_filter <- merged_freq$CHROM == "chr12" & merged_freq$N_CHR >= 140 & merged_freq$POS > start(SOK_12_peak_GR)[1] - 2e5 & merged_freq$POS < end(SOK_12_peak_GR)[1] + 2e5 & merged_freq$British_v_SOK > 0.55
tmp_freq <- merged_freq[which(peak1_site_filter),-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_SOK_peak_1_hm.pdf", col_lab_vec = rownames(tmp_freq))
tmp_freq <- merged_freq[which(peak1_site_filter),-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_SOK_peak_1_ext_hm.pdf", col_lab_vec = rownames(tmp_freq))


DAF_and_gene_plot(freq_df = merged_freq, site_filter = SOK_12_site_filter, plot_reg = c(start(SOK_12_peak_GR)[2] - 2e5, end(SOK_12_peak_GR)[2] + 2e5), full_gtf = SOK_12_peak_gtf, peak_GR = SOK_12_peak_GR, pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/chr12_SOK_peak_2_annot.pdf")
peak2_site_filter <- merged_freq$CHROM == "chr12" & merged_freq$N_CHR >= 140 & merged_freq$POS > start(SOK_12_peak_GR)[2] - 2e5 & merged_freq$POS < end(SOK_12_peak_GR)[2] + 2e5 & merged_freq$British_v_SOK > 0.55
tmp_freq <- merged_freq[which(peak2_site_filter),-hm_col_filter]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_SOK_peak_2_hm.pdf", col_lab_vec = rownames(tmp_freq))
tmp_freq <- merged_freq[which(peak2_site_filter),-hm_col_filter_ext]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large_ext, margins = c(12, 15), pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/Heatmaps/Chr_12_SOK_peak_2_ext_hm.pdf", col_lab_vec = rownames(tmp_freq))

DAF_and_gene_plot(freq_df = merged_freq, site_filter = SOK_12_site_filter, plot_reg = range(merged_freq$POS[SOK_12_site_filter]), peak_GR = SOK_12_peak_GR, pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/chr12_British_vs_SOK_annot.pdf")


#Collecting preliminary figures (regions from Jake)
reg1_chr10_plot_filter <- which(merged_freq$CHROM == "chr10" & merged_freq$POS > 21428922 & merged_freq$POS < 23428922 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.25)
tmp_freq <- merged_freq[reg1_chr10_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_10_22Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

reg2_chr12_plot_filter <- which(merged_freq$CHROM == "chr12" & merged_freq$POS > 14794449 & merged_freq$POS < 16794449 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.20)
tmp_freq <- merged_freq[reg2_chr12_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_12_16Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))


reg3_chr15_plot_filter <- which(merged_freq$CHROM == "chr15" & merged_freq$POS > 8027574 & merged_freq$POS < 10027574 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.30)
tmp_freq <- merged_freq[reg3_chr15_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_15_9Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

reg4_chr17_plot_filter <- which(merged_freq$CHROM == "chr17" & merged_freq$POS > 25581482 & merged_freq$POS < 27581482 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.25)
tmp_freq <- merged_freq[reg4_chr17_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_17_26Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

reg5_chr18_plot_filter <- which(merged_freq$CHROM == "chr18" & merged_freq$POS > 23353963 & merged_freq$POS < 25353963 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.30)
tmp_freq <- merged_freq[reg5_chr18_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_18_24Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

reg6_chr19_plot_filter <- which(merged_freq$CHROM == "chr19" & merged_freq$POS > 19543808 & merged_freq$POS < 21543808 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.30)
tmp_freq <- merged_freq[reg6_chr19_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_19_20Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

reg7_chr24_plot_filter <- which(merged_freq$CHROM == "chr24" & merged_freq$POS > 16639495 & merged_freq$POS < 18639495 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.20)
tmp_freq <- merged_freq[reg7_chr24_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_24_17Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

regX_chr20_plot_filter <- which( merged_freq$CHROM == "chr20" & merged_freq$POS > 12.5e6 & merged_freq$POS < 17e6 & merged_freq$N_CHR >= 150 & merged_freq$Baltic_Spr_v_slatter_like >= 0.20)
tmp_freq <- merged_freq[regX_chr20_plot_filter,-hm_col_filter]
plot_flip_vec <- tmp_freq$Baltic_Spr < 0.5
tmp_freq[plot_flip_vec,pool_order_vec_large] <- 1 - tmp_freq[plot_flip_vec,pool_order_vec_large]
rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
pool_freq_hm_v3(tmp_freq, pool_order_vec = pool_order_vec_large, margins = c(12, 15),
                pdf_file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/Chr_20_15Mb_hm.pdf", 
                col_lab_vec = rownames(tmp_freq))

#Frequency-based PCA
merged_freq_PCA <- merged_freq[sample.int(n = dim(merged_freq)[1], size = 1e6),ext_pool_names] #Small example
merged_freq_PCA <- merged_freq_PCA[,-grep("HWS|Pacific", names(merged_freq_PCA))]
merged_freq_PCA <- merged_freq_PCA[rowVars(as.matrix(merged_freq_PCA), na.rm = T) > 0 & rowSums(is.na(merged_freq_PCA)) < 1,]

mf_PCA_t <-t(merged_freq_PCA)
mf.pca <- prcomp(mf_PCA_t, scale. = T)
save(mf.pca, merged_freq_PCA, file ="~/Projects/Herring/data/slattersill_etc/WGS_PCA.RData")

#Version using selected markers from Han et al
selected_han_SNPs <- read.table("~/Projects/Herring/data/British_populations/gensinc_revision/pc_analysis/selectiveSNPs.pos", stringsAsFactors = F)
names(selected_han_SNPs) <- c("CHROM", "POS")
selected_han_SNPs$SNP_ID <- paste(selected_han_SNPs$CHROM, selected_han_SNPs$POS, sep = "_")
merged_freq_PCA_selected <- merged_freq[merged_freq$join_check_SNP_ID %in% selected_han_SNPs$SNP_ID,ext_pool_names]
merged_freq_PCA_selected <- merged_freq_PCA_selected[,-grep("HWS|Pacific", names(merged_freq_PCA_selected))]
merged_freq_PCA_selected <- merged_freq_PCA_selected[rowVars(as.matrix(merged_freq_PCA_selected), na.rm = T) > 0 & rowSums(is.na(merged_freq_PCA_selected)) < 1,]
mf_sel_PCA_t <-t(merged_freq_PCA_selected)
mf_sel.pca <- prcomp(mf_sel_PCA_t, scale. = T)


#mf.pca_df <- as.data.frame(mf.pca$x[, 1:2]) # Original, random SNPs
mf.pca_df <- as.data.frame(mf_sel.pca$x[, 1:2]) # Slected SNPs
mf.pca_df <- mf.pca_df[-grep("Baltic_Spr$|Baltic_Aut$|Atlantic_Spr$|Atlantic_Aut$|British$", rownames(mf.pca_df)),] #Removing averages


#mf.pca_df$plot_col = "grey70"
#mf.pca_df$plot_col[grep("Atlantic_Spr|HGS11", rownames(mf.pca_df))] <- "olivedrab1"
#mf.pca_df$plot_col[grep("Atlantic_Aut|HGS10|HGS16|HGS17|HGS19|HGS21|British", rownames(mf.pca_df))] <- "orange1"
#mf.pca_df$plot_col[grep("Baltic_Spr|TysklandS18|Gavle_regular", rownames(mf.pca_df))] <- "olivedrab4"
#mf.pca_df$plot_col[grep("Baltic_Summer", rownames(mf.pca_df))] <- "olivedrab4"
#mf.pca_df$plot_col[grep("Baltic_Aut|AutumnSpawner_regular", rownames(mf.pca_df))] <- "orange4"
#mf.pca_df$plot_col[grep("LandvikS17", rownames(mf.pca_df))] <- "olivedrab1"

#mf.pca_df$plot_col[grep("HWS|Pacific", rownames(mf.pca_df))] <- "black"
#mf.pca_df$plot_col[grep("Gavle_large|Kalmar$|Ostergotland|Slattersill|Sthlm|Blekinge|Finland", rownames(mf.pca_df))] <- "firebrick"
#mf.pca_df$label <- ""
#mf.pca_df$label[mf.pca_df$plot_col == "firebrick"] <- rownames(mf.pca_df)[mf.pca_df$plot_col == "firebrick"]
#mf.pca_df$plot_pos <- 3
#mf.pca_df$plot_pos[rownames(mf.pca_df) %in% c("Sthlm", "Kalmar")] <- 1

#Jake's colours
mf.pca_df$plot_col = "grey70"
mf.pca_df$plot_col[grep("Atlantic_Spr|HGS11", rownames(mf.pca_df))] <- "indianred3"
mf.pca_df$plot_col[grep("Atlantic_Aut|HGS10|HGS16|HGS17|HGS19|HGS21|British", rownames(mf.pca_df))] <- "black"
mf.pca_df$plot_col[grep("Baltic_Spr|TysklandS18|Gavle_regular", rownames(mf.pca_df))] <- "indianred3"
mf.pca_df$plot_col[grep("Baltic_Summer", rownames(mf.pca_df))] <- "indianred3"
mf.pca_df$plot_col[grep("Baltic_Aut|AutumnSpawner_regular", rownames(mf.pca_df))] <- "black"
mf.pca_df$plot_col[grep("LandvikS17", rownames(mf.pca_df))] <- "indianred3"
mf.pca_df$plot_col[grep("Gavle_large|Kalmar$|Ostergotland|Slattersill|Sthlm|Blekinge|Finland", rownames(mf.pca_df))] <- "limegreen"



#pdf(file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/WGS_PCA.pdf", height = 10, width = 10) #Original verision, random SNPs
pdf(file = "~/Projects/Herring/doc/Large_herring_WGS/draft/figure_sketches/WGS_PCA_selected.pdf", height = 10, width = 10) #Selected SNPs
plot(x = mf.pca_df$PC1, y = mf.pca_df$PC2, pch = 20, col = mf.pca_df$plot_col, cex = 2, xlab = "PC1", ylab = "PC2")
#text(x = mf.pca_df$PC1, y = mf.pca_df$PC2, labels = mf.pca_df$label, cex = 1, pos = mf.pca_df$plot_pos)
#legend(x = "topleft", 
#       legend = c("Large Baltic", "Baltic Spring", "Baltic Autumn", "Atlantic Spring", "Atlantic Autumn"),
#       col = c("firebrick", "olivedrab4", "orange4", "olivedrab1", "orange1"), pch = 20, cex = 1.2)
dev.off()



#Support functions
pool_freq_hm_v3 <- function(freq_df, pdf_file = NULL, pool_order_vec, hm_col = c("blue4", "gold1"), margins = c(12,12), reg_name = NULL, row_lab_vec = NULL, col_lab_vec = NULL, ...){
  if(!is.null(pdf_file)) {
    pdf(file = pdf_file, height = 10, width = 10)
  }
  par(xpd = NA)
  if(is.null(row_lab_vec)){
    row_lab_vec <- sub("[A-Z0-9]+_", "", colnames(freq_df)[pool_order_vec])
  }
  
  target_snp_id <- rownames(freq_df)
  if(is.null(col_lab_vec)){
    col_lab_vec <- rep("", length(target_snp_id))
  }
  
  cr <- colorRamp(hm_col)
  crp = colorRampPalette(colors = hm_col)(500)
  col_colors <- rep("black", dim(freq_df)[1])
  col_colors[freq_df$NonSyn] <- "darkorchid"
  heatmap(t(as.matrix(freq_df[,pool_order_vec])), scale = "none", Colv = NA, Rowv = NA, labCol = col_lab_vec, labRow = row_lab_vec, col = crp, margins = margins, ColSideColors = col_colors, ...) #ColSideColors = col_colors
  scale_x_vec <- seq(from = par("usr")[1], to = par("usr")[2]*0.5, length.out = 505)
  rect(ybottom = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.00, xleft = scale_x_vec[-(501:505)], ytop = par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.02, xright =  scale_x_vec[-(1:5)], col = rgb(cr((1:500)/500), maxColorValue=255), border = NA)
  text(x=scale_x_vec[c(10,253,505)], y = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.04, labels = paste0(c(0, 50, 100), "%"))
  if(!is.null(reg_name)) {
    text(x=scale_x_vec[253], y = par("usr")[3] + (par("usr")[4]-par("usr")[3])*0.05, labels = reg_name, cex = 1.3)
  }
  if(!is.null(pdf_file)) {
    dev.off()
  }
  #return(par("usr"))
}

DAF_and_gene_plot <- function(freq_df, site_filter, peak_GR, full_gtf = NULL, plot_reg, pdf_file  = "", DAF_col = "British_v_SOK", gene_cex = 1){
  genes_gtf <- full_gtf[full_gtf$type == "gene"]
  if(pdf_file != "") pdf(file = pdf_file, height = 7, width = 10)
  DAF_vec <- freq_df[,DAF_col]
  plot(y = DAF_vec[site_filter], x = freq_df$POS[site_filter], pch = 20, cex = 0.5, ylim = c(0,1), xlim = plot_reg, ylab = "Delta Allele Frequency", xlab = "Position")
  points(y = DAF_vec[site_filter & freq_df$missense_pos], x = freq_df$POS[site_filter & freq_df$missense_pos], pch = 20, cex = 1, col = "firebrick")
  segments(y0 = 1, x0 = start(peak_GR), x1 = end(peak_GR), lwd = 2, col = "darkorchid")
  abline(h = 0.55)
  if(!is.null(full_gtf)){
    gene_idx <-1
    for(target_gene in genes_gtf$gene_id){
      y_offset <-  c(0,0.05)[1 + (gene_idx %% 2)]
      target_cds <- full_gtf$type == "CDS" & full_gtf$gene_id == target_gene
      if(sum(target_cds) > 0)rect(ybottom = 0.86 + y_offset, ytop = 0.89 + y_offset, xleft = start(full_gtf[target_cds]), xright = end(full_gtf[target_cds]), col = "steelblue", border = NA)
      target_utr <- 1:length(full_gtf$type) %in% grep("utr",full_gtf$type) & full_gtf$gene_id == target_gene
      if(sum(target_utr) > 0) rect(ybottom = 0.87 + y_offset, ytop = 0.88 + y_offset, xleft = start(full_gtf[target_utr]), xright = end(full_gtf[target_utr]), col = "steelblue4", border = NA)
      target_gene_entry <- genes_gtf$gene_id == target_gene
      segments(y0 = 0.875 + y_offset, x0 = start(genes_gtf[target_gene_entry]), x1 = end(genes_gtf[target_gene_entry]), col = "grey30", lwd = 2)
      if(genes_gtf$core[gene_idx]){
        text(x = mid(genes_gtf[target_gene_entry]), y = c(0.83, 0.97)[1 + (gene_idx %% 2)], labels = genes_gtf$gene_name[target_gene_entry], cex = text_cex) 
      }
      gene_idx <- gene_idx + 1
    }
  }
  if(pdf_file != "") dev.off()
}
