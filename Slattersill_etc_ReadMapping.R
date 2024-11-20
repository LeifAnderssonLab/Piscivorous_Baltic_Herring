#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se



#Slåttersill etc data processing
#find /proj/snic2020-2-19/private/herring/data_deliveries/snpseq00388/DataDelivery_2023-06-20_09-28-45_snpseq00388/files/WE-3675/230609_A00605_0577_AH2VKNDSX7 -maxdepth 5 -type f -name '*.fastq.gz' | grep "R1_00" > /proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/mapping/Slatter_etc_R1_files.txt
#find /proj/snic2020-2-19/private/herring/data_deliveries/snpseq00388/DataDelivery_2023-06-20_09-28-45_snpseq00388/files/WE-3675/230609_A00605_0577_AH2VKNDSX7 -maxdepth 5 -type f -name '*.fastq.gz' | grep "R2_00" > /proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/mapping/Slatter_etc_R2_files.txt

##
R1_files <- read.table("~/Projects/Herring/data/slattersill_etc/mapping/Slatter_etc_R1_files.txt", stringsAsFactors = F)
R1_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_00[0-9].fastq.gz", "\\1", R1_files[,1])
R1_files[,"Short_ID"] <- sub("_S[0-9]+_L[0-9]+", "", R1_files[,"ID"])
R1_files[,"Short_ID"] <- sub("WE[-]3675[-]", "", R1_files[,"Short_ID"])
R2_files <- read.table("~/Projects/Herring/data/slattersill_etc/mapping/Slatter_etc_R2_files.txt", stringsAsFactors = F)
R2_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_00[0-9].fastq.gz", "\\1", R2_files[,1])

all(R2_files$ID %in% R1_files$ID)
#[1] TRUE
all(R1_files$ID %in% R2_files$ID)
#[1] TRUE
#All IDs accounted for

all(R1_files$ID == R2_files$ID)
#[1] TRUE
#And in order

R2_files <- R2_files[match(R1_files$ID, R2_files$ID),]
all(R1_files$ID == R2_files$ID)
#[1] 
all(R1_files$Ind == R2_files$Ind)
#[1] 
#
write(R1_files[,"Short_ID"], file = "~/Projects/Herring/data/slattersill_etc/mapping/Slatter_etc_sample_IDs.txt")

##Second batch
#find /proj/snic2020-2-19/private/herring/data_deliveries/snpseq00453/DataDelivery_2023-08-24_14-39-45_snpseq00453/files/WF-3715/ -maxdepth 5 -type f -name '*.fastq.gz' | grep "R1_00" > /proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/mapping/Switchers_etc_R1_files.txt
#find /proj/snic2020-2-19/private/herring/data_deliveries/snpseq00453/DataDelivery_2023-08-24_14-39-45_snpseq00453/files/WF-3715/ -maxdepth 5 -type f -name '*.fastq.gz' | grep "R2_00" > /proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/mapping/Switchers_etc_R2_files.txt

##
R1_files <- read.table("~/Projects/Herring/data/slattersill_etc/mapping/Switchers_etc_R1_files.txt", stringsAsFactors = F)
R1_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_00[0-9].fastq.gz", "\\1", R1_files[,1])
R1_files[,"Short_ID"] <- sub("_S[0-9]+_L[0-9]+", "", R1_files[,"ID"])
R1_files[,"Short_ID"] <- sub("WF[-]3715[-]", "", R1_files[,"Short_ID"])
R2_files <- read.table("~/Projects/Herring/data/slattersill_etc/mapping/Switchers_etc_R2_files.txt", stringsAsFactors = F)
R2_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_00[0-9].fastq.gz", "\\1", R2_files[,1])

all(R2_files$ID %in% R1_files$ID)
#[1] TRUE
all(R1_files$ID %in% R2_files$ID)
#[1] TRUE
#All IDs accounted for

all(R1_files$ID == R2_files$ID)
#[1] TRUE
#And in order
write(R1_files[,"Short_ID"], file = "~/Projects/Herring/data/slattersill_etc/mapping/Switchers_etc_sample_IDs.txt")

#Filtration etc

#zgrep "PASS" SlatterSill_SpawnSwitchers_WGS.filter.vcf.gz | wc -l
#35023126

Slatter_etc_metadata <- read.table(file = "~/Projects/Herring/data/slattersill_etc/Metadata_WGS_SlatterSill_SpawnSwitchers.txt", header = T, sep = "\t", stringsAsFactors = F)
write(x = Slatter_etc_metadata$Sample_ID[Slatter_etc_metadata$Classification %in% c("Large_Herring", "SlåtterSill")], file = "~/Projects/Herring/data/slattersill_etc/All_large_herring_inds.txt")

Slatter_etc_metadata_class <- read.table(file = "~/Projects/Herring/data/slattersill_etc/Metadata_WGS_SlatterSill_SimilarityGroups.txt", header = T, sep = "\t", stringsAsFactors = F) #All large herring samples

slatter_group_IDs <- unique(Slatter_etc_metadata_class$Group)

write(x = Slatter_etc_metadata_class$Sample_ID[grep("Gävleborg_[SI]", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Gavle_large_herring_inds.txt") #Large herring from the Gävle area,excluding original Slåttersill
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Gävleborg_[N]", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Gavle_regular_herring_inds.txt") #Control herring from the Gävle area
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Gävlebukten", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Slattersill_inds.txt") #Original Slåttersill
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Blekinge", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Blekinge_inds.txt") #Blekinge
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Finland", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Finland_inds.txt") #Finland
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Kalmar", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Kalmar_inds.txt") #Kalmar
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Östergötland", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Ostergotland_inds.txt") #Östergötland
write(x = Slatter_etc_metadata_class$Sample_ID[grep("Stockholm", Slatter_etc_metadata_class$Group)], file = "~/Projects/Herring/data/slattersill_etc/Stockholm_inds.txt") #Stockholm

#Spawnswitcher groups
write(x = Slatter_etc_metadata$Sample_ID[grep("GSAutumnSpawner", Slatter_etc_metadata$Classification)], file = "~/Projects/Herring/data/slattersill_etc/GS_AutumnSpawner_inds.txt") #Genetic spring-spawners spawning in autumn
write(x = Slatter_etc_metadata$Sample_ID[grep("^AutumnSpawner", Slatter_etc_metadata$Classification)], file = "~/Projects/Herring/data/slattersill_etc/AutumnSpawner_regular_inds.txt") #Regular autumn-spawners
write(x = Slatter_etc_metadata$Sample_ID[grep("GASpringSpawner", Slatter_etc_metadata$Classification)], file = "~/Projects/Herring/data/slattersill_etc/GA_SpringSpawner_inds.txt") #Genetic autumn-spawners spawning in spring
write(x = Slatter_etc_metadata$Sample_ID[grep("^SpringSpawner", Slatter_etc_metadata$Classification)], file = "~/Projects/Herring/data/slattersill_etc/SpringSpawner_regular_inds.txt") #Regular spring-spawners

#Getting to grips with surprising pattern
snp_chip_sanity_check <- read.table("~/Projects/Herring/data/slattersill_etc/Temp_Sanity_Check_SNPChip_AF_Chrom10_to_12.txt", header = T)
snp_chip_sanity_check$CHROM <- paste0("chr", snp_chip_sanity_check$CHROM)
sanity_check_merged <- inner_join(x = snp_chip_sanity_check, merged_freq, by = c("CHROM", "POS"))

#Looking at individual genotypes
sanity_check_ind <- read.table("~/Projects/Herring/data/slattersill_etc/Combined_Genotypes_WGS_SNPChip.txt", header = T, stringsAsFactors = F)


CEL_cols <- grep("a[.]CEL", names(sanity_check_ind), value = T)
Gav20_cor_df <- data.frame(row.names = CEL_cols)
Gav20_WGS_col <- grep("Gav20[.][0-9]+$", names(sanity_check_ind), value = T)
Upp05_cor_df <- data.frame(row.names = CEL_cols)
Upp05_WGS_col <- grep("Upp05[.][0-9]+$", names(sanity_check_ind), value = T)
Sth02_cor_df <- data.frame(row.names = CEL_cols)
Sth02_WGS_col <- grep("Sth02[.][0-9]+$", names(sanity_check_ind), value = T)
Gav19_cor_df <- data.frame(row.names = CEL_cols)
Gav19_WGS_col <- grep("Gav19[.][0-9]+$", names(sanity_check_ind), value = T)


#Gav20_3_correlations <- numeric()
#Gav20_74_correlations <- numeric()
for(WGS_col in Gav20_WGS_col){
  Gav20_cor_df[,WGS_col] <- NA
  for(CEL_col in CEL_cols){
    Gav20_cor_df[CEL_col,WGS_col]<- cor(x = sanity_check_ind[,WGS_col], y = sanity_check_ind[,CEL_col], use = "complete.obs")
    #Gav20_3_correlations[i] <- cor(x = sanity_check_ind$Gav20.3, y = sanity_check_ind[,CEL_cols[i]], use = "complete.obs")
    #names(Gav20_3_correlations)[i] <- CEL_cols[i]
    #Gav20_74_correlations[i] <- cor(x = sanity_check_ind$Gav20.74, y = sanity_check_ind[,CEL_cols[i]], use = "complete.obs")
    #names(Gav20_74_correlations)[i] <- CEL_cols[i]
  }
}
for(WGS_col in Gav19_WGS_col){
  Gav19_cor_df[,WGS_col] <- NA
  for(CEL_col in CEL_cols){
    Gav19_cor_df[CEL_col,WGS_col]<- cor(x = sanity_check_ind[,WGS_col], y = sanity_check_ind[,CEL_col], use = "complete.obs")
  }
}

for(WGS_col in Upp05_WGS_col){
  Upp05_cor_df[,WGS_col] <- NA
  for(CEL_col in CEL_cols){
    Upp05_cor_df[CEL_col,WGS_col]<- cor(x = sanity_check_ind[,WGS_col], y = sanity_check_ind[,CEL_col], use = "complete.obs")
  }
}
for(WGS_col in Sth02_WGS_col){
  Sth02_cor_df[,WGS_col] <- NA
  for(CEL_col in CEL_cols){
    Sth02_cor_df[CEL_col,WGS_col]<- cor(x = sanity_check_ind[,WGS_col], y = sanity_check_ind[,CEL_col], use = "complete.obs")
  }
}


#Other direction
WGS_cols <- grep("[A-Za-z0-9]+[.][0-9]+$", names(sanity_check_ind), value = T)
Gav20_WGS_cor_df <- data.frame(row.names = WGS_cols)
Gav20_CEL_col <- grep("Gav20_[0-9]+_a[.]CEL", names(sanity_check_ind), value = T)

for(CEL_col in Gav20_CEL_col){
  Gav20_WGS_cor_df[,CEL_col] <- NA
  for(WGS_col in WGS_cols){
    Gav20_WGS_cor_df[WGS_col,CEL_col]<- cor(x = sanity_check_ind[,WGS_col], y = sanity_check_ind[,CEL_col], use = "complete.obs")
  }
}

##bcftools commantds for renaming chromsosme ins SNP-chip vcf.
#bcftools annotate -o ./Merged_Gavlebukten_Slattersill_2023_2024_Sampling_chr.vcf.gz --rename-chrs ./chr_rename_list.txt ./Merged_Gavlebukten_Slattersill_2023_2024_Sampling_rename_chr.vcf.gz 
#bcftools view -i 'CHROM~"chr"' -o ./Merged_Gavlebukten_Slattersill_2023_2024_Sampling_chr_only.vcf.gz ./Merged_Gavlebukten_Slattersill_2023_2024_Sampling_chr.vcf.gz 
#bcftools index --tbi ./Merged_Gavlebukten_Slattersill_2023_2024_Sampling_chr_only.vcf.gz 


