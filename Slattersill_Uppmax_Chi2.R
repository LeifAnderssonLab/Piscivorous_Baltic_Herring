#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#For running on UPPMAX, allele-count based chi-squared tests
require(tidyverse)

#load count data
load("/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Han_et_al_counts.RData")
load("/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/merged_counts.RData")
#load("~/Projects/Herring/data/slattersill_etc/counts/Han_et_al_counts.RData")
#load("~/Projects/Herring/data/slattersill_etc/counts/merged_counts.RData")

#Define populations for contrasts
Baltic_spr_pops <- c("A_Kalix_Baltic_Spring", "B_Vaxholm_Baltic_Spring", "G_Gamleby_Baltic_Spring", "HGS1_Riga_Baltic_Spring", "HGS2_Riga_Baltic_Spring", "PB11_Kalmar_Baltic_Spring", "PB12_Karlskrona_Baltic_Spring", "PB1_HastKar_Baltic_Spring", "PB4_Hudiksvall_Baltic_Spring", "PB5_Galve_Baltic_Spring", "PN3_CentralBaltic_Baltic_Spring", "PB6_Galve_Baltic_Summer")
large_pops <- c("Blekinge", "Finland", "Gavle_large", "Kalmar", "Ostergotland", "Slattersill", "Sthlm")

#Merge counts
combined_counts <- inner_join(merged_count_df[,c("CHROM", "POS", large_pops)], Han_et_al_counts[,c("CHROM", "POS", Baltic_spr_pops)], by = c("CHROM", "POS"))

Large_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = large_pops, pop2 = Baltic_spr_pops)
save(Large_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Large_v_Spr_pop_chisq.RData")

Sthlm_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Sthlm", pop2 = Baltic_spr_pops)
save(Sthlm_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Sthlm_v_Spr_pop_chisq.RData")

Kalmar_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Kalmar", pop2 = Baltic_spr_pops)
save(Kalmar_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Kalmar_v_Spr_pop_chisq.RData")

Ostergotland_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Ostergotland", pop2 = Baltic_spr_pops)
save(Ostergotland_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Ostergotland_v_Spr_pop_chisq.RData")

Slattersill_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Slattersill", pop2 = Baltic_spr_pops)
save(Slattersill_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Slattersill_v_Spr_pop_chisq.RData")

Gavle_large_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Gavle_large", pop2 = Baltic_spr_pops)
save(Gavle_large_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Gavle_large_v_Spr_pop_chisq.RData")

Blekinge_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Blekinge", pop2 = Baltic_spr_pops)
save(Blekinge_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Blekinge_v_Spr_pop_chisq.RData")

Finland_v_Spr_pop_chisq <- pop_allele_chisq(count_df = combined_counts, pop1 = "Finland", pop2 = Baltic_spr_pops)
save(Finland_v_Spr_pop_chisq, file = "/proj/snic2020-2-19/private/herring/users/mats/slattersill_etc_2023/variants/counts/Finland_v_Spr_pop_chisq.RData")


#Support functions
pop_allele_chisq <- function(count_df, pop1, pop2){
  out_df <- count_df[,c("CHROM", "POS")]
  out_df[,"pop1_A"] <- 0
  out_df[,"pop1_D"] <- 0
  for(p1_col in pop1){
    out_df[,"pop1_A"] <- out_df[,"pop1_A"] + as.integer(sub("([0-9]+),([0-9]+)", "\\1", count_df[,p1_col]))
    out_df[,"pop1_D"] <- out_df[,"pop1_D"] + as.integer(sub("([0-9]+),([0-9]+)", "\\2", count_df[,p1_col]))
  }
  
  out_df[,"pop2_A"] <- 0
  out_df[,"pop2_D"] <- 0
  for(p2_col in pop2){
    out_df[,"pop2_A"] <- out_df[,"pop2_A"] + as.integer(sub("([0-9]+),([0-9]+)", "\\1", count_df[,p2_col]))
    out_df[,"pop2_D"] <- out_df[,"pop2_D"] + as.integer(sub("([0-9]+),([0-9]+)", "\\2", count_df[,p2_col]))
  }
  out_df[,"total"] <- out_df[,"pop1_A"] + out_df[,"pop1_D"] + out_df[,"pop2_A"] + out_df[,"pop2_D"]
  out_df[,"pop1_A_exp"] <- ((out_df[,"pop1_A"] + out_df[,"pop2_A"])/out_df[,"total"]) * (out_df[,"pop1_A"] + out_df[,"pop1_D"])
  out_df[,"pop1_D_exp"] <- ((out_df[,"pop1_D"] + out_df[,"pop2_D"])/out_df[,"total"]) * (out_df[,"pop1_A"] + out_df[,"pop1_D"])
  out_df[,"pop2_A_exp"] <- ((out_df[,"pop1_A"] + out_df[,"pop2_A"])/out_df[,"total"]) * (out_df[,"pop2_A"] + out_df[,"pop2_D"])
  out_df[,"pop2_D_exp"] <- ((out_df[,"pop1_D"] + out_df[,"pop2_D"])/out_df[,"total"]) * (out_df[,"pop2_A"] + out_df[,"pop2_D"])
  
  chisq_site_filter <- (out_df[,"pop1_A_exp"] > 0 | out_df[,"pop2_A_exp"] > 0) & 
    (out_df[,"pop1_D_exp"] > 0 | out_df[,"pop2_D_exp"] > 0) & 
    (out_df[,"pop1_A_exp"] + out_df[,"pop1_D_exp"]) > 10*length(pop1) & 
    (out_df[,"pop2_A_exp"] + out_df[,"pop2_D_exp"]) > 10*length(pop2)   
  
  
  out_df <-  out_df[chisq_site_filter,]
  
  #Continuity correction
  pop1_A_diff <- abs(out_df[,"pop1_A"] - out_df[,"pop1_A_exp"]) - 0.5
  pop1_A_diff[pop1_A_diff < 0] <- 0
  pop1_D_diff <- abs(out_df[,"pop1_D"] - out_df[,"pop1_D_exp"]) - 0.5
  pop1_D_diff[pop1_D_diff < 0] <- 0
  pop2_A_diff <- abs(out_df[,"pop2_A"] - out_df[,"pop2_A_exp"]) - 0.5
  pop2_A_diff[pop2_A_diff < 0] <- 0
  pop2_D_diff <- abs(out_df[,"pop2_D"] - out_df[,"pop2_D_exp"]) - 0.5
  pop2_D_diff[pop2_D_diff < 0] <- 0
  
  out_df[,"chisq_stat"] <-  ((pop1_A_diff^2)/out_df[,"pop1_A_exp"] + 
                               (pop1_D_diff^2)/out_df[,"pop1_D_exp"] +
                               (pop2_A_diff^2)/out_df[,"pop2_A_exp"] + 
                               (pop2_D_diff^2)/out_df[,"pop2_D_exp"]) 
  #Uncorrectred version
  #out_df[,"chisq_stat"] <-  (((out_df[,"pop1_A"] - out_df[,"pop1_A_exp"])^2)/out_df[,"pop1_A_exp"] + 
  #                             ((out_df[,"pop1_D"] - out_df[,"pop1_D_exp"])^2)/out_df[,"pop1_D_exp"] +
  #                             ((out_df[,"pop2_A"] - out_df[,"pop2_A_exp"])^2)/out_df[,"pop2_A_exp"] + 
  #                             ((out_df[,"pop2_D"] - out_df[,"pop2_D_exp"])^2)/out_df[,"pop2_D_exp"]) 
  
  out_df[,"chisq_p"] <- pchisq(q = out_df[,"chisq_stat"], df = 1, lower.tail = F)
  return(out_df)
}

