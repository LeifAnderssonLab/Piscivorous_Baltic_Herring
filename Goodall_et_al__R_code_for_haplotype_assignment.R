# The following will make haplotype assignments based on Chr 6,12,17,23

# 1. Prepare the data for Haplotype assignment #
################################################

# Set working directory
setwd("/pathtodata/Raw_Data_vcf/")

# Import the vcf file
library(vcfR)
raw_vcf <- vcfR::read.vcfR("Samples_to_haplotype_assign.vcf")

# Combine genotype data and fixed metadata
combined_data <- cbind(raw_vcf@fix, raw_vcf@gt)

# Subset columns 1:9 and 10:573 separately
subset_columns1 <- combined_data[, 1:9] # this should subset all the vcf metadata
subset_columns2 <- combined_data[, 10:(ncol(combined_data))] # This should just be the genotype data

# Apply the conversion function to columns 10:573
converted_matrix <- apply(subset_columns2, c(1, 2), function(x) {
  if (is.na(x)) {
    return(NA)  # Handle missing values separately
  } else if (x == "0/0") {
    return(0)
  } else if (x == "0/1") {
    return(1)
  } else if (x == "1/1") {
    return(2)
  } else {
    return(as.numeric(x))
  }
})

# Combine columns 1:9 and converted_matrix using cbind
converted_matrix <- cbind(subset_columns1, converted_matrix) # combine the reformatted genotypes with the original metadata

# Subset out the chromosomes with inversions
chr6_subset_data <- converted_matrix[converted_matrix[, 1] == "6", ]
chr12_subset_data <- converted_matrix[converted_matrix[, 1] == "12", ]
chr17_subset_data <- converted_matrix[converted_matrix[, 1] == "17", ]
chr23_subset_data <- converted_matrix[converted_matrix[, 1] == "23", ]

# Extract out the inversion region from Chr 6 
chr6_subset_data <- as.data.frame(chr6_subset_data)
chr6_subset_data_Inv_Only <- chr6_subset_data[chr6_subset_data$POS >= 22282765 & chr6_subset_data$POS <= 24868581, ]

# Extract out the inversion region from Chr 12 
chr12_subset_data <- as.data.frame(chr12_subset_data)
chr12_subset_data_Inv_Only <- chr12_subset_data[chr12_subset_data$POS >= 17826318 & chr12_subset_data$POS <= 25603093, ]

# Extract out the inversion region from Chr 17 
chr17_subset_data <- as.data.frame(chr17_subset_data)
chr17_subset_data_Inv_Only <- chr17_subset_data[chr17_subset_data$POS >= 25805445 & chr17_subset_data$POS <= 27568511, ]

# Extract out the inversion region from Chr 23 
chr23_subset_data <- as.data.frame(chr23_subset_data)
chr23_subset_data_Inv_Only <- chr23_subset_data[chr23_subset_data$POS >= 16226443 & chr23_subset_data$POS <= 17604273, ]



# 2. Haplotype assignment for Chr 6            #
################################################

# Chr 6 dataset
Chr6_R_Dataset <- read.delim2("/pathtodata/Chr6_R_Dataset.txt") # imports the reference North and South datasets for chr6 to compare against
Chr6_Haplo_Ref <- Chr6_R_Dataset[, 1:9]
Chr6_Haplo_Ref_Clean <- Chr6_Haplo_Ref[, -c(6,8)]

# Get the unique POS values from chr6_subset_data_Inv_Only
pos_subset_data <- unique(chr6_subset_data_Inv_Only$POS)

# Get the unique POS values from Chr6_Haplo_Ref_Clean
pos_ref_clean <- unique(Chr6_Haplo_Ref_Clean$POS)

# Find the matching POS values
matching_pos <- intersect(pos_subset_data, pos_ref_clean)

# Subset to only entries matching Chr6_REF
subsetted_data <- chr6_subset_data_Inv_Only[chr6_subset_data_Inv_Only$POS %in% matching_pos, ]

# Recreate the dataset with the required info
temp_1 <- subsetted_data[, 1:5]
temp_2 <- subsetted_data[, 10:573]
temp_3 <- Chr6_Haplo_Ref_Clean[, 6:7]

# Create the final working dataset for chromosome 6
Working_Chr6 <- cbind(temp_1, temp_3, temp_2)

# create empty data frame to store results
Chr6_Haplotype_Assign_Results <- data.frame(column_name = character(), sum_abs_diff = numeric(), n_diff_loci = numeric(), Haplo_Assign = character(), stringsAsFactors = FALSE)

# loop over columns and calculate sum of absolute differences
for (i in 8:571) {
  colname <- colnames(Working_Chr6)[i]
  coldata <- as.numeric(Working_Chr6[, colname])  # Convert coldata to numeric
  abs_diff <- sum(abs(as.numeric(Working_Chr6$SOUTH_Haplo[!is.na(Working_Chr6$SOUTH_Haplo)]) - coldata[!is.na(coldata)]))
  n_diff_loci <- abs_diff/2
  
  if (n_diff_loci < 18) {
    haplo_assign <- "SS"
  } else if (n_diff_loci > 34) {
    haplo_assign <- "NN"
  } else {
    haplo_assign <- "NS"
  }
  
  Chr6_Haplotype_Assign_Results[i-9,] <- c(colname, abs_diff, n_diff_loci, haplo_assign)
}

# print result to console
print(Chr6_Haplotype_Assign_Results)
write.table(Chr6_Haplotype_Assign_Results, "/pathtoexport/Chr6_Haplotype_Assign_Results.txt")

# 3. Haplotype assignment for Chr 12           #
################################################

# Chr 12 dataset
Chr12_R_Dataset <- read.delim2("/pathtodata/Chr12_R_Dataset.txt") # imports the reference North and South datasets for chr12 to compare against
Chr12_Haplo_Ref <- Chr12_R_Dataset[, 1:9]
Chr12_Haplo_Ref_Clean <- Chr12_Haplo_Ref[, -c(6,8)]

# Get the unique POS values from Chr12_subset_data_Inv_Only
pos_subset_data <- unique(chr12_subset_data_Inv_Only$POS)

# Get the unique POS values from Chr12_Haplo_Ref_Clean
pos_ref_clean <- unique(Chr12_Haplo_Ref_Clean$POS)

# Find the matching POS values
matching_pos <- intersect(pos_subset_data, pos_ref_clean)

# Subset to only entries matching Chr12_REF
subsetted_data <- chr12_subset_data_Inv_Only[chr12_subset_data_Inv_Only$POS %in% matching_pos, ]

# Recreate the dataset with the required info
temp_1 <- subsetted_data[, 1:5]
temp_2 <- subsetted_data[, 10:573]
temp_3 <- Chr12_Haplo_Ref_Clean[, 3:7]

# Create the final working dataset for chromosome 6
Working_Chr12 <- cbind(temp_1, temp_3, temp_2)


# create empty data frame to store results
Chr12_Haplotype_Assign_Results <- data.frame(column_name = character(), sum_abs_diff = numeric(), n_diff_loci = numeric(), Haplo_Assign = character(), stringsAsFactors = FALSE)

# loop over columns and calculate sum of absolute differences 
#### !! NOTE THERE IS AN ISSUE WITH THE VCFs REF/ALT POLARISATION, I HAVE REVERSED THE DESIGNATIONS HERE TO ACCOUNT FOR THIS ! #####
for (i in 11:574) {
  colname <- colnames(Working_Chr12)[i]
  coldata <- as.numeric(Working_Chr12[, colname])  # Convert coldata to numeric
  abs_diff <- sum(abs(as.numeric(Working_Chr12$SOUTH_Haplo[!is.na(Working_Chr12$SOUTH_Haplo)]) - coldata[!is.na(coldata)]))
  n_diff_loci <- abs_diff/2
  
  if (n_diff_loci < 38) {
    haplo_assign <- "NN"
  } else if (n_diff_loci > 76) {
    haplo_assign <- "SS"
  } else {
    haplo_assign <- "NS"
  }
  
  Chr12_Haplotype_Assign_Results[i-9,] <- c(colname, abs_diff, n_diff_loci, haplo_assign)
}

# print result to console
print(Chr12_Haplotype_Assign_Results)
write.table(Chr12_Haplotype_Assign_Results, "/pathtoexport/Chr12_Haplotype_Assign_Results.txt")


# 4. Haplotype assignment for Chr 17                #
#####################################################

# Chr 17 dataset
Chr17_R_Dataset <- read.delim2("/pathtodata/Chr17_R_Dataset.txt")
chr17_Haplo_Ref <- Chr17_R_Dataset[, 1:9]
chr17_Haplo_Ref_Clean <- chr17_Haplo_Ref[, -c(6,8)]

# Get the unique POS values from chr17_subset_data_Inv_Only
pos_subset_data <- unique(chr17_subset_data_Inv_Only$POS)

# Get the unique POS values from chr17_Haplo_Ref_Clean
pos_ref_clean <- unique(chr17_Haplo_Ref_Clean$POS)

# Find the matching POS values
matching_pos <- intersect(pos_subset_data, pos_ref_clean)

# Subset to only entries matching chr17_REF
subsetted_data <- chr17_subset_data_Inv_Only[chr17_subset_data_Inv_Only$POS %in% matching_pos, ]

# Recreate the dataset with the required info
temp_1 <- subsetted_data[, 1:5]
temp_2 <- subsetted_data[, 10:573]
temp_3 <- chr17_Haplo_Ref_Clean[, 3:7]

# Create the final working dataset for chromosome 6
Working_chr17 <- cbind(temp_1, temp_3, temp_2)


# create empty data frame to store results
chr17_Haplotype_Assign_Results <- data.frame(column_name = character(), sum_abs_diff = numeric(), n_diff_loci = numeric(), Haplo_Assign = character(), stringsAsFactors = FALSE)

# loop over columns and calculate sum of absolute differences 
#### !! NOTE THERE IS AN ISSUE WITH THE VCFs REF/ALT POLARISATION, I HAVE REVERSED THE DESIGNATIONS HERE TO ACCOUNT FOR THIS ! #####
for (i in 11:574) {
  colname <- colnames(Working_chr17)[i]
  coldata <- as.numeric(Working_chr17[, colname])  # Convert coldata to numeric
  abs_diff <- sum(abs(as.numeric(Working_chr17$SOUTH_Haplo[!is.na(Working_chr17$SOUTH_Haplo)]) - coldata[!is.na(coldata)]))
  n_diff_loci <- abs_diff/2
  
  if (n_diff_loci < 18) {
    haplo_assign <- "SS"
  } else if (n_diff_loci > 36) {
    haplo_assign <- "NN"
  } else {
    haplo_assign <- "NS"
  }
  
  chr17_Haplotype_Assign_Results[i-9,] <- c(colname, abs_diff, n_diff_loci, haplo_assign)
}

# print result to console
print(chr17_Haplotype_Assign_Results)
write.table(chr17_Haplotype_Assign_Results, "/pathtoexport/Chr17_Haplotype_Assign_Results.txt")


# 5. Haplotype assignment for Chr 23                #
#####################################################

# Chr 23 dataset
Chr23_R_Dataset <- read.delim2("/pathtodata/Chr23_R_Dataset.txt")
chr23_Haplo_Ref <- Chr23_R_Dataset[, 1:9]
chr23_Haplo_Ref_Clean <- chr23_Haplo_Ref[, -c(6,8)]

# Get the unique POS values from chr23_subset_data_Inv_Only
pos_subset_data <- unique(chr23_subset_data_Inv_Only$POS)

# Get the unique POS values from chr23_Haplo_Ref_Clean
pos_ref_clean <- unique(chr23_Haplo_Ref_Clean$POS)

# Find the matching POS values
matching_pos <- intersect(pos_subset_data, pos_ref_clean)

# Subset to only entries matching chr23_REF
subsetted_data <- chr23_subset_data_Inv_Only[chr23_subset_data_Inv_Only$POS %in% matching_pos, ]

# Recreate the dataset with the required info
temp_1 <- subsetted_data[, 1:5]
temp_2 <- subsetted_data[, 10:573]
temp_3 <- chr23_Haplo_Ref_Clean[, 3:7]

# Create the final working dataset for chromosome 6
Working_chr23 <- cbind(temp_1, temp_3, temp_2)

# create empty data frame to store results
chr23_Haplotype_Assign_Results <- data.frame(column_name = character(), sum_abs_diff = numeric(), n_diff_loci = numeric(), Haplo_Assign = character(), stringsAsFactors = FALSE)

# loop over columns and calculate sum of absolute differences 
#### !! NOTE THERE IS AN ISSUE WITH THE VCFs REF/ALT POLARISATION, I HAVE REVERSED THE DESIGNATIONS HERE TO ACCOUNT FOR THIS ! #####
for (i in 11:574) {
  colname <- colnames(Working_chr23)[i]
  coldata <- as.numeric(Working_chr23[, colname])  # Convert coldata to numeric
  abs_diff <- sum(abs(as.numeric(Working_chr23$SOUTH_Haplo[!is.na(Working_chr23$SOUTH_Haplo)]) - coldata[!is.na(coldata)]))
  n_diff_loci <- abs_diff/2
  
  if (n_diff_loci < 0.5) {
    haplo_assign <- "SS"
  } else if (n_diff_loci > 0.9) {
    haplo_assign <- "NN"
  } else {
    haplo_assign <- "NS"
  }
  
  chr23_Haplotype_Assign_Results[i-9,] <- c(colname, abs_diff, n_diff_loci, haplo_assign)
}

# print result to console
print(chr23_Haplotype_Assign_Results)
write.table(chr23_Haplotype_Assign_Results, "/pathtoexport/Chr23_Haplotype_Assign_Results.txt")


# 6. Visualising Haplotype assignments              #
#####################################################
library(pheatmap)

# Chrom 6 Visualized
#View(Working_Chr6)
datamatrix <- Working_Chr6[, -(1:7)]
datamatrix_2 <- datamatrix[, -(1:247)]
numeric_dataset <- sapply(datamatrix_2, function(x) as.numeric(as.character(x)))
is.numeric(numeric_dataset)

chr6_heatmap <- pheatmap(t(numeric_dataset), show_colnames = F, labels_col= Working_Chr6$ID, fontsize = 6, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color= NA, color= c("gold", "pink", "midnightblue"), show_rownames = T, cluster_cols = F, cluster_rows = F, legend_breaks = c(0, 1, 2), legend_labels = c("NN", "NS", "SS"))

# Chrom 12 Visualized
# NOTE: This one has an inverted colour scheme because of the polarisation issues for this chromosome
#View(Working_Chr12)
datamatrix <- Working_Chr12[, -(1:10)]
datamatrix_2 <- datamatrix[, -(1:247)]
numeric_dataset <- sapply(datamatrix_2, function(x) as.numeric(as.character(x)))
is.numeric(numeric_dataset)

chr12_heatmap <- pheatmap(t(numeric_dataset), show_colnames = F, labels_col= Working_Chr12$ID, fontsize = 6, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color= NA, color= c("midnightblue", "pink", "gold"), show_rownames = T, cluster_cols = F, cluster_rows = T, legend_breaks = c(0, 1, 2), legend_labels = c("SS", "NS", "NN"))

# Chrom 17 Visualized
#View(Working_chr17_PhaseCorrected)
datamatrix <- Working_chr17[, -(1:10)]
datamatrix_2 <- datamatrix[, -(1:247)]
numeric_dataset <- sapply(datamatrix_2, function(x) as.numeric(as.character(x)))
is.numeric(numeric_dataset)

chr17_heatmap <- pheatmap(t(numeric_dataset), show_colnames = F, labels_col= Working_chr17$ID, fontsize = 6, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color= NA, color= c("gold", "pink", "midnightblue"), show_rownames = T, cluster_cols = F, cluster_rows = F, legend_breaks = c(0, 1, 2), legend_labels = c("NN", "NS", "SS"))

# Chrom 23 Visualized
#View(Working_chr23)
datamatrix <- Working_chr23[, -(1:10)]
datamatrix_2 <- datamatrix[, -(1:247)]
numeric_dataset <- sapply(datamatrix_2, function(x) as.numeric(as.character(x)))
is.numeric(numeric_dataset)

# For some reason it gets weird about the labelling, given there is only one SNP. Should let it plot the labels itself
chr23_heatmap <- pheatmap(t(numeric_dataset), show_colnames = T, fontsize = 6, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", border_color= NA, color= c("gold", "pink", "midnightblue"), show_rownames = T, cluster_cols = F, cluster_rows = F, legend_breaks = c(0, 1, 2), legend_labels = c("NN", "NS", "SS"))



# 7. Automate the removal of inversions and replacement with haplo assignments              #
#############################################################################################

replace_matrix <- as.data.frame(converted_matrix)

# Remove rows based on specified coordinates for Chr 6
replace_matrix <- replace_matrix[!(replace_matrix$CHROM == "6" & replace_matrix$POS >= 22282765 & replace_matrix$POS <= 24868581), ]

# Remove rows based on specified coordinates for Chr 12
replace_matrix <- replace_matrix[!(replace_matrix$CHROM == "12" & replace_matrix$POS >= 17826318 & replace_matrix$POS <= 25603093), ]

# Remove rows based on specified coordinates for Chr 17
replace_matrix <- replace_matrix[!(replace_matrix$CHROM == "17" & replace_matrix$POS >= 25805445 & replace_matrix$POS <= 27568511), ]

# Remove rows based on specified coordinates for Chr 23
replace_matrix <- replace_matrix[!(replace_matrix$CHROM == "23" & replace_matrix$POS >= 16226443 & replace_matrix$POS <= 17604273), ]

str(replace_matrix)
str(Chr6_Haplotype_Assign_Results)

# Remove specific samples if needed
replace_matrix_2 <- replace_matrix[, -c(10:256)]


# Extract Haplo_Assign column from Chr6_Haplotype_Assign_Results
Chr6_Haplotype_Assign_Results_2 <- Chr6_Haplotype_Assign_Results[-c(1:245) ,]
haplo_assignments_chr6 <- Chr6_Haplotype_Assign_Results_2$Haplo_Assign

# Extract Haplo_Assign column from Chr6_Haplotype_Assign_Results
Chr12_Haplotype_Assign_Results_2 <- Chr12_Haplotype_Assign_Results[-c(1:248) ,]
haplo_assignments_chr12 <- Chr12_Haplotype_Assign_Results_2$Haplo_Assign

# Extract Haplo_Assign column from Chr6_Haplotype_Assign_Results
Chr17_Haplotype_Assign_Results_2 <- chr17_Haplotype_Assign_Results[-c(1:248) ,]
haplo_assignments_chr17 <- Chr17_Haplotype_Assign_Results_2$Haplo_Assign

# Extract Haplo_Assign column from Chr6_Haplotype_Assign_Results
Chr23_Haplotype_Assign_Results_2 <- chr23_Haplotype_Assign_Results[-c(1:248) ,]
haplo_assignments_chr23 <- Chr23_Haplotype_Assign_Results_2$Haplo_Assign

# Assign Haplo_Assign values to columns starting from the 10th column in replace_matrix
replace_matrix_2["Chr6_Haplo", 10:ncol(replace_matrix_2)] <- haplo_assignments_chr6 # will be added to the bottom of the dataframe
replace_matrix_2["Chr12_Haplo", 10:ncol(replace_matrix_2)] <- haplo_assignments_chr12 # will be added to the bottom of the dataframe
replace_matrix_2["Chr17_Haplo", 10:ncol(replace_matrix_2)] <- haplo_assignments_chr17 # will be added to the bottom of the dataframe
replace_matrix_2["Chr23_Haplo", 10:ncol(replace_matrix_2)] <- haplo_assignments_chr23 # will be added to the bottom of the dataframe


write.table(replace_matrix_2, "ALLSamples_Haplotype_Reduced_RAW.txt")

