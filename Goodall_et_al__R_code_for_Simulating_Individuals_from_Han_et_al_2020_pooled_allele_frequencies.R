# Import the dataset
####################################################################################################
## Importing the data
setwd("~/Documents/Herring_Array/Herring_VCF_SNP_DESIGN_ORIGINAL/") # change to your file location

Han_allele_freq <- read.delim("Han_allele_freq_TRIMMED_TO_INCLUDE_ONLY_ARRAY_SITES.txt") # pooled allele frequency data was filtered to only contain loci from the SNPChip design
colnames(Han_allele_freq) # get the name of all of the samples
Han_allele_freq$SNP_ID # get the name of all of the SNP_IDs in order


# Run the simulation
####################################################################################################
## Subset the Han et al. 2020 dataset to only include Atlantic herring
base_herring_matrix <- Han_allele_freq[c(4:38,45:58,60:63)] # subset out only the locations that you want

## Create an empty matrix that you can append the simulated data to
Simulated_Matrix <- matrix(nrow = length(base_herring_matrix[,1]), ncol = 1) #nrow has to match the number of alleles in base_herring_matrix and ncol=number of individuals to simulate
Simulated_Matrix

# Run a loop to simulate all data across the allele table, then rejoin it back into a table with the appropriate names

for(v in 1:ncol(base_herring_matrix)) {       # for-loop over columns

## Get some information on missing values in the allele values
sum(is.na(base_herring_matrix[ , v])) # tell me how many NAs present 
SUM.NAs <- sum(is.na(base_herring_matrix[ , v])) # store how many NAs present
which(is.na(base_herring_matrix[ , v])) # tell me which are NAs
Is_NA <- which(is.na(base_herring_matrix[ , v])) # store the NAs locations for later use
base_herring_matrix[ , v][is.na(base_herring_matrix[ , v])] <- 1 # recode all of the NAs to 1
sum(is.na(base_herring_matrix[ , v])) # check how many NAs remaining after the recode (HINT: should be zero)


## Simulate a 8 individual dataset using the reference allele values from the Han et al 2020 paper
Temp_Sample_Matrix <- matrix(nrow = length(base_herring_matrix[,1]), ncol = 8) #nrow has to match the number of alleles in base_herring_matrix[ , i] and ncol=number of individuals to simulate
for(j in 1:dim(Temp_Sample_Matrix)[2]){
  for(i in 1:length(base_herring_matrix[ , v])){
    p0 <- base_herring_matrix[ , v][i]^2
    p2 <- (1-base_herring_matrix[ , v][i])^2
    p1 <- 1 -(p0 + p2)
    Temp_Sample_Matrix[i,j] <- 2
    if(runif(1) <= p0){
      Temp_Sample_Matrix[i,j] <- 0
    }
    if (runif(1) <= p0 + p1 & runif(1) > p0){
      Temp_Sample_Matrix[i,j] <- 1
    }
  }
}

head(Temp_Sample_Matrix) # I want to see the outcome of the above function

## Re-Label the rows/columns
numbering<-1:8 # this has the match the number of simulated individuals you've used
features <- c(sprintf(paste0(colnames(base_herring_matrix)[v],"_",numbering))) # make a sequential list of sample names
colnames(Temp_Sample_Matrix) <- features # Assign the names across the new array
head(Temp_Sample_Matrix) # check they were implemented properly

rownames(Temp_Sample_Matrix)<-Han_allele_freq$SNP_ID # reassign the SNP_IDs to the simulated SNP
head(Temp_Sample_Matrix) # check they were implemented properly

## Recreate NAs at the appropriate sites
Temp_Sample_Matrix[(Is_NA),] <- NA
SUM.NAs * 8 # The original number of NAs times 8 (as we are adding 8 samples)
sum(is.na(Temp_Sample_Matrix)) # number of NAs in our simulated dataset, this must match the above
Is_NA # List out which SNP were NAs in the original dataset
which(is.na(Temp_Sample_Matrix)) # list which SNP are NAs in our simulated dataset, this should match the above
#View(Temp_Sample_Matrix)

Simulated_Matrix <- cbind(Simulated_Matrix, Temp_Sample_Matrix)

}


# Finalise the dataset
############################################################################################################################
## Generate the final working dataset

Final_Simulated_Matrix <- Simulated_Matrix[,-1]
View(Final_Simulated_Matrix)

