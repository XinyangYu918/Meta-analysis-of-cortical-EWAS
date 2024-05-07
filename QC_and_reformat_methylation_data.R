#Prep1: methylation data QC, at the CpG level

# Probe Quality Check
library(minfi)

# Prep1: Probe Quality Check
# Load files
load('./RGset.rda')
load('./Quan-norm.rda')

# Add SNP information to the data
objectWithSNPinfo <- addSnpInfo(object)

# Drop probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension
objectSNPQCed <- dropLociWithSnps(objectWithSNPinfo, snps = c("SBE", "CpG", "Probe"), maf = 0.05)
rm(objectWithSNPinfo)

# A detection p-value is returned for every genomic position in every sample
detP <- detectionP(RGset)
Match1 <- match(colnames(objectSNPQCed),colnames(detP))
Match2 <- match(rownames(objectSNPQCed),rownames(detP))
detPSNPQCed <- detP[Match2[!is.na(Match2)],Match1[!is.na(Match1)]]
rm(Match1, Match2, detP)

# Positions with non-significant p-values (typically >0.01) should not be trusted
failed <- detPSNPQCed >0.01
rm(detPSNPQCed)

# Get beta values
beta <- getBeta(objectSNPQCed)

# Drop probes that failed quality control via the detection p-value in greater than 20% of samples
failedCG02 <- rowMeans(failed) >0.2

# Get the list of non-variable CpG sites i.e. those where beta values for all samples are ≤20% or ≥80%
ProbeInvar <- (rowSums(beta<=0.2)==ncol(beta))|(rowSums(beta>=0.8)==ncol(beta)) 

# Mark probes with either all beta value <=0.2 or all beta value>=0.80 and generate a list of probes marked as invariant
# This list will be included in the output file that you will send us
ListInvarProbe <- rownames(beta)[which(ProbeInvar)] 
rm(ProbeInvar)

# Remove sex chromosome probes
keepIndex=!seqnames(objectSNPQCed)%in%c("chrX","chrY") 
rm(objectSNPQCed)

# Remove failed probes
# Combine with failed probes in >20% of samples
keepIndex <- keepIndex&(!failedCG02)
rm(failedCG02)

# Remove all probes with detected P-value >0.01
beta[failed] <- NA 
rm(failed)

# Remove marked failed probes 
betaQC <- beta[which(keepIndex),] 
rm(beta, keepIndex)

# Prep2: Reformat beta value for EWAS
# Load quantile normalized methylation data
Methy <- as.data.frame(t(betaQC))

# Add subject ID to the methylation data as the last column
# You may need to modify the codes if your methylation data are labeled using lab-specific codes (such as slide and array IDs).
# Ensure that the "Subject" refers to your participant ID, which should correspond with your MRI and covariate data.
Methy$Subject <- colnames(betaQC) 
rm(betaQC)
gc()

# Re-order the column and making the subject ID as the first column
Methy <- Methy[,c(ncol(Methy),1:(ncol(Methy)-1))] 
# Save probe names
MethyName <- colnames(Methy)[-c(1)] 

# Save QCed methylation data for EWAS (We recommend saving this data locally at your site.)
save(Methy, file = "Methy.RData")

# Prep3: format methylation covariates, including the first four PCs of beta values and the first 2 PCs of cell type
# Load principal components of beta value
load("./fast_svd.rda") 
# Generate a variable with first four principal components of beta value
PC_Beta <- as.data.frame(ss$v[,1:4]) 
# Add subject ID to the component variable
# Ensure that the "Subject" refers to your participant ID, which should correspond with your MRI and covariate data.
PC_Beta$Subject <- Methy$Subject  

# Load estimated cell type proportion
load("./cellcount.rda") 
# Generate principal components of cell count
tmp <- prcomp(cellcount) 
# Generate a variable with subject ID and first 2 PCs of cell count
pc_cell <- as.data.frame(tmp$x[,1:2]) 
# Add subject ID to the component variable
# Ensure that the "Subject" refers to your participant ID, which should correspond with your MRI and covariate data.
pc_cell$Subject <- rownames(pc_cell)
