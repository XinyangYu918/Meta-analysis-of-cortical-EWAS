# Section 2: EWAS with cortical thickness, and surface areas
# Provide your cohort’s name
cohort = "IMAGEN"

# Get cortical measures
source("./prep_cortical_measures.R")

# Save column names of cortical surface measures, including thickness and area
CorticalMeasureName <- colnames(CorticalMeasure)[-1]

# Read and merge covariates
RawCov <- read.csv("./CorticalCovariates.csv",header=T)

# Assume the IID that is equal to the subject ID from cortical and methylation data
RawCov <- RawCov[,-c(1)]

# Calculate age^2 and add to the last column of Raw_Cov
RawCov$Age_Square <- (RawCov$Age)^2 

# Re-name the ID column for easier merge
names(RawCov)[1] <- c("Subject")

# Merge all covariates into one dataset
Cov <- merge(RawCov,PC_Beta,by ="Subject",all=F)
Cov <- merge(Cov,pc_cell,by ="Subject",all=F)


# Save column names of covariates for EWAS with cortical measures
CovName <- colnames(Cov)[-c(1)] 

# Generate files with complete data across all methylation, cortical measures, and covariates with matched individuals
Data <- merge(Cov,CorticalMeasure, by="Subject", all=F) 
Match <- match(Methy$Subject,Data$Subject) 
Cov <- Data[Match[!is.na(Match)],2:(length(CovName)+1)]  
CorticalMeasure <- Data[Match[!is.na(Match)],-c(1:(length(CovName)+1))] 
Methy <- Methy[!is.na(Match),-c(1)]

# Generate a dataframe to save age information for further analysis
AgeInfo <- data.frame(min(Cov$Age),max(Cov$Age),mean(Cov$Age))
names(AgeInfo) <- c("MinAge","MaxAge","MeanAge")

# EWAS analysis of cortical measures
# Get the numbers of methylation probes, cortical measures and covariates
Num_Methy <- ncol(Methy) 
Num_Cov <- ncol(Cov)
Num_Cortical <- ncol(CorticalMeasure)

###########################################################################################################
###################Association analysis: the following section should be run by ALL cohorts ###############
###########################################################################################################
#### ANALYSIS 1: EWAS analysis of cortical measures: in the whole sample, controlling for ICV ####
# Prepare the output files for the cortical EWAS
Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)

Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)

Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)

Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

# Linear regression
Covar <- matrix(unlist(Cov), ncol=ncol(Cov), byrow=F)
for (i in 1:Num_Methy) {
  for (j in 1:Num_Cortical) {
    Out <- summary(lm(Methy[,i] ~ Covar + CorticalMeasure[,j]))
    Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
    Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
    Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
    Origin_N[i,j] <- Out$df[1]+Out$df[2]
  }
}

# Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov),collapse = "+"),"+CorticalMeasure[,j]")

# Save the output files
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures.RData",sep=""),compress=T)


###################################################################################
#######################Alternative approach with parallelization###################
###################################################################################
# Note: only works on linux systems
# EWAS analysis of cortical measures: in the whole sample, controlling for ICV
# Prepare the output files for the cortical EWAS
Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)

Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)

Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)

Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

numberOfProcesses <- 5

doOne <- function(i, j)
{
  tryCatch({
    Out <- summary(lm(Methy[,i] ~ Covar + CorticalMeasure[,j]))
    return(as.data.frame(t(c(Out$coefficients[nrow(Out$coefficients),1], Out$coefficients[nrow(Out$coefficients),2], Out$coefficients[nrow(Out$coefficients),4], Out$df[1]+Out$df[2], names(Methy)[i]))))
  }, error = function(e) {
    cat("Warning: LM failed for probe ", names(Methy)[i], " and cortical region ", names(CorticalMeasure)[j], "\n")
    print(e)
    return(t(c(NA, NA, NA, NA, NA)))
  }, finally={})
}

library(pbmcapply)
Covar <- matrix(unlist(Cov), ncol=ncol(Cov), byrow=F)
for (j in 1:Num_Cortical) {
  print(paste0("Running brain region (", j, "/", Num_Cortical, "): ", names(CorticalMeasure)[j]))
  res = pbmclapply (1:Num_Methy, doOne, j=j, mc.cores = numberOfProcesses)
  res = dplyr::bind_rows(res)
  rownames(res) = res[,5]
  res = res[colnames(Methy),]
  Origin_Beta[,j] = res[,1]
  Origin_SE[,j] = res[,2]
  Origin_P[,j] = res[,3]
  Origin_N[,j] = res[,4]
  rm(res); 
}

# Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov),collapse = "+"),"+CorticalMeasure[,j]")

# Save the output files
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures.RData",sep=""),compress=T)

##############################################################################################################
##############################################################################################################
#### ANALYSIS 2: EWAS analysis of cortical measures: in the whole sample, not controlling for ICV ####
# Prepare the output files for the cortical EWAS
Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)

Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)

Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)

Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

# Linear regression
Cov_noICV <- Cov[, !(names(Cov) %in% c("ICV"))]
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)

for (i in 1:Num_Methy){
  for (j in 1:Num_Cortical){
    Out <- summary(lm(Methy[,i]~Covar + CorticalMeasure[,j]))
    Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
    Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
    Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
    Origin_N[i,j] <- Out$df[1]+Out$df[2]
  }
}

# Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov_noICV),collapse = "+"),"+CorticalMeasure[,j]")

# Save the output file of EWAS with cortical measures, in the whole sample
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)


###################################################################################
#######################Alternative approach with parallelization###################
###################################################################################
# Note: only works on linux systems
# EWAS analysis of cortical measures: in the whole sample, not controlling for ICV
# Prepare the output files for the cortical EWAS
Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_Beta) <- colnames(CorticalMeasure)
rownames(Origin_Beta) <- colnames(Methy)

Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_SE) <- colnames(CorticalMeasure)
rownames(Origin_SE) <- colnames(Methy)

Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_P) <- colnames(CorticalMeasure)
rownames(Origin_P) <- colnames(Methy)

Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
colnames(Origin_N) <- colnames(CorticalMeasure)
rownames(Origin_N) <- colnames(Methy)

numberOfProcesses <- 5

doOne <- function(i, j)
{
  tryCatch({
    Out <- summary(lm(Methy[,i] ~ Covar + CorticalMeasure[,j]))
    return(as.data.frame(t(c(Out$coefficients[nrow(Out$coefficients),1], Out$coefficients[nrow(Out$coefficients),2], Out$coefficients[nrow(Out$coefficients),4], Out$df[1]+Out$df[2], names(Methy)[i]))))
  }, error = function(e) {
    cat("Warning: LM failed for probe ", names(Methy)[i], " and cortical region ", names(CorticalMeasure)[j], "\n")
    print(e)
    return(t(c(NA, NA, NA, NA, NA)))
  }, finally={})
}

library(pbmcapply)
Cov_noICV <- Cov[, !(names(Cov) %in% c("ICV"))]
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)
for (j in 1:Num_Cortical) {
  print(paste0("Running brain region (", j, "/", Num_Cortical, "): ", names(CorticalMeasure)[j]))
  res = pbmclapply (1:Num_Methy, doOne, j=j, mc.cores = numberOfProcesses)
  res = dplyr::bind_rows(res)
  rownames(res) = res[,5]
  res = res[colnames(Methy),]
  Origin_Beta[,j] = res[,1]
  Origin_SE[,j] = res[,2]
  Origin_P[,j] = res[,3]
  Origin_N[,j] = res[,4]
  rm(res); 
}

# Save model description
model <- paste0("Methy[,i] ~ ",paste0(names(Cov_noICV),collapse = "+"),"+CorticalMeasure[,j]")

# Save the output file of EWAS with cortical measures, in the whole sample
save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,AgeInfo,file=paste("./Output_of_",cohort,"_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)

##############################################################################################################
##############################################################################################################
#### ANALYSIS 3: EWAS analysis of cortical measures: stratified by sex ####
# Please make sure Sex was specified as follows: Males=1, Females=2
Num_Methy <- ncol(Methy) 
Num_Cov <- ncol(Cov_noICV)
Num_Cortical <- ncol(CorticalMeasure)

Gender <- colnames(Cov_noICV)=="Sex"
Category <- c("Male","Female")

# Save age information for each sex
MaleCov <- subset(Cov_noICV, Sex==1)
MaleAgeInfo <- data.frame(min(MaleCov$Age),max(MaleCov$Age),mean(MaleCov$Age))
names(MaleAgeInfo) <- c("MaleMinAge","MaleMaxAge","MaleMeanAge")
FemaleCov <- subset(Cov_noICV, Sex==2)
FemaleAgeInfo <- data.frame(min(FemaleCov$Age),max(FemaleCov$Age),mean(FemaleCov$Age))
names(FemaleAgeInfo) <- c("FemaleMinAge","FemaleMaxAge","FemaleMeanAge")

# Linear regression
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)
for (k in 1:2) {
  Ind <- Covar[,Gender]==(k) 
  Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_Beta) <- colnames(CorticalMeasure)
  rownames(Origin_Beta) <- colnames(Methy)
  Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_SE) <- colnames(CorticalMeasure)
  rownames(Origin_SE) <- colnames(Methy)
  Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_P) <- colnames(CorticalMeasure)
  rownames(Origin_P) <- colnames(Methy)
  Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_N) <- colnames(CorticalMeasure)
  rownames(Origin_N) <- colnames(Methy) 
  for (i in 1:Num_Methy) {
    for (j in 1:Num_Cortical) {
      # Remove gender from covariates
      Out <- summary(lm(Methy[Ind,i]~Covar[Ind,!Gender]+CorticalMeasure[Ind,j]))         
      Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
      Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
      Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
      Origin_N[i,j] <- Out$df[1]+Out$df[2]
    }
  }
  # Save model description
  model <- paste0("Methy[Ind,i] ~ ",paste0(names(Cov_noICV)[!Gender],collapse = "+"),"+CorticalMeasure[Ind,j]")
  
  # Save output files for each gender separately
  save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,MaleAgeInfo,FemaleAgeInfo,file=paste("./Output_of_",cohort,"_",Category[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
}

###################################################################################
#######################Alternative approach with parallelization###################
###################################################################################
  doOne <- function(i, j)
  {
    tryCatch({
      Out <- summary(lm(Methy[Ind,i]~Covar[Ind,!Gender]+CorticalMeasure[Ind,j]))
      return(as.data.frame(t(c(Out$coefficients[nrow(Out$coefficients),1], Out$coefficients[nrow(Out$coefficients),2], Out$coefficients[nrow(Out$coefficients),4], Out$df[1]+Out$df[2], names(Methy)[i]))))
    }, error = function(e) {
      cat("Warning: LM failed for gender ", k, " and probe ", names(Methy)[i], " and cortical region ", names(CorticalMeasure)[j], "\n")
      print(e)
      return(t(c(NA, NA, NA, NA, NA)))
    }, finally={})
  }

Gender <- colnames(Cov_noICV)=="Sex"
Category <- c("Male","Female")

Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)

for (k in 1:2) {
    Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_Beta) <- colnames(CorticalMeasure)
    rownames(Origin_Beta) <- colnames(Methy)
    Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_SE) <- colnames(CorticalMeasure)
    rownames(Origin_SE) <- colnames(Methy)
    Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_P) <- colnames(CorticalMeasure)
    rownames(Origin_P) <- colnames(Methy)
    Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_N) <- colnames(CorticalMeasure)
    rownames(Origin_N) <- colnames(Methy)
    
    Ind <- Covar[,Gender]==(k)
    for(j in 1:Num_Cortical){
      print(paste0("Running sex: ", Category[k], " brain region (", j, "/", Num_Cortical, "): ", names(CorticalMeasure)[j]))
      res = pbmclapply(1:Num_Methy, doOne, j=j, mc.cores = numberOfProcesses, mc.style = "ETA")
      res = dplyr::bind_rows(res)
      rownames(res) = res[,5]
      res = res[colnames(Methy),]
      Origin_Beta[,j] = res[,1]
      Origin_SE[,j] = res[,2]
      Origin_P[,j] = res[,3]
      Origin_N[,j] = res[,4]
      rm(res)
    }
   
    # Save model description
    model <- paste0("Methy[Ind,i] ~ ",paste0(names(Cov_noICV)[!Gender],collapse = "+"),"+CorticalMeasure[Ind,j]")
    
    # Save output file for each gender separately
    save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,MaleAgeInfo,FemaleAgeInfo,file=paste("./Output_of_",cohort,"_",Category[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
  }
  
  


###########################################################################################################
##########Cohorts with patients should also run the following section, in addition to the above ###########
###########################################################################################################
# Make sure the AffectionStatus is coded as 1 (Case) and 0 (Control)
CaseCov <- subset(Cov, AffectionStatus==1)
CaseAgeInfo <- data.frame(min(CaseCov$Age),max(CaseCov$Age),mean(CaseCov$Age))
names(CaseAgeInfo) <- c("CaseMinAge","CaseMaxAge","CaseMeanAge")
ControlCov <- subset(Cov, AffectionStatus==0)
ControlAgeInfo <- data.frame(min(ControlCov$Age),max(ControlCov$Age),mean(ControlCov$Age))
names(ControlAgeInfo) <- c("ControlMinAge","ControlMaxAge","ControlMeanAge")
  
# Check if there is a column named "AffectionStatus", i.e. a case-control study
Affect <- colnames(Cov_noICV)=="AffectionStatus"
Status <- c("Case","Control")

Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)

if (sum(Affect)==1) {  
  for (k in 1:2) {
    #AffectionStatus should be coded as 1 (Case) and 0 (Control)
    Ind <- Covar[,Affect]==(2-k) 
    Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_Beta) <- colnames(CorticalMeasure)
    rownames(Origin_Beta) <- colnames(Methy)
    Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_SE) <- colnames(CorticalMeasure)
    rownames(Origin_SE) <- colnames(Methy)
    Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_P) <- colnames(CorticalMeasure)
    rownames(Origin_P) <- colnames(Methy)
    Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
    colnames(Origin_N) <- colnames(CorticalMeasure)
    rownames(Origin_N) <- colnames(Methy) 
    for (i in 1:Num_Methy) {
      for (j in 1:Num_Cortical) {
        # Remove affection status from covariates
        Out <- summary(lm(Methy[Ind,i]~Covar[Ind,!Affect]+CorticalMeasure[Ind,j]))         
        Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
        Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),2]
        Origin_P[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
        Origin_N[i,j] <- Out$df[1]+Out$df[2]
        }
    }
    # Save model description
    model <- paste0("Methy[Ind,i] ~ ",paste0(names(Cov_noICV)[!Affect],collapse = "+"),"+CorticalMeasure[Ind,j]")
    
    # Save output file, for cases and controls separately    
    save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,CaseAgeInfo,ControlAgeInfo,file=paste("./Output_of_",cohort,"_",Status[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
  }
}
      
###################################################################################
#######################Alternative approach with parallelization###################
###################################################################################
doOne <- function(i, j)
  {
  tryCatch({
    Out <- summary(lm(Methy[Ind,i]~Covar[Ind,!Affect]+CorticalMeasure[Ind,j]))
    return(as.data.frame(t(c(Out$coefficients[nrow(Out$coefficients),1], Out$coefficients[nrow(Out$coefficients),2], Out$coefficients[nrow(Out$coefficients),4], Out$df[1]+Out$df[2], names(Methy)[i]))))
    }, error = function(e) {
      cat("Warning: LM failed for affection status ", k, " and probe ", names(Methy)[i], " and cortical region ", names(CorticalMeasure)[j], "\n")
      print(e)
      return(t(c(NA, NA, NA, NA, NA)))
    }, finally={})
  }

Affect <- colnames(Cov_noICV)=="AffectionStatus"
Status <- c("Case","Control")
      
Covar <- matrix(unlist(Cov_noICV), ncol=ncol(Cov_noICV), byrow=F)
      
for (k in 1:2) {
  Origin_Beta <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_Beta) <- colnames(CorticalMeasure)
  rownames(Origin_Beta) <- colnames(Methy)
  Origin_SE <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_SE) <- colnames(CorticalMeasure)
  rownames(Origin_SE) <- colnames(Methy)
  Origin_P <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_P) <- colnames(CorticalMeasure)
  rownames(Origin_P) <- colnames(Methy)
  Origin_N <- matrix(data=NA,nrow=Num_Methy,ncol=Num_Cortical,byrow=F,dimnames=NULL)
  colnames(Origin_N) <- colnames(CorticalMeasure)
  rownames(Origin_N) <- colnames(Methy)
  Ind <- Covar[,Affect]==(2-k)
  
  for(j in 1:Num_Cortical){
    print(paste0("Running affection status: ", Status[k], " brain region (", j, "/", Num_Cortical, "): ", names(CorticalMeasure)[j]))
    res = pbmclapply(1:Num_Methy, doOne, j=j, mc.cores = numberOfProcesses, mc.style = "ETA")
    res = dplyr::bind_rows(res)
    rownames(res) = res[,5]
    res = res[colnames(Methy),]
    Origin_Beta[,j] = res[,1]
    Origin_SE[,j] = res[,2]
    Origin_P[,j] = res[,3]
    Origin_N[,j] = res[,4]
    rm(res)
    }
  # Save model description
  model <- paste0("Methy[Ind,i] ~ ",paste0(names(Cov_noICV)[!Affect],collapse = "+"),"+CorticalMeasure[Ind,j]")
  
  # Save output files, for cases and controls separately  
  save(Origin_Beta,Origin_SE,Origin_P,model,Origin_N,ListInvarProbe,CaseAgeInfo,ControlAgeInfo,file=paste("./Output_of_",cohort,"_",Status[k],"_Individuals_EWAS_with_CorticalMeasures_noICV.RData",sep=""),compress=T)
}

    