#### Use bacon to remove inflation and bias for each cohort and save corrected data for meta-analysis ####
library(bacon)
library(BiocParallel)
register(MulticoreParam(1, log=TRUE))

# Define a function to perform BACON correction 
processBACON <- function(dataframe,dataframe_name) {
  # Extract columns
  es <- dataframe[, c(2:19)]
  se <- dataframe[, c(20:37)]
  pval <- dataframe[, c(38:55)]
  n <- dataframe[, c(56:73)]
  
  # Perform bacon correction
  bc <- bacon(NULL, es, se, na.exclude = TRUE)
  
  # Extract bacon-corrected values
  es_corrected <- es(bc)
  se_corrected <- se(bc)
  pval_corrected <- pval(bc)
  colnames(pval_corrected) <- colnames(pval)
  
  # Additional bacon outputs
  inflation <- inflation(bc)
  bias <- bias(bc)
  
  # Merge corrected data
  Bacon_corrected <- data.frame(dataframe[1], es_corrected, se_corrected, pval_corrected, n)
  
  # Identify output file path
  path <- "for meta-analysis/cortical sumstats/data/Bacon_corrected_data_noICV/"
  
  inflation_file_name <- paste0(path, dataframe_name, "_noICV_inflation.txt")
  bias_file_name <- paste0(path, dataframe_name, "_noICV_bias.txt")
  corrected_sumstats_file_name <- paste0(path, dataframe_name, "_noICV_Bacon_corrected_sumstats.txt")
  
  # Save outputs to files
  capture.output(print(inflation), file = inflation_file_name)
  capture.output(print(bias), file = bias_file_name)
  write.table(Bacon_corrected, file = corrected_sumstats_file_name, row.names = FALSE, quote = FALSE)
}

# Define function to remove cross-reactive probes for 450k or EPIC arrays
processData <- function(is_epic_array = FALSE) {
  # Load the cross-reactive probes
  if (is_epic_array) {
    cross_reactive <- read.csv("for meta-analysis/Cross-reactive probes of EPIC array.csv")
    probe_id_column <- "ProbeID"
  } else {
    cross_reactive <- read.csv("for meta-analysis/CrossReactive_Probes.csv")
    probe_id_column <- "TargetID"
  }
  
  # Assuming Origin_Beta, Origin_SE, Origin_P, and Origin_N are loaded in the environment
  Origin_Beta <- Origin_Beta[!(rownames(Origin_Beta) %in% cross_reactive[[probe_id_column]]),]
  Origin_SE <- Origin_SE[!(rownames(Origin_SE) %in% cross_reactive[[probe_id_column]]),]
  Origin_P <- Origin_P[!(rownames(Origin_P) %in% cross_reactive[[probe_id_column]]),]
  Origin_N <- Origin_N[!(rownames(Origin_N) %in% cross_reactive[[probe_id_column]]),]
  
  # Merge all data together
  All <- data.frame(cbind(rownames(Origin_Beta), Origin_Beta, Origin_SE, Origin_P, Origin_N))
  col_names <- c("CpG", sapply(1:72, function(i) paste(c("MeanFrontalThickness", "MeanTemporalThickness", "MeanOccipitalThickness",
                                                         "MeanParietalThickness", "MeanCingulateThickness", "MeanInsulaThickness",
                                                         "LeftHemisphereThickness", "RightHemisphereThickness", "MeanThickness",
                                                         "TotalFrontalArea", "TotalTemporalArea", "TotalOccipitalArea", "TotalParietalArea",
                                                         "TotalCingulateArea", "TotalInsulaArea", "LeftHemisphereArea", "RightHemisphereArea",
                                                         "TotalArea"), rep(c("beta", "SE", "P", "N"), each = 18), sep = "_")[i]))
  colnames(All) <- col_names
  
  # Convert all characters to numeric
  for (i in 2:73) {
    All[, i] <- as.numeric(All[, i])
  }
  
  return(All)
}

# for 450k cohorts
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/ENIGMA-Epigenetics2/BIPSYD/Output_of_CTAAC_EWAS_with_CorticalMeasures_noICV.RData")
results1 <- processData(is_epic_array = FALSE)
processBACON(results1, "CTAAC")

# for EPIC cohorts
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/ENIGMA-Epigenetics2/BIPSYD/Output_of_IMAGEN_FU2_EWAS_with_CorticalMeasures_noICV.RData")
results1 <- processData(is_epic_array = TRUE)
processBACON(results1, "IMAGEN_FU2")
