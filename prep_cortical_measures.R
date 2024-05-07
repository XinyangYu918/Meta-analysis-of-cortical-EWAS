#!/usr/bin/env Rscript
#
# Prepare cortical measures (average thickness and total area) for ENIGMA EWAS
# based on four input files named (lh|rh).aparc.(thickness|area).tsv

library(tidyverse)


CorticalMeasure <- lapply(c("thickness", "area"), function(m) {
  
  # desired measures: AVERAGE thickness and TOTAL area
  rowFun <- function(x) { if (m == "thickness") { rowMeans(x) } else { rowSums(x) }}
  
  # calculate aggregate measures WITHIN left and right hemisphere, respectively
  agg1 <- lapply(c("lh", "rh"), function(h) {
    
    RawMRI <- read.csv(file = paste0(h, ".aparc.", m, ".csv"), header = TRUE)
    Agg <- data.frame(Subject=RawMRI[1])
    
    Agg$Frontal <- rowFun(RawMRI[,c(paste0(h, "_superiorfrontal_", m),
                                    paste0(h, "_rostralmiddlefrontal_", m),
                                    paste0(h, "_caudalmiddlefrontal_", m),
                                    paste0(h, "_parsopercularis_", m),
                                    paste0(h, "_parstriangularis_", m),
                                    paste0(h, "_parsorbitalis_", m),
                                    paste0(h, "_lateralorbitofrontal_", m),
                                    paste0(h, "_medialorbitofrontal_", m),
                                    paste0(h, "_precentral_", m),
                                    paste0(h, "_paracentral_", m),
                                    paste0(h, "_frontalpole_", m))])
    
    Agg$Temporal <- rowFun(RawMRI[,c(paste0(h, "_superiortemporal_", m),
                                     paste0(h, "_bankssts_", m),
                                     paste0(h, "_fusiform_", m),
                                     paste0(h, "_transversetemporal_", m),
                                     paste0(h, "_entorhinal_", m),
                                     paste0(h, "_temporalpole_", m),
                                     paste0(h, "_parahippocampal_", m))])
    
    Agg$Occipital <- rowFun(RawMRI[,c(paste0(h, "_lateraloccipital_", m),
                                      paste0(h, "_lingual_", m),
                                      paste0(h, "_cuneus_", m),
                                      paste0(h, "_pericalcarine_", m))])
    
    Agg$Parietal <- rowFun(RawMRI[,c(paste0(h, "_superiorparietal_", m),
                                     paste0(h, "_inferiorparietal_", m),
                                     paste0(h, "_supramarginal_", m),
                                     paste0(h, "_postcentral_", m),
                                     paste0(h, "_precuneus_", m))])
    
    Agg$Cingulate <- rowFun(RawMRI[,c(paste0(h, "_rostralanteriorcingulate_", m),
                                      paste0(h, "_caudalanteriorcingulate_", m),
                                      paste0(h, "_posteriorcingulate_", m),
                                      paste0(h, "_isthmuscingulate_", m))])
    
    Agg$Insula <- RawMRI[[paste0(h, "_insula_", m)]]
    
    # Left/Right prefix
    pref <- ifelse(h == "lh", "Left", "Right")
    Agg <- rename_with(Agg, .fn = ~paste0(pref, .x), .cols = -Subject)
    
    return(Agg)
  }) %>% Reduce(function(x, y) merge(x, y, by = "Subject"), .)
  
  # aggregate measures ACROSS left and right hemisphere
  Agg2 <- data.frame(Subject=agg1$Subject)
  Agg2$Frontal     <- rowFun(agg1[,c("LeftFrontal",   "RightFrontal")])
  Agg2$Temporal    <- rowFun(agg1[,c("LeftTemporal",  "RightTemporal")])
  Agg2$Occipital   <- rowFun(agg1[,c("LeftOccipital", "RightOccipital")])
  Agg2$Parietal    <- rowFun(agg1[,c("LeftParietal",  "RightParietal")])
  Agg2$Cingulate   <- rowFun(agg1[,c("LeftCingulate", "RightCingulate")])
  Agg2$Insula      <- rowFun(agg1[,c("LeftInsula",    "RightInsula")])
  
  # Mean/Total prefix
  pref <- ifelse(m == "thickness", "Mean", "Total")
  Agg2 <- rename_with(Agg2, .fn = ~paste0(pref, .x), .cols = -Subject)
  
  Agg2$LeftHemisphere  <- rowFun(agg1[,c("LeftFrontal",    "LeftTemporal",
                                         "LeftOccipital",  "LeftParietal",
                                         "LeftCingulate",  "LeftInsula")])
  Agg2$RightHemisphere <- rowFun(agg1[,c("RightFrontal",   "RightTemporal",
                                         "RightOccipital", "RightParietal",
                                         "RightCingulate", "RightInsula")])
  
  if (m == "thickness") {
    Agg2$Mean <- rowMeans(Agg2[,c("LeftHemisphere", "RightHemisphere")])
  } else {
    Agg2$Total <- rowSums(Agg2[,c("LeftHemisphere", "RightHemisphere")])
  }
  
  # Thickness/Area suffix
  suff <- ifelse(m == "thickness", "Thickness", "Area")
  Agg2 <- rename_with(Agg2, .fn = ~paste0(.x, suff), .cols = -Subject)
  
  return(Agg2)
}) %>% Reduce(function(x, y) merge(x, y, by = "Subject"), .)


# Save the merged cortical thickness and surface area measures for later use
write.csv(CorticalMeasure, file = "CorticalMeasure.csv", na = "", row.names = FALSE)
