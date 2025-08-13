# Meta-analysis of Brain EWAS

This repository contains code for conducting a **meta-analysis epigenome-wide association study (EWAS)** investigating DNA methylation (DNAm) signatures associated with **cortical morphology**, focusing on two key measures:

- **Cortical thickness (CT)**
- **Cortical surface area (SA)**

The workflow is developed as part of the **ENIGMA-Epigenetics Working Group**, enabling large-scale, harmonized analyses across international cohorts.

## Overview

The scripts are designed to:

- Perform cohort-level quality control (QC) on DNAm and cortical measures, including aggregation of global and regional CT and SA metrics  
- Conduct cohort-level EWAS on cortical measures  
- Harmonise and quality-check EWAS results across multiple cohorts (QCEWAS)  
- Run fixed-effects meta-analyses using **METAL** on discovery cohorts, with replication testing in independent cohorts and combined datasets  
- Perform functional and genomic enrichment analyses, as well as phenome-wide association studies (**PheWAS**)   
- Conduct **Mendelian randomisation** analyses to investigate potential causal pathways from meQTLs → CT- or SA-associated DNAm → neuropsychiatric disorders or educational attainment  
