#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: 01-Clean-metadata.R
# Description: Clean and tidy the metadata for ENA: ERP106487, 
#E-MTAB-5464: RNA sequencing of purified intestinal epithelial cells from
# paediatric biopsies including Inflammatory Bowel Disease and 
# healthy control (Homo spaiens).
#
# Author: Mike Martinez
# Project: Howell Analysis
# Env: /Users/mike/Desktop/HowellData/envs/pixi.toml
# Date created: January 27th, 2026
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Libraries
library(dplyr)
library(tidyverse)

#----- Set project directory
wd <- "/Users/mike/Desktop/HowellData/"

#----- Specify directories
inputDir <- paste0(wd, "data/")
outputRoot <- paste0(wd, "outputs/")
outputDir <- paste0(outputRoot, "01-Howell-EDA/")
figDir <- paste0(outputDir, "figures/")

#----- Create directories
dirs <- c(outputRoot, outputDir, figDir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE COUNTS AND METADATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Counts
counts <- read.csv(paste0(inputDir, "Howell-RawCounts-matrix.tsv"), sep = "\t")
counts$Gene <- paste(counts$Gene.ID, counts$Gene.Name, sep = " - ")
rownames(counts) <- counts$Gene
counts$Gene.ID <- NULL
counts$Gene.Name <- NULL
counts$Gene <- NULL

#----- Read in metadata
meta <- read.csv(paste0(inputDir, "Howell-metadata-matrix.tsv"), sep = "\t")
rownames(meta) <- meta$Run
meta <- meta[,c(1, 2, 4, 6, 8, 10, 14, 16, 18, 20)]
meta <- meta[,c(1:8)]
colnames(meta) <- c("Sample", "Age", "Cell_Type", "Inflammation", "Disease", "Patient_Number", "Region", "Sex")

#----- Arrange by patient number
meta <- meta %>%
  arrange(Patient_Number)

#----- Clean up metadata
meta$Cell_Type <- recode(meta$Cell_Type, "intestinal epithelial cell" = "IECs")
meta$Inflammation <- recode(meta$Inflammation, 
                                "no inflammation" = "No",
                                "inflammation" = "Yes")
meta$Disease <- recode(meta$Disease, 
                       "Crohn's disease" = "CD",
                       "ulcerative colitis" = "UC",
                       "normal" = "Normal")
meta$Region <- recode(meta$Region, 
                      "ascending colon" = "Ascending_Colon",
                      "sigmoid colon" = "Sigmoid_Colon",
                      "terminal ileum" = "Terminal_Ileum")
meta$Sex <- recode(meta$Sex,
                   "male" = "M",
                   "female" = "F")
rownames(meta) <- meta$Sample

#----- Re-order the counts based on metadata
counts <- counts[,rownames(meta)]
all(colnames(counts) == rownames(meta))

#----- Save
write.csv(meta, file = paste0(outputDir, "Cleaned-metadata.csv"))
write.csv(counts, file = paste0(outputDir, "Cleaned-counts.csv"))



