#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: 02-EDA.R
# Description: Perform exploratory data analysis on the 
# Howell data set.
#
# Author: Mike Martinez
# Project: Howell Data
# Env: /Users/mike/Desktop/HowellData/envs/pixi.toml
# Date created: Tuesday January 27th. 2026
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Libraries
library(sessioninfo)
library(dplyr)
library(DESeq2)
library(RGenEDA)
library(viridis)
library(ggplot2)

#----- Set project directory
wd <- "/Users/mike/Desktop/HowellData/"

#----- Specify directories
inputDir <- paste0(wd, "outputs/01-Howell-EDA/")
outputDir <- paste0(wd, "outputs/02-EDA/")
figDir <- paste0(outputDir, "figures/")

#----- Create directories
dirs <- c(outputDir, figDir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE PROCESSED COUNTS AND METADATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Counts
counts <- read.csv(paste0(inputDir, "Cleaned-counts.csv"))
rownames(counts) <- counts$X
counts$X <- NULL

#----- Metadata
meta <- read.csv(paste0(inputDir, "Cleaned-metadata.csv"))
rownames(meta) <- meta$X

#----- Sanity check
all(colnames(counts) == rownames(meta))

#----- Convert numeric metadata to factor
meta$Patient_Number <- factor(meta$Patient_Number)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE DESEQ2 OBJECT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Deseq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~Disease)

#----- Set reference levels
dds$Disease <- relevel(dds$Disease, ref = "Normal")

#----- Run DESeq2
dds <- DESeq(dds)

#----- Rlog transform the data and extract normalized matrix
rld <- rlog(dds)
mat <- assay(rld)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE GENEDA OBJECT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Create a color vector for age
age_levels <- c("6 year", "7 year", "8 year", "9 year",
                "10 year", "11 year", "12 year", "13 year",
                "14 year", "15 year")
meta$Age <- factor(meta$Age, levels = age_levels)
ageCols <- setNames(
  viridis(length(age_levels), option = "C", begin = 0.15, end = 0.9),
  age_levels
)

#----- Get number of patients
num_patients <- length(unique(meta$Patient_Number))

#----- Create color vector list for RGenEDA
colors <- list(
  Cell_Type    = c("IECs" = "#4C72B0"),
  Age          = ageCols,
  Inflammation = c("No"  = "#9E9E9E",
                   "Yes" = "#C44E52"),
  Disease      = c("Normal" = "#0072B2",
                   "UC"     = "#E69F00",
                   "CD"     = "#009E73"),
  Region       = c("Ascending_Colon" = "#4C72B0",
                   "Sigmoid_Colon"   = "#C44E52",
                   "Terminal_Ileum"  = "#55A868"),
  Sex          = c("M" = "#4C72B0",
                   "F" = "#DD8452"),
  Patient_Number = setNames(rainbow(num_patients), sort(unique(meta$Patient_Number)))
)

#----- Initialize GenEDA object with normalized counts and metadata
obj <- GenEDA(
  normalized = mat,
  metadata = meta)

#----- Plot normalized count distributions
PlotCountDist(obj, split_by = "Disease")

#----- Plot Eucliden distances
hm <- PlotDistances(
  obj,
  meta_cols = c("Sex", "Region", 
                "Inflammation", "Disease"),
  palettes = colors,
  return = "plot"
)
hm$heatmap
GenSave(hm, paste0(figDir, "EucDists.png"), width = 12, height = 10)


#----- Find HVGs
PlotHVGVariance(obj)
obj <- FindVariableFeatures(obj, 3000)
ggsave(paste0(figDir, "HVG_Elbowplot.png"), width = 6, height = 6)

#----- Run PCA
obj <- RunPCA(obj)
PlotScree(obj)
ggsave(paste0(figDir, "Screeplot.png"), width = 6, height = 6)
#PC1       PC2       PC3       PC4       PC5 
#"35.59 %" "19.54 %"  "5.99 %"  "4.41 %"  "3.97 %" 

#----- Plot eigencorrelations
ec <- PlotEigenCorr(obj, 
                    meta_cols = c("Sex", "Region", "Patient_Number", 
                                  "Inflammation", "Disease", "Age"),
                    num_pcs = 5)
ec$plot
ggsave(paste0(figDir, "Eigencorr.png"), width = 6, height = 6)

#----- Explore PCA in more depth
pcaDF <- ExtractPCA(obj)

ggplot(pcaDF, aes(x = PC1, y = PC2, fill = Region)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold")) +
  labs(x = "PC1: 33.59%", y = "PC2: 19.54%") +
  scale_fill_manual(
    values = colors[["Region"]],
    labels = c("Ascending Colon", "Sigmoid Colon", "Terminal Ileum")) +
  facet_grid(~Inflammation)
ggsave(paste0(figDir, "PCA.png"), width = 8, height = 8)




