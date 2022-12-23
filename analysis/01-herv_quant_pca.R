################################################################################
################################################################################
################################################################################
##############################  RiboSeq Bulk / SC PCA ##########################

#################################### SETUP #####################################

library(tidyverse)
library(tibble)
library(gridExtra)
library(cowplot)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(PCAtools)
library(DESeq2)
library(scopetools)
library(tximport)
library(tximportData)


################################## LOAD DATA ###################################

# Read in sample metadata
samples <- read.csv("/fsx/projects/te_riboseq_atlas/config/samples_riboseq.csv", 
                    sep = ",")

samples$tissue_cell_type <- paste0(samples$Tissue, samples$Cell_type)

# sample names
sample_names <- samples$BioSample

# set rownames 
row.names(samples) <- samples$BioSample

# Load files
t_files <- file.path("/fsx/projects/te_riboseq_atlas/results/salmon/all_salmon", 
                     paste0(sample_names, 
                            '_',
                            'quant.sf'))
names(t_files) <- samples$bulk_RNAseq

counts <- scopetools::load_tsv_counts(t_files, 
                                      colnames = names(t_files), 
                                      gene_id_column = 1, 
                                      count_column = 5,
                                      header=TRUE)


samples$tissue_cell_type <- paste0(samples$Tissue, samples$Cell_type)
samples$cell_treatment <- paste0(samples$Cell_type, samples$TREATMENT)
samples_fibroblasts <- samples[samples$tissue_cell_type == "Fibroblast",]
samples_other <- samples[samples$TREATMENT == "untreated",]


counts_fibroblasts <- counts[,rownames(samples_fibroblasts)]
counts_other <- counts[,rownames(samples_other)]


################################### DESEQ2 #####################################

# All
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts),
                                      colData = samples,
                                      design = ~ 1)

dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)


# Fibroblasts
dds_fibro <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts_fibroblasts),
                                            colData = samples_fibroblasts,
                                            design = ~ 1)

dds_fibro <- DESeq2::DESeq(dds_fibro, parallel=T)
tform_fibro <- DESeq2::varianceStabilizingTransformation(dds_fibro, blind=FALSE)

# Other
dds_other <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts_other),
                                            colData = samples_other,
                                            design = ~ 1)

dds_other <- DESeq2::DESeq(dds_other, parallel=T)
tform_other <- DESeq2::varianceStabilizingTransformation(dds_other, blind=FALSE)

#################################### PCA  ######################################

removeVar <- 0.1
pca.obj <- 
  PCAtools::pca(assay(tform), metadata=samples, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
            removeVar*100, length(pca.obj$xvars)))

varline <- 50
varline.x <- min(which(cumsum(pca.obj$variance) >= varline))

horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
elbow <- PCAtools::findElbowPoint(pca.obj$variance)

PCAtools::screeplot(pca.obj,
                    axisLabSize = 6,
                    components = getComponents(pca.obj, 1:30),
                    hline=varline, vline=c(varline.x, horn$n, elbow)
) +
  geom_label(aes(x=varline.x+1, y=50, 
                 label = paste0(varline, '% var'), vjust = -1)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1))


cat(sprintf('%d PCs for Elbow method\n', elbow)) # 3
cat(sprintf('%d PCs for Horn method\n', horn$n)) # 4
cat(sprintf('%d PCs needed to explain %d percent of variation\n', varline.x, varline)) # 1

################################### PCA ALL ####################################

biplot(pca.obj, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       shape = "TREATMENT",
       colby = "tissue_cell_type",
       colkey = c("HUVEC" = "#E64B35B2", "Brain" = "#4DBBD5B2",
                  "ES" = "#00A087B2", "fat" = "#3C5488B2",
                  "HA_EC" = "#F39B7FB2", "HCAEC" = "#8491B4B2",
                  "Hepatocytes" = "#91D1C2B2", "VSMC" = "#DC0000B2", 
                  "Fibroblast" = "grey"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       colLegendTitle = "Tissue / Cell Type",
       shapeLegendTitle = "Treatment")

biplot(pca.obj, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse=TRUE,
       ellipseLevel = 0.95,
       shape = "TREATMENT",
       colby = "tissue_cell_type",
       colkey = c("HUVEC" = "#E64B35B2", "Brain" = "#4DBBD5B2",
                  "ES" = "#00A087B2", "fat" = "#3C5488B2",
                  "HA_EC" = "#F39B7FB2", "HCAEC" = "#8491B4B2",
                  "Hepatocytes" = "#91D1C2B2", "VSMC" = "#DC0000B2", 
                  "Fibroblast" = "grey"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       colLegendTitle = "Tissue / Cell Type",
       shapeLegendTitle = "Treatment")

biplot(pca.obj, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse=TRUE,
       ellipseLevel = 0.95,
       xlim = c(-300, 450),
       ylim = c(-400, 400),
       shape = "TREATMENT",
       colby = "tissue_cell_type",
       colkey = c("HUVEC" = "#E64B35B2", "Brain" = "#4DBBD5B2",
                  "ES" = "#00A087B2", "fat" = "#3C5488B2",
                  "HA_EC" = "#F39B7FB2", "HCAEC" = "#8491B4B2",
                  "Hepatocytes" = "#91D1C2B2", "VSMC" = "#DC0000B2", 
                  "Fibroblast" = "grey"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       colLegendTitle = "Tissue / Cell Type",
       shapeLegendTitle = "Treatment")

################################## PCA FIBRO ###################################

pca.obj.fibro <- 
  PCAtools::pca(assay(tform_fibro), metadata=samples_fibroblasts, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
            removeVar*100, length(pca.obj$xvars)))

biplot(pca.obj.fibro, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse=TRUE,
       ellipseLevel = 0.95,
       shape = "source_name",
       colby = "source_name",
       xlim = c(-200, 400),
       ylim = c(-150, 150),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       colLegendTitle = "Treatment",
       shapeLegendTitle = "Treatment")

biplot(pca.obj.fibro, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       shape = "source_name",
       colby = "source_name",
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       xlim = c(-200, 400),
       ylim = c(-150, 150),
       colLegendTitle = "Treatment",
       shapeLegendTitle = "Treatment")

################################## PCA OTHER ###################################

pca.obj.other <- 
  PCAtools::pca(assay(tform_other), metadata=samples_other, removeVar=removeVar)
cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
            removeVar*100, length(pca.obj$xvars)))

biplot(pca.obj.other, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3,
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "tissue_cell_type",
       colkey = c("HUVEC" = "#E64B35B2", "Brain" = "#4DBBD5B2",
                  "ES" = "#00A087B2", "fat" = "#3C5488B2",
                  "HA_EC" = "#F39B7FB2", "HCAEC" = "#8491B4B2",
                  "Hepatocytes" = "#91D1C2B2", "VSMC" = "#DC0000B2",
                  "Fibroblast" = "grey"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       xlim = c(-200, 300),
       ylim = c(-150, 250),
       legendPosition = "right",
       colLegendTitle = "Tissue / Cell Type")

biplot(pca.obj.other, 
       x="PC1",
       y="PC2",
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3,
       encircle = FALSE,
       sizeLoadingsNames = 3,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse=TRUE,
       ellipseType = "t",
       ellipseLevel = 0.95,
       xlim = c(-300, 400),
       ylim = c(-350, 400),
       colby = "tissue_cell_type",
       colkey = c("HUVEC" = "#E64B35B2", "Brain" = "#4DBBD5B2",
                  "ES" = "#00A087B2", "fat" = "#3C5488B2",
                  "HA_EC" = "#F39B7FB2", "HCAEC" = "#8491B4B2",
                  "Hepatocytes" = "#91D1C2B2", "VSMC" = "#DC0000B2",
                  "Fibroblast" = "grey"),
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       legendPosition = "right",
       colLegendTitle = "Tissue / Cell Type") 
