# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(bioplotr)

# Import data
setwd('./Documents/QMUL/PCC/AlexDaCosta')
clin <- fread('./Data/clinical.csv')
t2g <- readRDS('./Data/Hs91.t2g.rds')
files <- file.path('./Data/Counts', clin$Idx, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Status)  

# Filter, transform counts
dds <- estimateSizeFactors(dds)
keep <- rowSums(cpm(counts(dds)) > 1) >= 3  # Or 12?
mat <- counts(dds, normalized = TRUE)[keep, ]
mat <- cpm(mat, log = TRUE, prior.count = 1)
colnames(mat) <- clin$Sample

# EDA
plot_mv(mat, trans = 'sqrt')
plot_drivers(mat, clin[, .(Sample, Mutation, Status)], n.pc = 5, alpha = 0.05)
lcpm <- expression(log[2]*'-CPM')
plot_density(mat, group = list(Status = clin$Status), xlab = lcpm) 
plot_box(mat, group = list(Status = clin$Status), ylab = lcpm)
plot_similarity(mat, group = list(Status = clin$Status))
plot_pca(mat, group = list(Status = clin$Status), label = TRUE)
plot_tsne(mat, group = list(Status = clin$Status), label = TRUE)

### SOM ###
library(kohonen)
library(limma)
library(doMC)
registerDoMC(4)

# Build SOM
des <- model.matrix(~ Status, data = clin)
y <- plot_som(mat, top = 10000, design = des, coef = 2, stat = 't',
              title = 'SOM: Healthy vs. Disease')


