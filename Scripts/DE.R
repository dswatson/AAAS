# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(bioplotr)

# Import data
setwd('./Documents/QMUL/AAAS')
clin <- fread('./Data/clinical.csv', stringsAsFactors = TRUE) %>%
  mutate(Status = relevel(Status, ref = 'Healthy'))
t2g <- readRDS('./Data/Hs91.t2g.rds')
files <- file.path('./Data/Counts', clin$Idx, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)

# Differential expression
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Status)  
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, betaPrior = TRUE)
res <- results(dds, coef = 2, tidy = TRUE) %>%
  na.omit(.) %>%
  mutate(AvgExpr = log10(baseMean)) %>% 
  rename(EnsemblID = row,
             logFC = log2FoldChange,
           p.value = pvalue,
           q.value = padj) %>%
  arrange(p.value) %>%
  select(EnsemblID, AvgExpr, logFC, p.value, q.value)


