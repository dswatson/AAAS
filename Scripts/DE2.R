# Load libraries
library(data.table)
library(tximport)
library(DESeq2)
library(tidyverse)

# Import data
setwd('./Documents/QMUL/AAAS')
anno <- fread('./Data/Hs.anno.csv')
clin <- fread('./Data/clinical.csv', stringsAsFactors = TRUE) %>%
  mutate(Status = relevel(Status, ref = 'Healthy')) %>%
  filter(Mutation != 'Heterozygous_ex2')
t2g <- readRDS('./Data/Hs91.t2g.rds')
files <- file.path('./Data/Counts', clin$Idx, 'abundance.tsv')
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, importer = fread)

# Differential expression
dds <- DESeqDataSetFromTximport(txi, colData = clin, design = ~ Status)  
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, betaPrior = TRUE)
results(dds, coef = 2, tidy = TRUE) %>%
  na.omit(.) %>%
  mutate(AvgExpr = log10(baseMean)) %>% 
  rename(EnsemblID = row,
         logFC = log2FoldChange,
         p.value = pvalue,
         q.value = padj) %>%
  arrange(p.value) %>%
  inner_join(anno, by = 'EnsemblID') %>%
  select(EnsemblID, GeneSymbol, Description, 
         AvgExpr, logFC, p.value, q.value) %>%
  fwrite('./Results/Disease_vs_WT.csv')


