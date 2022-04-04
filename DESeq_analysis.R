
#### intalling and loading libraries
library(BiocParallel)
register(MulticoreParam(4))

library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

library(Biobase)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(DESeq2)
require(DOSE)
library(pathview)

# https://www.youtube.com/channel/UCnL-Lx5gGlW01OkskZL7JEQ
# https://www.youtube.com/watch?v=5tGCBW3_0IA
# https://www.youtube.com/watch?v=tlf6wYJrwKY
# https://youtu.be/UFB993xufUU
# https://pubmed.ncbi.nlm.nih.gov/31589133/

#### reading count matrix
counts = fread("data/GSE124326_count_matrix.txt")
colnames(counts) <- sapply(colnames(counts), function(x){sub(".counts","",x)})
counts[,1:3] %>% head


#### reading phenotypic data
pheno = fread("data/GSE124326_pheno.txt")
pheno %>% head

levels(factor(pheno$lithium)) # 0: non-lithium user, 1: lithium user
levels(factor(pheno$diagnosis)) # BP1: patients treated with lithum, BP2: patients not treated with lithium 


#### Checking and mapping counts and phenotypic data

phenoNames = pheno$sample
pheno = pheno[,-c("sample")]

rownames(pheno) <- phenoNames

countsNames <- counts$gene
counts = counts[,-c("gene")]
rownames(counts) <- countsNames

#Remove the samples in counts that are missing in pheno 
missingNames <- setdiff(colnames(counts), rownames(pheno)) # elements of x not in y
missingNames

counts <- counts[, !c(..missingNames)]
rownames(counts) <- countsNames

#check samples of counts match with those phenotypic information
all(colnames(counts) == rownames(pheno))

# combining lithium and non-lithium treated BP patients
pheno$diagnosis[pheno$diagnosis=="BP1"]<-"BP"
pheno$diagnosis[pheno$diagnosis=="BP2"]<-"BP"

# Converting to factors
pheno$diagnosis <- factor(pheno$diagnosis)
pheno$lithium <- factor(pheno$lithium)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno,
                              design = ~ diagnosis)

mcols(dds) <- DataFrame(mcols(dds), data.frame(gene=rownames(counts)))
mcols(dds)
dds <- dds[rowSums(counts(dds)) >= 20,]

rnames <- rowData(dds)$gene
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05)
res$gene <- rnames
res <- res[order(res$padj),]
res <- na.omit(res)

summary(res)

ggplot(as(res, "data.frame"), aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "Royalblue", boundary = 0)

plotMA(dds, ylim = c( -2, 2))

plotDispEsts(dds)

write.csv(res, "data/deseq_out.csv")

#res[(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),]