
#### intalling and loading libraries
library(BiocParallel)
register(MulticoreParam(4))

library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

library(Biobase)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75)
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(DESeq2)
require(DOSE)
library(pathview)


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

############ DEG analysis

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

# PlotMA log2 fold changes attributable to a given variable 
# over the mean of normalized counts for all the samples 
plotMA(res, ylim = c( -2, 2))

plotDispEsts(dds)

write.csv(res, "data/deseq_out.csv")

############ significant DEG 

res_filtered <- res[(abs(res$log2FoldChange) > 1 & res$padj < 0.05),]

ggplot(as(res_filtered, "data.frame"), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()

top25_genes <- res_filtered[order(abs(res_filtered$log2FoldChange), decreasing = T),][1:25,]$gene
top25_genes

write.csv(res_filtered, "data/deseq_out_filtered.csv")

############ GO enrichment analysis

res_filtered$gene <- gsub("\\..*", "", res_filtered$gene)
res_filtered$SYMBOL <- mapIds(EnsDb.Hsapiens.v75,
                              keys = res_filtered$gene,
                              column = "SYMBOL",
                              keytype = "GENEID",
                              multiVals = "first")

res_filtered$PATH <- mapIds(org.Hs.eg.db,
                              keys = res_filtered$SYMBOL,
                              column = "PATH",
                              keytype = "SYMBOL",
                              multiVals = "first")

res_filtered$UNIPROT <- mapIds(org.Hs.eg.db,
                            keys = res_filtered$SYMBOL,
                            column = "UNIPROT",
                            keytype = "SYMBOL",
                            multiVals = "first")

humanGeneUniverse <- as.character(unique(select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), column="SYMBOL")$SYMBOL))
geneList <- factor(as.integer(humanGeneUniverse %in% res_filtered$SYMBOL))
names(geneList) <- res_filtered$SYMBOL

GOdata <- new("topGOdata", 
              ontology="BP", 
              allGenes=geneList, 
              nodeSize=5,
              annotationFun=annFUN.org, 
              mapping="org.Hs.eg.db", 
              ID="symbol")

GOresult <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultTable <- GenTable(GOdata, GOresult, topNodes = 50, numChar=200)
resultTable

pValue.classic <- score(GOresult)
head(pValue.classic)

showSigOfNodes(GOdata, score(GOresult), firstSigNodes = 5, useInfo = 'all')
