# intalling and loading libraries
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

# reading count matrix
counts = fread("data/GSE124326_count_matrix.txt")
colnames(counts) <- sapply(colnames(counts), function(x){sub(".counts","",x)})
counts[,1:3] %>% head


# reading phenotypic data
pheno = fread("data/GSE124326_pheno.txt")
pheno %>% head

levels(factor(pheno$lithium)) # 0: non-lithium user, 1: lithium user
levels(factor(pheno$diagnosis)) # BP1: patients treated with lithum, BP2: patients not treated with lithium 


#  Checking and mapping counts and phenotypic data

phenoNames = pheno$sample
pheno = pheno[,-c("sample")]

rownames(pheno) <- phenoNames

countsNames <- counts$gene
counts = counts[,-c("gene")]
rownames(counts) <- countsNames

# Remove the samples in counts that are missing in pheno 
missingNames <- setdiff(colnames(counts), rownames(pheno)) # elements of x not in y
missingNames

counts <- counts[, !c(..missingNames)]
rownames(counts) <- countsNames

# check samples of counts match with those phenotypic information
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
dds <- DESeq(dds, parallel = T)
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

res_filtered <- res[(abs(res$log2FoldChange) >= 1 & res$padj <= 0.05),]

ggplot(as(res_filtered, "data.frame"), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point()

top25_genes <- res_filtered[order(abs(res_filtered$log2FoldChange), decreasing = T),][1:25,]$gene
top25_genes

write.csv(res_filtered, "data/deseq_out_filtered.csv")

############ GO enrichment analysis

res$gene <- gsub("\\..*", "", res$gene)
res$SYMBOL <- mapIds(EnsDb.Hsapiens.v75,
                     keys = res$gene,
                     column = "SYMBOL",
                     keytype = "GENEID",
                     multiVals = "first")

res$ENTREZID <- mapIds(org.Hs.eg.db,
                     keys = res$SYMBOL,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

res$PATH <- mapIds(org.Hs.eg.db,
                   keys = res$SYMBOL,
                   column = "PATH",
                   keytype = "SYMBOL",
                   multiVals = "first")

res$UNIPROT <- mapIds(org.Hs.eg.db,
                      keys = res$SYMBOL,
                      column = "UNIPROT",
                      keytype = "SYMBOL",
                      multiVals = "first")

humanGeneUniverse <- as.character(unique(select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), column="SYMBOL")$SYMBOL))

geneList <- factor(as.integer(humanGeneUniverse %in% res$SYMBOL))
names(geneList) <- res$SYMBOL

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

top_50 <- resultTable$GO.ID
showSigOfNodes(GOdata, score(GOresult), firstSigNodes = 5, useInfo = 'all')

############ GS enrichment analysis

df <- res

original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$SYMBOL

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe having only the genes that were successfully mapped
df2 = df[which(dedup_ids$SYMBOL %in% df$SYMBOL),]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


enrichplot::dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# pathways ridge plot (Helpful to interpret up/down-regulated pathways.)
enrichplot::ridgeplot(kk2) + labs(x = "enrichment distribution")

# top 2 enriched KEGG pathways diagram
top2_id <- kk2[order(kk2$enrichmentScore, decreasing = T),][1:2,]$ID

pathview(gene.data=kegg_gene_list, pathway.id=top2_id[1], species = kegg_organism)
pathview(gene.data=kegg_gene_list, pathway.id=top2_id[2], species = kegg_organism)

### Prepare dataset
res_filtered <- read.csv("data/deseq_out_filtered.csv")
res_filtered <- subset(res_filtered, select = -c(X))

rownames(counts) <- gsub("\\..*", "", rownames(counts))

counts$gene <- rownames(counts)
pheno$sample <- rownames(pheno)

top_100_genes <- res_filtered[1:100,]$gene
top_100 <- res_filtered[1:100,]

top_100_counts <- counts[which(top_100_genes %in% counts$gene),]
rownames(top_100_counts) <- top_100_counts$gene

dds100 = dds[rownames(top100,)]

############ Classification
library(MLSeq)


# Train Test Split
n <- ncol(dds100) # number of samples
p <- nrow(dds100) # number of features=
class = data.frame(Bipolar=as.numeric(dds100$Bipolar), row.names = colnames(dds100))

# number of samples for test set (30% test, 70% train).
nTest <- ceiling(n*0.3)
set.seed(2)
ind <- sample(n, nTest, FALSE)

# train and test split
data.train <- as.matrix(counts_100[,-ind]+1)
data.test <- as.matrix(as.matrix(counts_100[,ind]+1))
classtr <- factor(class[-ind,])
classtr <- data.frame(condition=classtr, row.names=colnames(data.train))
classts <- factor(class[ind,])
classts <- DataFrame(condition=classts, row.names=colnames(data.test))

# Design = ~ 1 indicates that there is no need to create a differential testing design
# like before
data.trainS4 <- DESeqDataSetFromMatrix(countData=data.train, colData=classtr,
                                       design= ~ 1)
data.testS4 <- DESeqDataSetFromMatrix(countData=data.test, colData=classts,
                                      design= ~ 1)
######## svm radial
# The control variable sets the method for cross validation, number of partitions and how
# many times to repeat the process. classProbs gives us the probabilities that an observation
# belongs to different classes
ctrl <- trainControl(method = "repeatedcv",
                     number=2,
                     repeats=2,
                     classProbs = FALSE)
fit.svm = classify(data=data.trainS4, method='svmRadial', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.svm)
confusionMat(fit.svm)

availableMethods()

######### lda
fit.lda = classify(data=data.trainS4, method='lda', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.lda)
confusionMat(fit.lda)

######### glm
fit.glm = classify(data=data.trainS4, method='glm', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.glm)
confusionMat(fit.glm)

######### neural networks
fit.nnet = classify(data=data.trainS4, method='nnet', preProcessing = 'deseq-rlog',
                    control = ctrl)
show(fit.nnet)
confusionMat(fit.nnet)

######### Decision trees
fit.rpart = classify(data=data.trainS4, method='rpart', preProcessing = 'deseq-rlog',
                     control = ctrl)
show(fit.rpart)
confusionMat(fit.rpart)

######### Gradient Boost Machine
fit.gbm = classify(data=data.trainS4, method='gbm', preProcessing = 'deseq-rlog',
                   control = ctrl)
show(fit.gbm)
confusionMat(fit.gbm)


################# Parameter Tuning
set.seed(2)
# Here tuneLength sets the number of levels for parameter tuning
# ref="1" means that our reference class in the target variable (Bipolar) is 1 which
# equivalent for being a patient of bipolar disorder
######### Support vector machines with radial basis function kernel
fit.svm <- classify(data = data.trainS4, method = "svmRadial",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.svm)
######### LDA
fit.lda <- classify(data = data.trainS4, method = "lda",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.lda)
######### Generalized Linear Model
fit.glm <- classify(data = data.trainS4, method = "glm",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.glm)
######### Neural Network
fit.nnet <- classify(data = data.trainS4, method = "nnet",
                     preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                     control = ctrl)
show(fit.nnet)
######### Decision Tree with caret
fit.rpart <- classify(data = data.trainS4, method = "rpart",
                      preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                      control = ctrl)
show(fit.rpart)
######### Gradient Boost Machines
fit.gbm <- classify(data = data.trainS4, method = "gbm",
                    preProcessing = "deseq-vst", ref = "1", tuneLength = 10,
                    control = ctrl)
show(fit.gbm)

