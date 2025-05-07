library("DESeq2")
library("RColorBrewer")
library("pheatmap")

# Create empty pdf 
pdf("deseq2.pdf")

# Path to directory of htseq outputs 
args <- commandArgs(trailingOnly = TRUE)
directory <- './7_deseq2_inputs/'

# Retrieve htseq count files
sampleFiles <- grep("treated",list.files(directory),value=TRUE)

# Retrieve condition (treated/nontreated) for each file
sampleCondition <- sub(".*_(treated|nontreated).*", "\\1", sampleFiles)

# Create a table of sample metadata 
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
# Level                          
sampleTable$condition <- factor(sampleTable$condition)

# Read htseq datasets 
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

# Specify the reference level:
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "nontreated")

# DE analysis
ddsHTSeq <- DESeq(ddsHTSeq, parallel=TRUE)
res <- results(ddsHTSeq)

# Shrinkage
resLFC <- lfcShrink(ddsHTSeq, coef="condition_treated_vs_nontreated", type="apeglm")

# Export isomiRs with 
## adjusted p-values were less than 0.05 and abs(logFC) > 1
DEisomiRs_idx <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
DEisomiRs <- row.names(res)[DEisomiRs_idx]
writeLines(c(DEisomiRs), "DE isomiRs.txt")


vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=FALSE)
# PCA
plotPCA(vsd, intgroup=c("condition"))

# Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Dispersion estimate
plotDispEsts(ddsHTSeq)

# MA plot
plotMA(resLFC, ylim=c(-2,2))

# P value
h1 <- hist(res$pvalue, breaks=0:50/50, plot=FALSE)
plot(h1, main="P-value Histogram", xlab="P-value", ylab="Frequency", col='powderblue')

dev.off()
