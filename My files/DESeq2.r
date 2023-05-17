library("DESeq2")
library("apeglm")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")

cts <- as.matrix(read.csv("/Users/christopherburtner/Documents/Molecular Longevity Lab/RNAseq/my R analysis/genecounts.csv", row.names = "gene_id"))
coldata <- read.csv("/Users/christopherburtner/Documents/Molecular Longevity Lab/RNAseq/my R analysis/samples.csv", row.names = 1)
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = cts, colData= coldata, design = ~ condition)
dds
dds <- DESeq(dds)

#log2 fold change shrinkage eliminates genes with very few reads that are very unlikely to be statistically different. This reduces the number of tests in multiple comparison. It is advised you use lfcShrink for results.
res <- results(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")

get_significant_genes <- function(df, pvalue = 0.05, l2fc = 1) {
  sig <- df[df$padj <= pvalue, ]
  up <- sig[sig$log2FoldChange >= l2fc, ]
  down <- sig[sig$log2FoldChange <= -l2fc, ]
  return(list("up"=up, "down"=down))
}

write.csv(na.omit(res), file="shrunk_deseq_results_non_nan.csv", quote=FALSE)
sig_genes <- get_significant_genes(na.omit(res), 1)
write.csv(sig_genes$up, file="upregulated.csv", quote=FALSE)
write.csv(sig_genes$down, file="downregulated.csv", quote=FALSE)

# variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Create heat map from data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatMap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA graphing function
pcaData <- plotPCA(vsd, returnData=TRUE)

# get percent variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

# pca code
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition, name = name)) +
    stat_ellipse(aes(group = condition, linetype = condition), type = "t", level = 0.95, size = 1.25, show.legend = FALSE) +
    geom_point(size = 3, show.legend = TRUE) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + theme_bw() + ggtitle("PCA") +
    guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) 
