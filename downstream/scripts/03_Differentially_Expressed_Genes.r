library(DESeq2)
library(AnnotationDbi)
library(EnhancedVolcano)
library(pheatmap)
library(org.Hs.eg.db)
library(tidyverse)



# ------------------------------------------------
#            Differentially Expressed Genes (DEGs)
# ------------------------------------------------

dds <- readRDS("dds.rds")

if (!exists("dds")) {
  stop("ERROR: dds object not loaded. Check readRDS() step.")
}

if (is.null(dds)) {
  stop("ERROR: dds is NULL")
}

message("Loaded DESeq2 object with ",
        nrow(dds), " genes and ",
        ncol(dds), " samples")

vsdata <- vst(dds, blind = FALSE)
colData <- as.data.frame(colData(dds))["condition"]

# 01-Get the normalized results form DESeq object and get Symbol Gene

siOGT1 <- results(dds, contrast = c("condition", "siOGT1", "siNC"), alpha = 0.05)
siOGT2 <- results(dds, contrast = c("condition", "siOGT2", "siNC"), alpha = 0.05)

summary(siOGT1)
summary(siOGT2)

siOGT1_res <- as.data.frame(siOGT1)
siOGT1_res <- siOGT1_res[!is.na(siOGT1_res$padj),]
siOGT1_res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(siOGT1_res),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

siOGT2_res <- as.data.frame(siOGT2)
siOGT2_res <- siOGT2_res[!is.na(siOGT2_res$padj),]
siOGT2_res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(siOGT2_res),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

stopifnot(all(rownames(vsdata) %in% rownames(dds)))

# 02-MA Plot 

pdf("../plots/DE/01_MA_siOGT.pdf")

plotMA(
  siOGT1,
  ylim = c(-7, 7),
  alpha = 0.05,
  main = "MA plot – siOGT1 vs Control"
)
abline(h = c(-1, 1), col = "grey50", lty = 2)
plotMA(siOGT2, 
       ylim = c(-7,7),
       alpha = 0.05,
       main = "MA plot – siOGT2 vs Control" 
)
abline(h = c(-1,1), col="grey50", lty = 2)

dev.off()

# 03-VolcanoPlot_vs01

pdf("../plots/DE/02_EnhancedVolcano_siOGT.pdf")

EnhancedVolcano(
  siOGT1_res, 
  lab = siOGT1_res$symbol, 
  x = "log2FoldChange", 
  y = "padj",
  title = "Volcano plot – siOGT1 vs siNC",
  subtitle = "All genes tested (DESeq2 results)"
)

EnhancedVolcano(
  siOGT2_res, 
  lab = siOGT2_res$symbol,
  x = "log2FoldChange", 
  y = "padj",
  title = "Volcano plot – siOGT2 vs siNC",
  subtitle = "All genes tested (DESeq2 results)"
)

dev.off()

# 03-Significant genes

siOGT1_res$threshold <- abs(siOGT1_res$log2FoldChange) > 0.5 & siOGT1_res$padj < 0.05

siOGT2_res$threshold <- abs(siOGT2_res$log2FoldChange) > 0.5 & siOGT2_res$padj < 0.05 

# 04-VolcanoPlot_vs02 for significant genes

v1 <- ggplot(siOGT1_res, aes(x = log2FoldChange,
                     y = -log10(padj),
                     colour = threshold)) + 
  geom_point(alpha = 0.6, size = 1.6, na.rm = TRUE) + 
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  labs(
    title = "Volcano plot – siOGT1 vs siNC",
    subtitle = "All genes tested; DEGs highlighted (adjusted p-value < 0.05, |log2FC| > 0.5)",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    colour = "Significant"
  ) + 
  coord_cartesian() + 
  theme_classic() + 
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

v2 <- ggplot(siOGT2_res, aes(x = log2FoldChange,
                     y = -log10(padj),
                     colour = threshold)) + 
  geom_point(alpha = 0.6, size = 1.6, na.rm = TRUE) + 
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  labs(
    title = "Volcano plot – siOGT2 vs siNC",
    subtitle = "All genes tested; DEGs highlighted (adjusted p-value < 0.05, |log2FC| > 0.5)",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    colour = "Significant"
  ) + 
  coord_cartesian() + 
  theme_classic() + 
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

pdf("../plots/DE/03_VolcanoPlot_Highlighted_DEG.pdf")
print(v1)
print(v2)
dev.off()

# 05-Find DEGs

# siOGT1
DEG_siOGT1 <- siOGT1_res[siOGT1_res$threshold == TRUE,]

DEG_siOGT1_up <- rownames(DEG_siOGT1[DEG_siOGT1$log2FoldChange > 0, ])

DEG_siOGT1_down <- rownames(DEG_siOGT1[DEG_siOGT1$log2FoldChange < 0, ])

# siOGT2
DEG_siOGT2 <- siOGT2_res[siOGT2_res$threshold == TRUE,]

DEG_siOGT2_up <- rownames(DEG_siOGT2[DEG_siOGT2$log2FoldChange > 0, ])

DEG_siOGT2_down <- rownames(DEG_siOGT2[DEG_siOGT2$log2FoldChange < 0, ])


# Common DEGs

DEG_common <- intersect(rownames(DEG_siOGT1), rownames(DEG_siOGT2))
DEG_common_up <- intersect(DEG_siOGT1_up, DEG_siOGT2_up)
DEG_common_down <- intersect(DEG_siOGT1_down, DEG_siOGT2_down)

# 06-Heatmap of DEGs 

DEG_matrix1 <- assay(vsdata)[rownames(DEG_siOGT1),]
DEG_matrix2 <- assay(vsdata)[rownames(DEG_siOGT2), ]
DEG_common_matrix <- assay(vsdata)[DEG_common,]

stopifnot(rownames(colData) == colnames(vsdata))

# siOGT1

pdf("../plots/DE/04_Heatmaps_DEGs.pdf", width = 10, height = 10)

pheatmap(
  DEG_matrix1, 
  scale="row", 
  color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
  annotation_col = colData,
  show_rownames = FALSE,
  main = "siOGT1: DEGs padj < 0.05 & |log2FC| > 0.5"
  )

# siOGT2

pheatmap(
  DEG_matrix2,
  scale="row", 
  color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
  annotation_col = colData,
  show_rownames = FALSE,
  main = "siOGT2: DEGs padj < 0.05 & |log2FC| > 0.5"
  )

# Common DEGs

pheatmap(
  DEG_common_matrix, 
  scale="row", 
  color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
  annotation_col = colData,
  show_rownames = FALSE,
  main = "Common Genes for siOGT1 and siOGT2: DEGs padj < 0.05 & |log2FC| > 0.5"
  )

dev.off()

# 07-Save DEG common 

writeLines(DEG_common_up, "../objects/DEG_common_up.txt") #for ORA
writeLines(DEG_common_down, "../objects/DEG_common_down.txt") #for ORA

message("===================================")
cat("Differential expression analysis completed successfully!\n")
cat(
  "siOGT1-upregulated: ", length(DEG_siOGT1_up), "\n",
  "siOGT1-downregulated: ", length(DEG_siOGT1_down), "\n",
  "Total DEG for siOGT1: ", length(DEG_siOGT1_up) + length(DEG_siOGT1_down), "\n",
  "\n",
  "siOGT2-upregulated: ", length(DEG_siOGT2_up), "\n",
  "siOGT2-downregulated: ", length(DEG_siOGT2_down), "\n",
  "Total DEG for siOGT2: ", length(DEG_siOGT2_up) + length(DEG_siOGT2_down), "\n",
  "\n",
  "DEG common upregulated: ", length(DEG_common_up), "\n",
  "DEG common downregulated: ", length(DEG_common_down), "\n",
  "Total common DEG: ", length(DEG_common_up) + length(DEG_common_down)
)
