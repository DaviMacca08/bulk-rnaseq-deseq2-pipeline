library(DESeq2)
library(ggplot2)
library(GEOquery)



# ------------------------------------------
#            Building DESeq Object
# ------------------------------------------

# 01-Read counts.txt from featureCounts

counts <- read.table("counts.txt", sep = "\t", header = T, skip = 1)
counts <- counts[, -c(2:6)]
rownames(counts) <- counts$Geneid 
counts <- counts[, -1]  
colnames(counts) <- gsub(".*SRR([0-9]+)_.*", "SRR\\1", colnames(counts))

gse <- getGEO("GSE316021", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

#check
message("Counts matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

stopifnot(ncol(counts) > 1)
stopifnot(nrow(counts) > 1000)

stopifnot(!any(duplicated(rownames(counts))))
stopifnot(!any(duplicated(colnames(counts))))

stopifnot(!any(is.na(counts)))
stopifnot(all(is.finite(as.matrix(counts))))

# 02-Set the correct samples names (siOGT --> siRNA, siNC --> Negative CTRL) 

#extract sample names
title <- pheno$title

sample_info <- gsub("KGN cells, ", "", title)
sample_info <- gsub(", 48h", "", sample_info)
sample_info <- gsub("-", "_", sample_info)
sample_info <- rev(sample_info)

colnames(counts) <- sample_info

colData <- data.frame(
  condition = factor(gsub("_.*", "", sample_info)),
  replicate = gsub(".*_", "", sample_info)
)

rownames(colData) <- sample_info

#check
message("Conditions:")
print(table(colData$condition))

stopifnot(all(colnames(counts) == rownames(colData)))
stopifnot(length(unique(colData$condition)) >= 2)

# 03-Build a colData with samples info and check the names and orders

stopifnot(all(colnames(counts) %in% rownames(colData)))
stopifnot(all(colnames(counts) == rownames(colData)))

# 04-Build a DESeq Dataset object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

#low expression filtering
dds <- dds[rowSums(counts(dds)) > 10, ]

# 05-Set the factor level and run DESeq function

dds$condition <- relevel(dds$condition, ref = "siNC")
dds <- DESeq(dds)
res <- results(dds)

message("DE genes (padj < 0.05): ", sum(res$padj < 0.05, na.rm = TRUE))

stopifnot(!is.null(dispersions(dds)))

# 06-Variance stabilization using vst()

vsdata <- vst(dds, blind = FALSE)

pdf("../../../downstream/plots/QC/01_PCA_VstData.pdf", width = 8, height = 6)
pca <- plotPCA(vsdata, intgroup = "condition")
print(
  pca + ggtitle("PCA - KGN cells, siOGT vs Control")
)
dev.off()

pdf("../../../downstream/plots/QC/02_DispEsts_dds.pdf", width = 8, height = 6)
plotDispEsts(dds,  main = "Dispersion Estimates - KGN cells, siOGT vs Control") 
dev.off()

# 07-Save output and final messages 

saveRDS(dds, "../../../downstream/objects/dds.rds")
writeLines(rownames(vsdata), "../../../downstream/objects/background_genes.txt") #for ORA
write.csv(as.data.frame(res), "../../../downstream/objects/DESeq2_results.csv")

message("===================================")
cat("DESeq2 analysis completed successfully!\n")
cat("DESeq2 object built with ", nrow(dds), " genes and ", ncol(dds), " samples")
cat("Condition: ", paste(levels(dds$condition), collapse = " vs "))

