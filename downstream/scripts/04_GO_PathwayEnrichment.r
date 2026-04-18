library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(enrichplot)



# ------------------------------------------------
#            Gene Onthology (GO) on common Genes
# ------------------------------------------------ 

## FUNCTIONS

# 01-Load data

clean_gene_list <- function(x) {
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}
# 03-mapIds()

mapping <- function(genes){
  
  if (is.null(genes) || length(genes) == 0){
    cat("There is a problem with genes list.\n")
    return(NULL)
  }
  
  id <- mapIds(
    org.Hs.eg.db,
    keys = genes,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  id <- id[!is.na(id)]
  id <- unique(id)
  
  return(id)
}

# 02-enrichGO()

go_terms <- function(genes, onto, background){
  
  if (is.null(genes) || is.null(background)) {
    cat("There is a problem with genes or background.\n")
    return(NULL)
  }
  
  terms <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = onto,
    universe = background,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  return(terms)
}

# 03-enrichKEGG()

kegg_func <- function(genes, background){
  
  if (is.null(genes) || is.null(background)) {
    cat("There is a problem with genes or background.\n")
    return(NULL)
  }
  
  pathways <- enrichKEGG(
    gene = genes,
    organism = "hsa",
    universe = background, 
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
  )
  
  return(pathways)
}

# 04-enrichPathway()

reactome_func <- function(genes, background){
  
  if (is.null(genes) || is.null(background)) {
    cat("There is a problem with genes or background.\n")
    return(NULL)
  }
  pathways <- enrichPathway(
    gene = genes,
    organism = "human",
    universe = background,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
    )
  
  return(pathways)
  
}

# 05-barplot function

plotting1 <- function(df, title){
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  Top <- df |>
    dplyr::group_by(Ontology) |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = 20) |>
    dplyr::ungroup() |> 
    dplyr::arrange(Ontology, desc(Count)) |>
    dplyr::ungroup()
  
  Top$Description <- factor(
    Top$Description,
    levels = rev(unique(Top$Description))
  )
  
  ggplot(Top, aes(
    x = Count,
    y = Description
  )) +
    geom_col(aes(fill = -log10(p.adjust)), width = 0.7) +
    facet_wrap(~Ontology, scales = "free_y", ncol = 1) +
    scale_fill_gradient(
      low = "steelblue",
      high = "firebrick",
      name = "-log10(FDR)"
    ) +
    theme_classic() +
    labs(
      title = title,
      x = "Gene count",
      y = NULL,
      fill = "-log10(FDR)"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
}

plotting2 <- function(df, title){
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  Top <- df |>
    dplyr::group_by(Database) |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = 20) |>
    dplyr::arrange(Database, desc(Count)) |>
    dplyr::ungroup()
  
  Top$Description <- factor(
    Top$Description,
    levels = rev(unique(Top$Description))
  )
  
  ggplot(Top, aes(
    x = Count,
    y = Description
  )) +
    geom_col(aes(fill = -log10(p.adjust)), width = 0.7) +
    facet_wrap(~Database, scales = "free_y", ncol = 1) +
    scale_fill_gradient(
      low = "steelblue",
      high = "firebrick",
      name = "-log10(FDR)"
    ) +
    theme_classic() +
    labs(
      title = title,
      x = "Gene count",
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
}

# ------------------------------------------------
#            Start Analysis
# ------------------------------------------------ 

# 01-Load list of common DEGs 

Gene_ID_bg <- clean_gene_list(readLines("background_genes.txt"))
DEG_common_up <- clean_gene_list(readLines("DEG_common_up.txt"))
DEG_common_down <- clean_gene_list(readLines("DEG_common_down.txt"))


cat("Background genes:", length(Gene_ID_bg), "\n")
cat("DEG up:", length(DEG_common_up), "\n")
cat("DEG down:", length(DEG_common_down), "\n")
cat("Tot DEGs:", length(DEG_common_up) + length(DEG_common_down))

# 02-Take EntrezID

background <- mapping(Gene_ID_bg)

genes_upregulated <- mapping(DEG_common_up)

genes_downregulated <- mapping(DEG_common_down)

# Check
if(all(c(
  !is.null(background),
  !is.null(genes_upregulated),
  !is.null(genes_downregulated)
))) message("All gene lists correctly mapped to EntrezID")

# 03-Search GO terms
# Up-regulated genes 

GO_up_BP <- go_terms(genes_upregulated, "BP", background)
GO_up_MF <- go_terms(genes_upregulated, "MF", background)
GO_up_CC <- go_terms(genes_upregulated, "CC", background)

cat("BP terms up-regulated genes:", nrow(as.data.frame(GO_up_BP)), "\n")
cat("MF terms up-regulated genes:", nrow(as.data.frame(GO_up_MF)), "\n")
cat("CC terms up-regulated genes:", nrow(as.data.frame(GO_up_CC)), "\n")

# Down-regulated genes

GO_down_BP <- go_terms(genes_downregulated, "BP", background)
GO_down_MF <- go_terms(genes_downregulated, "MF", background)
GO_down_CC <- go_terms(genes_downregulated, "CC", background)

cat("BP terms for down-regulated genes:", nrow(as.data.frame(GO_down_BP)), "\n")
cat("MF terms for down-regulated genes:", nrow(as.data.frame(GO_down_MF)), "\n")
cat("CC terms for down-regulated genes:", nrow(as.data.frame(GO_down_CC)), "\n")

## Plotting
# Biological Process

GO_down_BP <- setReadable(
  GO_down_BP,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

GO_down_BP <- pairwise_termsim(GO_down_BP)

pdf("../plots/ORA/01-GO_analysis_BiologicalProcess.pdf", width = 12, height = 9)
dotplot(GO_down_BP, showCategory = 20, title = "GO Biological Process enrichment analysis (common downregulated genes)")
emapplot(GO_down_BP, layout = "fr") + ggtitle("GO Biological Process similarity network (common downregulated genes)")
cnetplot(GO_down_BP, showCategory = 10) + ggtitle("Gene–GO term interaction network (common downregulated genes)")
dev.off()

# Molecular Functions

GO_down_MF <- setReadable(
  GO_down_MF,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

GO_down_MF <- pairwise_termsim(GO_down_MF)

pdf("../plots/ORA/02-GO_analysis_MolecularFunctions.pdf", width = 12, height = 9)
dotplot(GO_down_MF, showCategory = 20, title = "GO Molecular Functions enrichment analysis (common downregulated genes)")
emapplot(GO_down_MF, layout = "fr") + ggtitle("GO Molecular Functions similarity network (common downregulated genes)")
cnetplot(GO_down_MF, showCategory = 10) + ggtitle("Gene–GO term interaction network (common downregulated genes)")
dev.off()

# Cellular Components 

GO_down_CC <- setReadable(
  GO_down_CC,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

GO_down_CC <- pairwise_termsim(GO_down_CC)

pdf("../plots/ORA/03-GO_analysis_CellularComponents.pdf", width = 12, height = 9)
dotplot(GO_down_CC, showCategory = 20, title = "GO Cellular Components enrichment analysis (common downregulated genes)")
emapplot(GO_down_CC, layout = "fr") + ggtitle("GO Cellular Components similarity network (common downregulated genes)")
cnetplot(GO_down_CC, showCategory = 10) + ggtitle("Gene–GO term interaction network (common downregulated genes)")
dev.off()

# barplot of all ontology

GO_down_BP@result$Ontology <- "BP"
GO_down_MF@result$Ontology <- "MF"
GO_down_CC@result$Ontology <- "CC"

GO_all <- bind_rows(
  as.data.frame(GO_down_BP),
  as.data.frame(GO_down_CC),
  as.data.frame(GO_down_MF)
)

plot_all_GO <- plotting1(GO_all, "GO enrichment analysis of differentially expressed genes (common downregulated genes)")
ggsave("../plots/ORA/04_GO_analyisis_complete.pdf", plot = plot_all_GO, width = 14, height = 14)



# ------------------------------------------------
#            Pathway Enrichment
# ------------------------------------------------ 

# 03-Pathway Enrichment using KEGG and Reactome databases
# KEGG database

kegg_up <- kegg_func(genes_upregulated, background)
kegg_down <- kegg_func(genes_downregulated, background)

cat("KEGG pathways for up-regulated genes:", nrow(as.data.frame(kegg_up)), "\n")
cat("KEGG pathways for down-regulated genes:", nrow(as.data.frame(kegg_down)), "\n")

# Reactome database

reactome_up <- reactome_func(genes_upregulated, background)
reactome_down <- reactome_func(genes_downregulated, background)

cat("Reactome pathways for up-regulated genes:", nrow(as.data.frame(reactome_up)), "\n")
cat("Reactome pathways for down-regulated genes:", nrow(as.data.frame(reactome_down)), "\n")

# 04-Visualization of pathway enrichment
# KEGG

kegg_down <- setReadable(
  kegg_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_down <- pairwise_termsim(kegg_down)

pdf("../plots/ORA/05-KEGG_pathways.pdf", width = 10, height = 8)
dotplot(kegg_down, showCategory = 20, title = "KEGG pathway enrichment analysis (common downregulated genes)")
emapplot(kegg_down, layout = "fr") + ggtitle("KEGG pathway similarity network (common downregulated genes)")
cnetplot(kegg_down, showCategory = 10) + ggtitle("Gene–KEGG pathway interaction network (common downregulated genes)")
dev.off()

# Reactome

reactome_down <- setReadable(
  reactome_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

reactome_down <- pairwise_termsim(reactome_down)

pdf("../plots/ORA/06-Reactome_pathways.pdf", width = 14, height = 9)
dotplot(reactome_down, showCategory = 20, title = "Reactome pathway enrichment analysis (common downregulated genes)")
emapplot(reactome_down) + ggtitle("Reactome pathway similarity network (common downregulated genes)")
cnetplot(reactome_down, showCategory = 10) + ggtitle("Gene–Reactome pathway interaction network (common downregulated genes)")
dev.off()

# barplot of both databases

kegg_down@result$Database <- "KEGG"
reactome_down@result$Database <- "Reactome" 

db_all <- bind_rows(
  as.data.frame(kegg_down),
  as.data.frame(reactome_down)
)

plot_all_db <- plotting2(db_all, "Integrated pathway enrichment analysis of common downregulated genes (KEGG and Reactome)")
ggsave("../plots/ORA/07_PathwayEnrichment_complete.pdf", plot = plot_all_db, width = 14, height = 14)

message("===================================")
cat("Over representation analysis completed successfully!\n")


