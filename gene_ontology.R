library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)
library(dplyr)

if(exists("snakemake")) {
  file <- snakemake@input[["de_file"]]
  output_dir <- dirname(snakemake@output[["upregulated_GO"]])
} else {
  file <- "differential_expression_results.csv"
  output_dir <- "results/gene_ontology"
}


perform_GO_analysis <- function(gene_list, direction) {
  if (any(grepl("_", gene_list))) {
    gene_symbols <- sapply(strsplit(gene_list, "_"), function(x) {
      if (length(x) >= 2) x[2] else x[1]
    })
  } else {
    gene_symbols <- gene_list
  }
  
  head(gene_symbols)
  
  entrez_ids <- mapIds(
    org.Dm.eg.db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "ENTREZID"
  )
  
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) == 0) {
    message(paste("No ENTREZ IDs found for", direction))
    return(NULL)
  }
  
  go_enrich <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Dm.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  return(go_enrich)
}


de_data <- read.csv(file, row.names = 1)
head(de_data)
gene_names <- rownames(de_data)
head(gene_names)

# If log2FoldChange column exists adjust to logFC for edgeR compatibility
if("log2FoldChange" %in% colnames(de_data)) {
  colnames(de_data)[colnames(de_data) == "log2FoldChange"] <- "logFC"
}

up_condition1 <- gene_names[de_data$logFC > 0]
up_condition2 <- gene_names[de_data$logFC < 0]

analyses <- list(
  "upregulated" = up_condition1,
  "downregulated" = up_condition2
)

for (analysis_name in names(analyses)) {
  genes <- analyses[[analysis_name]]
  
  if (length(genes) == 0) {
    message(paste("No genes found for", analysis_name))
    next
  }
  
  go_result <- perform_GO_analysis(genes, analysis_name)
  
  if (!is.null(go_result)) {
    result_file <- file.path(output_dir, paste0(analysis_name, "_GO.csv"))
    write.csv(go_result@result, result_file, row.names = FALSE)
    
    tryCatch({
      if (nrow(go_result) > 0) {
        dot_plot <- dotplot(go_result, showCategory=15) + 
          ggtitle(paste("GO Analysis:", analysis_name))
        ggsave(file.path(output_dir, paste0(analysis_name, "_dotplot.png")), 
                dot_plot, width=10, height=8)
      }
    }, error = function(e) {
      message(paste("Error generating plot for", analysis_name, e$message))
    })
  }
}
