library(edgeR)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(org.Dm.eg.db)
library(AnnotationDbi)

set.seed(123)

if(exists("snakemake")) {
  counts_file <- snakemake@input[["counts"]]
  sample_info_file <- snakemake@input[["sample_info"]]
  output_dir <- dirname(snakemake@output[["normalized_counts"]])
} else {
  counts_file <- "counts.txt"
  sample_info_file <- "samples.csv"
  output_dir <- "results/differential_expression"
}

if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

count_data <- read.delim(counts_file, comment.char="#")
count_matrix <- count_data %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(6:ncol(.))

sample_info <- read.csv(sample_info_file)
sample_info$condition <- factor(sample_info$condition)
colnames(count_matrix) <- sample_info$sample

flybase_ids <- rownames(count_matrix)
id_map <- mapIds(org.Dm.eg.db,
                 keys = flybase_ids,
                 column = "SYMBOL",
                 keytype = "FLYBASE",
                 multiVals = "first")

symbol_count_matrix <- count_matrix
rownames(symbol_count_matrix) <- ifelse(is.na(id_map), flybase_ids, id_map)

symbol_count_matrix <- symbol_count_matrix %>%
  rownames_to_column("symbol") %>%
  group_by(symbol) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("symbol")

group <- sample_info$condition
dge <- DGEList(counts = symbol_count_matrix, group = group)
keep <- rowSums(cpm(dge) > 1) >= 3
table(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

conditions <- levels(group)

# Normalized counts (CPM)
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))



# PCA plot
logCPM <- cpm(dge, log=TRUE, prior.count=2)
pca <- prcomp(t(logCPM))
pca_data <- as.data.frame(pca$x)
pca_data$condition <- group
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))[1:2]

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pca_data)), show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA Plot")

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# Perform all pairwise comparisons
all_comparisons <- combn(conditions, 2, simplify = FALSE)

# Create a list to store all results
all_results <- list()

for (i in seq_along(all_comparisons)) {
  comp <- all_comparisons[[i]]
  comp_name <- paste(comp[1], "vs", comp[2], sep = "_")
  
  # Create contrast
  contrast_vec <- rep(0, ncol(design))
  contrast_vec[which(colnames(design) == paste0("group", comp[1]))] <- 1
  contrast_vec[which(colnames(design) == paste0("group", comp[2]))] <- -1
  
  # Test for differential expression
  qlf <- glmQLFTest(fit, contrast = contrast_vec)
  topTagsRes <- topTags(qlf, n = nrow(dge))
  res_df <- as.data.frame(topTagsRes) %>% rownames_to_column("gene")
  
  # Add comparison info
  res_df$comparison <- comp_name
  res_df$condition1 <- comp[1]
  res_df$condition2 <- comp[2]
  
  # Store results
  all_results[[comp_name]] <- res_df
  
  # Create comparison-specific directory
  comp_dir <- file.path(output_dir, comp_name)
  if (!dir.exists(comp_dir)) {
    dir.create(comp_dir, recursive = TRUE)
  }
  
  # Save comparison-specific results
  write.csv(res_df, file.path(comp_dir, paste0("differential_expression_", comp_name, ".csv")), row.names = FALSE)
  
  # Volcano plot for this comparison
  volcano_data <- res_df %>%
    mutate(significant = FDR < 0.05 & abs(logFC) > 1)
  
  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_text_repel(data = subset(volcano_data, FDR < 0.05 & abs(logFC) > 2),
                    aes(label = gene), max.overlaps = 15) +
    theme_bw() +
    labs(x = "log2 Fold Change", y = "-log10(p-value)",
         title = paste("Volcano Plot:", comp_name)) +
    theme(legend.position = "none")
  
  ggsave(file.path(comp_dir, paste0("volcano_plot_", comp_name, ".png")), volcano_plot, width = 8, height = 6, dpi = 300)
  
  # MA plot for this comparison
  ma_plot <- ggplot(res_df, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = FDR < 0.05), alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    geom_hline(yintercept = 0, color = "blue") +
    theme_bw() +
    labs(x = "Mean Expression (logCPM)", y = "log2 Fold Change",
         title = paste("MA Plot:", comp_name)) +
    theme(legend.position = "none")
  
  ggsave(file.path(comp_dir, paste0("ma_plot_", comp_name, ".png")), ma_plot, width = 8, height = 6, dpi = 300)
}

# Combine all results into a single dataframe
combined_results <- do.call(rbind, all_results)

# Save combined results
write.csv(combined_results, file.path(output_dir, "differential_expression_results.csv"), row.names = FALSE)

# Summary venn diagram
library(VennDiagram)
venn_list <- lapply(all_results, function(df) {
  sig_genes <- df %>% filter(FDR < 0.05 & abs(logFC) > 1) %>% pull(gene)
  return(sig_genes)
})
venn_plot <- venn.diagram(
  x = venn_list,
  category.names = names(venn_list),
  main = "Venn Diagram of Significant Genes",
  filename = NULL
)
png(file.path(output_dir, "venn_diagram.png"), width = 800, height = 800, res = 300)
grid.draw(venn_plot)
dev.off()

# Summary heatmap of top 20 differentially expressed genes per comparison
top_genes <- combined_results %>%
  group_by(comparison) %>%
  top_n(20, wt = -FDR) %>%
  ungroup() %>%
  distinct(gene) %>%
  pull(gene)
heatmap_data <- logCPM[top_genes, ]
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename = file.path(output_dir, "heatmap_top_genes.png"),
         width = 10, height = 10)
