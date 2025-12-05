library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(apeglm)
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
rownames(symbol_count_matrix) <- ifelse(is.na(id_map), 
                                       flybase_ids, 
                                       id_map)

symbol_count_matrix <- symbol_count_matrix %>%
  rownames_to_column("symbol") %>%
  group_by(symbol) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("symbol")


dds <- DESeqDataSetFromMatrix(countData = symbol_count_matrix,
                              colData = sample_info,
                              design = ~ condition)

dds <- DESeq(dds)

contrast_levels <- levels(sample_info$condition)
if(length(contrast_levels) == 2) {
  res <- results(dds, contrast = c("condition", contrast_levels[2], contrast_levels[1]))
  
  resLFC <- lfcShrink(dds, coef = paste0("condition_", contrast_levels[2], "_vs_", contrast_levels[1]), type = "apeglm")
} else {
  res <- results(dds)
  resLFC <- lfcShrink(dds, type = "apeglm")
}

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(normalized_counts, file.path(output_dir, "normalized_counts.csv"))

res_df <- as.data.frame(resLFC) %>% rownames_to_column("gene")
write.csv(res_df, file.path(output_dir, "differential_expression_results.csv"), row.names = FALSE)


vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA Plot")

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

volcano_data <- as.data.frame(resLFC) %>%
  rownames_to_column("gene") %>%
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1)

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = subset(volcano_data, padj < 0.05 & abs(log2FoldChange) > 2),
                  aes(label = gene), max.overlaps = 15) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot of Differential Expression") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "volcano_plot.png"), volcano_plot, width = 8, height = 6, dpi = 300)


top_genes <- rownames(resLFC)[order(resLFC$padj)][1:20]


heatmap_data <- normalized_counts[top_genes, ]
heatmap_data <- log2(heatmap_data + 1)

annotation_df <- data.frame(Condition = sample_info$condition)
rownames(annotation_df) <- colnames(heatmap_data)

color_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(255)

heatmap_plot <- pheatmap(heatmap_data,
                         color = color_palette,
                         scale = "row",  
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         annotation_col = annotation_df,
                         show_rownames = TRUE,
                         fontsize_row = 8,
                         main = "Top 20 Differentially Expressed Genes")

png(file.path(output_dir, "heatmap_top_genes.png"), width = 8, height = 6, units = "in", res = 300)
print(heatmap_plot)
dev.off()


ma_plot <- ggplot(as.data.frame(resLFC), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_x_log10() +
  scale_color_manual(values = c("gray", "red")) +
  geom_hline(yintercept = 0, color = "blue") +
  theme_bw() +
  labs(x = "Mean Expression", y = "log2 Fold Change",
       title = "MA Plot of Differential Expression") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "ma_plot.png"), ma_plot, width = 8, height = 6, dpi = 300)