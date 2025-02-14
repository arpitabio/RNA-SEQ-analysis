setwd("D:\\DESEQ2_Project")
getwd()
gene_expr <- read.csv("GSE239514_HTSeq_count_matrix.csv", row.names = 1)
gene_expr
metadata <- read.csv("Metadata.csv")
metadata
install.packages("BiocManager")
BiocManager::install("DESeq2")

remove.packages("cli")
install.packages("glue")
install.packages("cli")

rownames(gene_expr)
rownames(metadata)
colnames(gene_expr)
colnames(metadata)

rownames(metadata) <- metadata$Samples
gene_expr <- gene_expr[, metadata$Samples]
all(colnames(gene_expr) == rownames(metadata))  # Should return TRUE

all(metadata$Samples %in% colnames(gene_expr))
all(colnames(gene_expr) %in% metadata$Samples)

install.packages("EnhancedVolcano")
BiocManager::install("EnhancedVolcano")

library(DESeq2)
# Create DESeq2 dataset object
metadata$Condition <- factor(metadata$Condition)

dds <- DESeqDataSetFromMatrix(countData = gene_expr, 
                              colData = metadata, 
                              design = ~ Condition)

dds <- dds[rowSums(counts(dds)) > 10, ]  # Keep genes with counts > 1

dds <- DESeq(dds)
res <- results(dds)
summary(res)

plotMA(res)
plotMA(res, main = "DESeq2 MA Plot")

library(ggplot2)
# Create a volcano plot
volcano_data <- as.data.frame(res)
volcano_data$significant <- volcano_data$padj < 0.05
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle("Volcano Plot")

library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res), 
                x = 'log2FoldChange', y = 'pvalue', 
                pCutoff = 0.05, FCcutoff = 1)

write.csv(as.data.frame(res), file = "DESeq2_results_2.csv")

# Extract significant genes
res_sig <- res[!is.na(res$padj) & !is.na(res$log2FoldChange) & 
                 res$padj < 0.1 & abs(res$log2FoldChange) > 0, ]
dim(res_sig)  # Check the number of significant genes
head(res_sig)
res_sig
library(EnhancedVolcano)

EnhancedVolcano(res_sig,
                lab = rownames(res_sig),  # Ensure labels match the dataset
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-15, 15),  # Adjust x-axis range for better visibility
                ylim = c(0, max(-log10(res_sig$padj), na.rm = TRUE) + 1),  # Dynamic y-axis
                title = 'Differential Expression Analysis',
                subtitle = 'Volcano Plot (padj < 0.1)',
                pCutoff = 0.1,  # Adjusted p-value cutoff
                FCcutoff = 1,  # Log2 fold change cutoff
                pointSize = 2.5,  # Adjusted for clarity
                labSize = 4.0,  # Label size
                colAlpha = 0.75,  # Point transparency
                col = c("grey30", "blue", "red", "purple"),  # NS, log2FC, p-value, both
                legendLabels = c("NS", "Log2 FC", "p-value", "Significant"),
                legendPosition = 'right',
                drawConnectors = TRUE,  # Connect significant labels
                widthConnectors = 0.5)



#upregulated genes
upregulated_genes <- res_sig[res_sig$log2FoldChange > 0, ]

#Downregulated genes
downregulated_genes <- res_sig[res_sig$log2FoldChange < 0, ]

head(upregulated_genes)
head(downregulated_genes)

write.csv(upregulated_genes, file = "D:\\DESEQ2_Project\\upregulated_genes.csv", row.names = TRUE)
write.csv(downregulated_genes, file = "D:\\DESEQ2_Project\\downregulated_genes.csv", row.names = TRUE)

#Coloring the genes to differentiate
# Check if gene_colors has the same length as res_sig
length(gene_colors)  # Should be equal to nrow(res_sig)

# Create an empty color vector for all genes (initially gray)
gene_colors <- rep("gray", nrow(res_sig))
# Assign blue color to upregulated genes
gene_colors[rownames(res_sig) %in% rownames(upregulated_genes)] <- "green"
# Assign pink color to downregulated genes
gene_colors[rownames(res_sig) %in% rownames(downregulated_genes)] <- "red"
names(gene_colors) <- rownames(res_sig)

# Created the EnhancedVolcano plot of DEGs with custom labels and colors
EnhancedVolcano(res_sig, 
                lab = rownames(res_sig),  # Display gene names
                x = 'log2FoldChange', 
                y = 'padj', 
                pCutoff = 0.05, 
                FCcutoff = 1, 
                title = 'DEGs Volcano Plot',
                subtitle = 'Upregulated (Green) vs Downregulated (Red)',
                colCustom = gene_colors)  # Custom colors  # Apply custom colors


library(ggplot2)

# Select top 10 genes with the highest absolute log2FoldChange
top_genes <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE), ][1:10, ]

# Create a horizontal bar plot
ggplot(top_genes, aes(x = reorder(rownames(top_genes), log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", width = 0.6) +  # Bar plot with width adjustment
  coord_flip() +  # Flip to horizontal
  scale_fill_manual(values = c("blue", "red"), labels = c("Downregulated", "Upregulated")) + 
  theme_minimal() +
  labs(title = "Top 10 Differentially Expressed Genes",
       x = "Gene",
       y = "log2 Fold Change",
       fill = "Regulation") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")


BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)
rownames(res_sig)

valid_genes <- bitr(rownames(res_sig), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
valid_genes

# Run GO enrichment using ENTREZ IDs
go_results <- enrichGO(gene          = valid_genes$ENTREZID,  # Use ENTREZ IDs
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",  # Specify ENTREZ ID as key type
                       ont           = "BP",        # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1)

# Convert results to a data frame
go_results_df <- as.data.frame(go_results)

# Print summary
summary(go_results)

# Check if results are empty
if (nrow(go_results_df) == 0) {
  print("No enriched GO terms found. Try relaxing the p-value cutoff.")
} else {
  head(go_results_df)  # View top results
}


# Visualize results
barplot(go_results, showCategory = 10)

dotplot(go_results)


#Pathway Analysis
# Run KEGG enrichment using ENTREZ IDs
kegg_results <- enrichKEGG(gene = valid_genes$ENTREZID, 
                           organism = "hsa", 
                           pvalueCutoff = 0.05)

# View the results
summary(kegg_results)
kegg_df <- as.data.frame(kegg_results)

barplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")

# Install pathview if not installed
if (!requireNamespace("pathview", quietly = TRUE)) {
  install.packages("pathview")
}
library(pathview)

# Choose a KEGG pathway to visualize (use the first significant one)
kegg_pathway_id <- kegg_df$ID[1]  # Select first enriched pathway

# Run pathview to visualize gene expression in the pathway
pathview(gene.data = valid_genes$ENTREZID, 
         pathway.id = kegg_pathway_id, 
         species = "hsa",  # Human (Homo sapiens)
         out.suffix = "KEGG_Pathway")


write.csv(kegg_df, file = "D:\\DESEQ2_Project\\kegg_enrichment_results.csv", row.names = TRUE)  # Includes row names












