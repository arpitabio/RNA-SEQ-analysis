# **Differential Gene Expression Analysis using DESeq2**
## **Overview**
This repository contains an analysis of RNA-seq data using the DESeq2 package in R. The primary objective was to identify differentially expressed genes (DEGs) and perform functional enrichment analysis as part of an assignment. The workflow includes preprocessing, differential expression analysis, visualization, and pathway enrichment analysis.
## **Data Used**
- **Gene expression data**: `GSE239514_HTSeq_count_matrix.csv`
- **Metadata file**: `Metadata.csv`
## **Analysis Workflow**
### 1️⃣ Data Preparation
- Read gene expression and metadata files.
- Ensured proper formatting and alignment of metadata with gene expression data.

### 2️⃣ Differential Expression Analysis
- Conducted using `DESeq2`, filtering genes with low counts.
- Generated results including log2 fold change and adjusted p-values.

### 3️⃣ Visualization of Differentially Expressed Genes (DEGs)
- **Volcano Plots**:
  - All genes in the dataset (**padj < 0.05**).
  - DEGs with a stricter significance threshold (**padj < 0.1**).
  - Upregulated genes in **green**, downregulated genes in **red**.
- **Bar Plot**:
  - Top 10 DEGs ranked by absolute log2 fold change.
  - Downregulated genes shown in **blue**, upregulated in **red**.

### 4️⃣ Functional Enrichment Analysis
- **Gene Ontology (GO) Analysis**:
  - Focused on **Biological Processes (BP)**.
  - Bar plot displaying the **top 10 enriched GO terms**.
- **KEGG Pathway Enrichment Analysis**:
  - Identified significantly enriched pathways.
  - Bar plot showing the **top 10 pathways**.
  - Pathway visualization using `pathview`.
 ## **Files Included**
- `DESeq2_results_2.csv` → Full differential expression results.
- `upregulated_genes.csv` → List of upregulated genes (**padj < 0.1, log2FC > 0**).
- `downregulated_genes.csv` → List of downregulated genes (**padj < 0.1, log2FC < 0**).
- `kegg_enrichment_results.csv` → KEGG pathway enrichment results.

## **Visualizations Included**
📌 **Volcano Plots**:
- ✅ Volcano plot of **all genes** (padj < 0.05).
- ✅ Volcano plot of **DEGs** (padj < 0.1).
- ✅ Volcano plot of **upregulated (green) & downregulated (red) genes**.

📌 **Bar Plots**:
- ✅ **Top 10 DEGs** (blue: downregulated, red: upregulated).
- ✅ **GO enrichment analysis** (bar plot of the top 10 biological processes).
- ✅ **KEGG pathway enrichment analysis** (bar plot of the top 10 pathways).

## **Tools & Packages Used**
- **R Packages**:  
  `DESeq2`, `ggplot2`, `EnhancedVolcano`, `clusterProfiler`, `org.Hs.eg.db`, `pathview`
- **Platform**:  
  RStudio

## **Acknowledgment**
This analysis was conducted as part of an **assignment**, and the dataset used is publicly available.

