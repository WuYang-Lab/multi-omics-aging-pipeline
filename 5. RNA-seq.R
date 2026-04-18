library(DESeq2)
library(data.table)
setwd("D:\\1.work\\1.Aging\\GTEx_DEG\\gene_read_counts\\")
file_list <- list.files(pattern = "*.gct")
Age <- read.table("GTEx_Age_diff.txt", header = TRUE)

## format gene expression matrix
# 1. update colnames
for (file_name in file_list) {
  gene_data <- fread(file_name, header = TRUE)
  new_colnames <- sapply(colnames(gene_data), function(x) {
    parts <- strsplit(x, "-")[[1]]
    paste(parts[1:2], collapse = "-")
  })
  colnames(gene_data) <- new_colnames
  
# 2. Processing of gene symbols
  gene_data$`Description-NA` <- as.character(gene_data$`Description-NA`)
  gene_data$`Description-NA` <- make.unique(gene_data$`Description-NA`)
  gene_data <- gene_data[, -c(1, 2)]
  write.table(gene_data, file = "temp_gene.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  gene_data <- fread("temp_gene.txt", header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  rownames(gene_data) <- gene_data[, 1]
  gene_data <- gene_data[, -1]
  
# 3. Add case control group
  gene_data_colnames <- colnames(gene_data)
  Age_filtered <- Age[Age$SUBJID %in% gene_data_colnames, ]
  Age_filtered$Group <- ifelse(Age_filtered$AGE <= 25, 'control', 'case')
  rownames(Age_filtered) <- Age_filtered$SUBJID
  Age_filtered$Group <- factor(Age_filtered$Group, ordered = FALSE)
  
## Diff analysis  
  new_sample_order <- Age_filtered$SUBJID 
  gene_data <- gene_data[, new_sample_order] 
  Age_filtered$Group <- relevel(Age_filtered$Group, ref = "control") 
  dds <- DESeqDataSetFromMatrix(gene_data, colData = Age_filtered, design = ~ Group)
  dds_DE <- DESeq(dds)
  res_DE <- results(dds_DE)
  gene_info <- res_DE[, c("log2FoldChange", "pvalue")]
  output_file_name <- paste0(sub(".gct", "", file_name), "_results.txt")
  write.table(gene_info, file = output_file_name, sep = "\t", quote = FALSE, col.names = NA)
}
