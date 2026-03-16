#Step3-1  TCGA data
rm(list = ls())
gc()
getwd()

library(dplyr)
library(TCGAbiolinks)

setwd("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA")

# 读取sample_sheet和临床数据
sample_sheet <- read_tsv("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/gdc_sample_sheet.2025-09-01.tsv")
clinical_data <- read_tsv("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/clinical.tsv")
head(sample_sheet)
head(clinical_data)
print(colnames(sample_sheet))
print(colnames(clinical_data))

# 使用列名提取样本barcode
sample_barcodes <- sample_sheet$`Sample ID`

# 下载数据
batch_size <- 20
batches <- split(sample_barcodes, 
                 ceiling(seq_along(sample_barcodes) / batch_size))
if(!dir.exists("TCGA_PRAD_Batches")) {
  dir.create("TCGA_PRAD_Batches")
}

# 分批下载
for(i in seq_along(batches)) {
  cat("正在下载第", i, "批，样本数:", length(batches[[i]]), "\n")
  
  query <- GDCquery(
    project = "TCGA-PRAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    barcode = batches[[i]]
  )
  
  batch_dir <- paste0("TCGA_PRAD_Batches/batch_", i)
  GDCdownload(query, directory = batch_dir, method = "api")
  
  cat("第", i, "批下载完成\n")
  
  if(i < length(batches)) {
    Sys.sleep(5)
  }
}

# 合并所有有效数据
tsv_files <- list.files(pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
cat("找到", length(tsv_files), "个TSV文件\n")
sample_sheet$file_uuid <- sapply(strsplit(sample_sheet$`File Name`, "\\."), function(x) x[1])
file_to_sample <- setNames(sample_sheet$`Sample ID`, sample_sheet$file_uuid)

read_tcga_tsv <- function(file_path) {
  all_lines <- readLines(file_path)
  data_lines <- all_lines[!grepl("^#", all_lines) & all_lines != ""]
  con <- textConnection(data_lines)
  data <- read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  close(con)
  return(data)
}

# ID映射
test_data <- read_tcga_tsv(tsv_files[1])
gene_info <- test_data[, c("gene_id", "gene_name", "gene_type")]
counts_matrix <- matrix(0, nrow = nrow(gene_info), ncol = length(tsv_files))
rownames(counts_matrix) <- gene_info$gene_id

file_uuids <- sapply(strsplit(basename(tsv_files), "\\."), function(x) x[1])
sample_ids <- ifelse(file_uuids %in% names(file_to_sample),
                     file_to_sample[file_uuids],
                     file_uuids)

colnames(counts_matrix) <- sample_ids

#整合文件
for (i in seq_along(tsv_files)) {
  cat("处理文件", i, "/", length(tsv_files), ":", basename(tsv_files[i]), "\n")
  data <- read_tcga_tsv(tsv_files[i])
  if (!identical(data$gene_id, gene_info$gene_id)) {
    stop(paste("文件", tsv_files[i], "的基因顺序不一致"))
  }
  counts_matrix[, i] <- as.numeric(data$unstranded)
}

final_counts <- cbind(gene_info, as.data.frame(counts_matrix))
colnames <- colnames(final_counts)
correct_colnames <- gsub("\\.", "-", colnames)
correct_colnames_short <- substr(correct_colnames, 1, 12)
colnames(final_counts) <- correct_colnames_short

write.csv(final_counts, "/home/datahup/xzp/paper/PCa/data/AUcell-Model/TCGA_PRAD.csv", row.names = FALSE)
cat("基因数:", nrow(final_counts), "\n")
cat("样本数:", ncol(final_counts) - 3, "\n")


#数据处理
library(dplyr)
library(ggplot2)
library(edgeR)
library(pheatmap)
library(reshape2)
setwd("/home/datahup/xzp/paper/PCa/data/AUcell-Model")

tcga_data <- read.csv("TCGA_PRAD.csv", stringsAsFactors = FALSE)

cat("数据维度:", dim(tcga_data), "\n")
cat("前几列名:", colnames(tcga_data)[1:10], "\n")

#构建counts表达矩阵
gene_info <- tcga_data[, 1:3]  # 前3列是基因信息
gene_names <- tcga_data$gene_name
counts_matrix <- as.matrix(tcga_data[5:nrow(tcga_data), 4:ncol(tcga_data)])
rownames(counts_matrix) <- gene_names[5:length(gene_names)]
colnames(counts_matrix) <- colnames(tcga_data)[4:ncol(tcga_data)]

counts_matrix <- apply(counts_matrix, 2, as.numeric)
rownames(counts_matrix) <- gene_names[5:length(gene_names)]

cat("基因表达矩阵维度:", dim(counts_matrix), "\n")
cat("前5个行名（基因名）:", head(rownames(counts_matrix), 5), "\n")
cat("前5个列名（样本ID）:", head(colnames(counts_matrix), 5), "\n")

# 过滤低表达基因
keep <- rowSums(cpm(counts_matrix) > 1) >= 0.1 * ncol(counts_matrix)
counts_matrix_filtered <- counts_matrix[keep, ]
gene_info_filtered <- gene_info[5:nrow(gene_info), ][keep, ]
gene_names_filtered <- gene_names[5:length(gene_names)][keep]

cat("过滤后基因数量:", nrow(counts_matrix_filtered), "\n")
cat("过滤掉基因数量:", sum(!keep), "\n")
cat("保留基因比例:", round(mean(keep) * 100, 2), "%\n")

saveRDS(counts_matrix, "TCGA_PRAD_counts_matrix.rds")
gene_mapping <- gene_info[, c("gene_id", "gene_name")]
saveRDS(gene_mapping, "TCGA_gene_mapping.rds")

#构建TPM表达矩阵
library(edgeR)
counts_matrix <- readRDS("TCGA_PRAD_counts_matrix.rds")
cat("Counts矩阵维度:", dim(counts_matrix), "\n")
rpm <- t(t(counts_matrix) / colSums(counts_matrix)) * 1e6  # 计算RPM
logtpm_matrix <- log2(rpm + 1) 
saveRDS(logtpm_matrix, "TCGA_PRAD_tpm.rds")



# 读取已有的表达数据
tcga_data <- read.csv("TCGA_PRAD.csv", stringsAsFactors = FALSE)
cat("表达数据维度:", dim(tcga_data), "\n")

# 提取表达数据中的样本ID
expression_samples <- colnames(tcga_data)[4:ncol(tcga_data)]  # 跳过前3列基因信息列
cat("表达数据中的样本数量:", length(expression_samples), "\n")
cat("前10个表达数据样本ID:\n")
print(head(expression_samples, 10))

# 读取新的临床数据
clinical_data <- read_tsv("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/clinical_data.tsv")
cat("\n临床数据维度:", dim(clinical_data), "\n")
cat("临床数据列名:\n")
print(colnames(clinical_data))

# 查看临床数据中的ID列名，寻找可能的样本ID列
clinical_colnames <- colnames(clinical_data)
id_columns <- clinical_colnames[grepl("id|ID|barcode|Barcode|sample|Sample|patient|Patient", clinical_colnames, ignore.case = TRUE)]
cat("\n临床数据中可能的ID列:", id_columns, "\n")

# 简化版ID匹配
cat("=== 直接ID格式转换匹配 ===\n")

# 转换临床数据ID格式
clinical_sample_ids <- clinical_data$`Sample ID`
converted_clinical_ids <- gsub("-", ".", gsub("-01$", "", clinical_sample_ids))

cat("转换后的临床ID前10个:\n")
print(head(converted_clinical_ids, 10))
cat("表达数据ID前10个:\n")
print(head(expression_samples, 10))

# 检查匹配
matches <- converted_clinical_ids %in% expression_samples
cat("\n匹配结果:\n")
cat("匹配样本数:", sum(matches), "\n")
cat("匹配比例:", round(sum(matches)/length(clinical_sample_ids)*100, 2), "%\n")

if(sum(matches) > 0) {
  # 提取匹配的数据
  matched_clinical_ids <- clinical_sample_ids[matches]
  matched_expression_ids <- converted_clinical_ids[matches]
  
  # 表达数据
  matched_expression <- tcga_data[, c(1:3, which(colnames(tcga_data) %in% matched_expression_ids))]
  
  # 临床数据
  matched_clinical <- clinical_data[matches, ]
  
  # 保存文件
  write.csv(matched_expression, "TCGA_PRAD_matched_expression.csv", row.names = FALSE)
  write.csv(matched_clinical, "TCGA_PRAD_matched_clinical.csv", row.names = FALSE)
  
  # 创建元数据
  meta_data <- data.frame(
    Sample_ID_Expression = matched_expression_ids,
    Sample_ID_Clinical = matched_clinical_ids,
    Patient_ID = matched_clinical$`Patient ID`,
    Diagnosis_Age = matched_clinical$`Diagnosis Age`,
    T_Stage = matched_clinical$`AJCC Pathologic T-Stage`,
    Gleason_Grade = matched_clinical$`Primary Gleason Grade`,
    Survival_Status = matched_clinical$`Overall Survival Status`
  )
  
  write.csv(meta_data, "TCGA_PRAD_sample_metadata.csv", row.names = FALSE)
  
  cat("\n保存的文件:\n")
  cat("- 表达数据:", ncol(matched_expression)-3, "个样本\n")
  cat("- 临床数据:", nrow(matched_clinical), "个样本\n")
  cat("- 元数据:", nrow(meta_data), "个样本\n")
  
} else {
  cat("没有找到匹配的样本\n")
}