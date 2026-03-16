# Step 3-2: 建模数据处理
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/TCGA GEO-Model')
library(survival)
library(survminer)
library(glmnet)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(sva)

# 1. 加载TCGA数据
logtpm_matrix <- readRDS("/home/datahup/xzp/paper/PCa/data/AUcell-Model/TCGA_PRAD_tpm.rds")
clinical_data <- read.delim("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/clinical.tsv")

emt_related_genes <- read.csv("/home/datahup/xzp/paper/PCa/data/AUcell/EMT_Related_Genes.csv")
significant_genes <- emt_related_genes$gene
available_genes <- significant_genes[significant_genes %in% rownames(logtpm_matrix)]
cat("可用的显著基因数量:", length(available_genes), "\n")

# 2. 定义标签：T1/T2 vs T3/T4
grading_data <- clinical_data %>%
  mutate(
    patient_id = substr(cases.submitter_id, 1, 12),
    patient_id = gsub("\\.", "-", patient_id),
    pathologic_t = diagnoses.ajcc_pathologic_t
  ) %>%
  dplyr::select(patient_id, pathologic_t) %>%
  distinct(patient_id, .keep_all = TRUE)

# 清理T分期数据
grading_data <- grading_data %>%
  mutate(
    pathologic_t_clean = case_when(
      grepl("^T1", pathologic_t) | pathologic_t %in% c("T1", "T1a", "T1b", "T1c") ~ "T1",
      grepl("^T2", pathologic_t) | pathologic_t %in% c("T2", "T2a", "T2b", "T2c") ~ "T2", 
      grepl("^T3", pathologic_t) | pathologic_t %in% c("T3", "T3a", "T3b") ~ "T3",
      grepl("^T4", pathologic_t) | pathologic_t %in% c("T4") ~ "T4",
      pathologic_t == "'--" ~ NA_character_,
      is.na(pathologic_t) ~ NA_character_,
      TRUE ~ pathologic_t
    ),
    t_stage_label = case_when(
      pathologic_t_clean %in% c("T1", "T2") ~ 0,  # 非扩散
      pathologic_t_clean %in% c("T3", "T4") ~ 1,   # 扩散
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(pathologic_t_clean))

# 提取其他临床变量
clinical_extended <- clinical_data %>%
  mutate(
    patient_id = substr(cases.submitter_id, 1, 12),
    patient_id = gsub("\\.", "-", patient_id),
    pathologic_n = diagnoses.ajcc_pathologic_n,
    pathologic_m = diagnoses.ajcc_pathologic_m,
    gleason_score = diagnoses.gleason_score
  ) %>%
  dplyr::select(patient_id, pathologic_n, pathologic_m, gleason_score) %>%
  distinct(patient_id, .keep_all = TRUE)

grading_data <- grading_data %>%
  left_join(clinical_extended, by = "patient_id") %>%
  mutate(
    # 统一Gleason评分，将'--'转换为NA
    gleason_total = as.numeric(ifelse(gleason_score == "'--", NA, gleason_score)),
    # 清理N分期和M分期
    pathologic_n_clean = ifelse(pathologic_n == "'--", NA, pathologic_n),
    pathologic_m_clean = ifelse(pathologic_m == "'--", NA, pathologic_m)
  )

# 3. 准备TCGA表达矩阵和临床信息整合
tcga_emt_expr <- logtpm_matrix[available_genes, ]
tcga_expr_df <- as.data.frame(t(tcga_emt_expr))
tcga_expr_df$patient_id <- substr(rownames(tcga_expr_df), 1, 12)
tcga_expr_df$patient_id <- gsub("\\.", "-", tcga_expr_df$patient_id)

tcga_final_data <- tcga_expr_df %>%
  inner_join(grading_data, by = "patient_id") %>%
  filter(!is.na(t_stage_label))

cat("TCGA最终数据集样本数量:", nrow(tcga_final_data), "\n")
cat("TCGA标签分布:\n")
print(table(tcga_final_data$t_stage_label))

# 4. 加载和处理GSE62116数据
# 临时切换到GSE目录读取数据，然后返回主目录
gse_dir <- '/home/datahup/xzp/paper/PCa/data/raw_data_GEO/GSE62116'
gse62116_expr <- read.csv(file.path(gse_dir, "GSE62116_expression_matrix_annotated.csv"), row.names = 1)
gse62116_clinical <- read.csv(file.path(gse_dir, "GSE62116_clinical_simple.csv"))

# 为GSE62116创建标签和统一Gleason评分，并提取N分期信息
gse62116_clinical_processed <- gse62116_clinical %>%
  mutate(
    t_stage_label = case_when(
      grepl("T1|T2", pstage) ~ 0,
      grepl("T3|T4", pstage) ~ 1,
      pstage == "TxN+" ~ 1,
      TRUE ~ NA_real_
    ),
    sample_id = geo_accession,
    dataset = "GSE62116",
    # 将pathgs统一为gleason_total
    gleason_total = as.numeric(pathgs),
    # 从pstage中提取N分期信息
    pathologic_n_clean = case_when(
      grepl("N0", pstage) ~ "N0",    # 淋巴结阴性
      grepl("N\\+", pstage) ~ "N+",  # 淋巴结阳性
      grepl("N1", pstage) ~ "N+",    # N1视为阳性
      grepl("N2", pstage) ~ "N+",    # N2视为阳性
      grepl("N3", pstage) ~ "N+",    # N3视为阳性
      pstage == "TxN+" ~ "N+",       # 明确标注的淋巴结阳性
      TRUE ~ NA_character_           # 其他情况为NA
    ),
    # M分期信息（GSE62116通常没有M分期）
    pathologic_m_clean = NA
  ) %>%
  filter(!is.na(t_stage_label))

# 检查提取的N分期分布
cat("GSE62116 N分期分布:\n")
print(table(gse62116_clinical_processed$pathologic_n_clean, useNA = "always"))

# 去除重复样本
gse62116_clinical_dedup <- gse62116_clinical_processed[!duplicated(gse62116_clinical_processed$sample_id), ]

# 准备GSE62116表达矩阵
gse62116_available_genes <- available_genes[available_genes %in% rownames(gse62116_expr)]
gse62116_emt_expr <- gse62116_expr[gse62116_available_genes, ]
gse62116_expr_df <- as.data.frame(t(gse62116_emt_expr))
gse62116_expr_df$sample_id <- rownames(gse62116_expr_df)

gse62116_final_data <- gse62116_expr_df %>%
  inner_join(gse62116_clinical_dedup, by = "sample_id") %>%
  filter(!is.na(t_stage_label))

cat("GSE62116最终数据集样本数量:", nrow(gse62116_final_data), "\n")
cat("GSE62116标签分布:\n")
print(table(gse62116_final_data$t_stage_label))

# 5. 合并两个数据集
# 找到共同基因
clinical_cols <- c("patient_id", "sample_id", "pathologic_t", "pathologic_t_clean", 
                   "pathologic_n", "pathologic_m", "gleason_score", "t_stage_label",
                   "dataset", "geo_accession", "pathgs", "preop_psa", "age", "pstage",
                   "pathologic_n_clean", "pathologic_m_clean", "gleason_total", "data_source")

common_genes <- intersect(colnames(tcga_final_data), colnames(gse62116_final_data))
common_genes <- common_genes_corrected[!common_genes_corrected %in% clinical_cols]

cat("修正后的共同基因数量:", length(common_genes), "\n")
cat("修正后的共同基因:\n")
print(common_genes)

# 准备合并的数据 - 包含统一的临床信息
tcga_for_merge <- tcga_final_data %>%
  mutate(
    data_source = "TCGA",
    sample_id = patient_id
  ) %>%
  dplyr::select(all_of(common_genes), t_stage_label, gleason_total, pathologic_n_clean, pathologic_m_clean, data_source, sample_id)

gse62116_for_merge <- gse62116_final_data %>%
  mutate(
    data_source = "GSE62116"
  ) %>%
  dplyr::select(all_of(common_genes), t_stage_label, gleason_total, pathologic_n_clean, pathologic_m_clean, data_source, sample_id)

combined_data <- bind_rows(tcga_for_merge, gse62116_for_merge)

cat("合并后总样本数量:", nrow(combined_data), "\n")
cat("合并后标签分布:\n")
print(table(combined_data$t_stage_label, combined_data$data_source))

# 检查临床变量分布
cat("Gleason评分统计:\n")
gleason_summary <- combined_data %>%
  group_by(data_source) %>%
  summarise(
    n = n(),
    mean_gleason = mean(gleason_total, na.rm = TRUE),
    sd_gleason = sd(gleason_total, na.rm = TRUE),
    missing_gleason = sum(is.na(gleason_total))
  )
print(gleason_summary)

cat("N分期分布:\n")
n_stage_summary <- combined_data %>%
  group_by(data_source, pathologic_n_clean) %>%
  summarise(n = n()) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))
print(n_stage_summary)

# 6. 准备lasso建模数据
x <- as.matrix(combined_data[, common_genes])
y <- combined_data$t_stage_label
batch <- combined_data$data_source

cat("预测变量维度:", dim(x), "\n")
cat("响应变量分布:\n")
print(table(y))

# 7. 批处理校正
if(length(unique(batch)) > 1) {
  # 转置矩阵用于ComBat
  x_t <- t(x)
  x_combat_t <- ComBat(dat = x_t, batch = batch)
  x_combat <- t(x_combat_t)
  cat("批处理校正完成\n")
} else {
  x_combat <- x
  cat("只有一个批次，无需校正\n")
}

# 8. 数据分割和样本匹配信息保存
set.seed(1)
train_index <- sample(1:nrow(combined_data), 0.7 * nrow(combined_data))

x_train <- x_combat[train_index, ]
x_test <- x_combat[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

# 创建详细的样本匹配信息
sample_matching_info <- data.frame(
  sample_id = combined_data$sample_id,
  data_source = combined_data$data_source,
  t_stage_label = combined_data$t_stage_label,
  gleason_total = combined_data$gleason_total,
  pathologic_n = combined_data$pathologic_n_clean,
  pathologic_m = combined_data$pathologic_m_clean,
  row_index = 1:nrow(combined_data),
  is_training = 1:nrow(combined_data) %in% train_index,
  train_set_index = ifelse(1:nrow(combined_data) %in% train_index, 
                           match(1:nrow(combined_data), train_index), NA),
  test_set_index = ifelse(!(1:nrow(combined_data) %in% train_index), 
                          match(1:nrow(combined_data), setdiff(1:nrow(combined_data), train_index)), NA)
)

cat("训练集样本数:", length(y_train), "\n")
cat("测试集样本数:", length(y_test), "\n")
cat("训练集标签分布:\n")
print(table(y_train))
cat("测试集标签分布:\n")
print(table(y_test))

# 9. 保存训练集和测试集的样本ID列表
training_samples <- sample_matching_info %>% 
  filter(is_training) %>% 
  dplyr::select(sample_id, data_source, t_stage_label, gleason_total, pathologic_n, pathologic_m)

testing_samples <- sample_matching_info %>% 
  filter(!is_training) %>% 
  dplyr::select(sample_id, data_source, t_stage_label, gleason_total, pathologic_n, pathologic_m)

cat("训练集样本信息:\n")
print(head(training_samples))
cat("测试集样本信息:\n")
print(head(testing_samples))

# 10. 保存数据到当前目录('/home/datahup/xzp/paper/PCa/data/TCGA GEO-Model')
model_data <- list(
  x_train = x_train,
  x_test = x_test,
  y_train = y_train,
  y_test = y_test,
  combined_data = combined_data,
  common_genes = common_genes,
  batch_info = batch,
  sample_matching_info = sample_matching_info,
  training_samples = training_samples,
  testing_samples = testing_samples,
  train_indices = train_index,
  test_indices = setdiff(1:nrow(combined_data), train_index)
)

saveRDS(model_data, file = "lasso_model_data_preprocessed.rds")

# 单独保存样本匹配信息为CSV文件，方便查看
write.csv(sample_matching_info, "sample_matching_info.csv", row.names = FALSE)
write.csv(training_samples, "training_samples.csv", row.names = FALSE)
write.csv(testing_samples, "testing_samples.csv", row.names = FALSE)

# 保存合并后的完整数据
write.csv(combined_data, "combined_data_with_clinical_info.csv", row.names = FALSE)

# 保存基因列表
write.csv(data.frame(gene = common_genes), "emt_genes_list.csv", row.names = FALSE)
