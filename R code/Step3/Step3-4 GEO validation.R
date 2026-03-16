# Step 3-4: GEO validation
rm(list = ls())
gc()

library(tidyverse)
library(pROC)
library(ggplot2)
library(caret)
library(sva)
setwd('/home/datahup/xzp/paper/PCa/data/TCGA GEO-Model')

# 1. 加载模型和选中的基因
selected_genes_df <- read.csv("seed_125_lasso_results/core_emt_genes.csv", stringsAsFactors = FALSE)
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

# 加载训练模型
final_model <- readRDS("seed_125_lasso_results/final_model.rds")
model_data <- readRDS("lasso_model_data_preprocessed.rds")

# 获取训练数据用于批次校正
x_train <- model_data$x_train
common_genes <- model_data$common_genes

# 2. 加载GEO验证数据
load("/home/datahup/xzp/paper/PCa/data/AUcell-Model/annotation_df.RData")
load("/home/datahup/xzp/paper/PCa/data/AUcell-Model/count_matrix.RData")
load("/home/datahup/xzp/paper/PCa/data/AUcell-Model/clinical_df.RData")

cat("GEO验证数据信息:\n")
cat("count_matrix维度 (samples x genes):", dim(count_matrix), "\n")
cat("clinical_df样本数:", nrow(clinical_df), "\n")
cat("count_matrix列名（基因符号）示例:\n")
print(head(colnames(count_matrix)))

# 3. 数据预处理
# 将counts转换为log2(CPM+1)
counts_per_million <- function(counts){
  lib_sizes <- rowSums(counts, na.rm = TRUE)  # 行为样本
  # 避免除以0
  lib_sizes[lib_sizes == 0] <- NA
  counts / lib_sizes * 1e6
}

cpm_mat <- counts_per_million(as.matrix(count_matrix))
logtpm_GEO <- log2(cpm_mat + 1)

cat("已将counts转换为log2(CPM+1)\n")
cat("转换后表达矩阵维度 (samples x genes):", dim(logtpm_GEO), "\n")

# 4. 提取选中的基因
# 直接使用基因符号匹配
intersect_genes <- intersect(selected_genes, colnames(logtpm_GEO))
cat("在GEO数据中找到的选中基因数量:", length(intersect_genes), "\n")
cat("找到的基因:", paste(intersect_genes, collapse = ", "), "\n")

if(length(intersect_genes) == 0) {
  cat("错误：没有找到任何选中的基因在GEO数据中\n")
  cat("GEO数据中的基因符号示例:", head(colnames(logtpm_GEO)), "\n")
  stop("基因匹配失败")
}

selected_coef_present <- selected_coef[intersect_genes]

# 提取GEO数据中的选中基因表达
X_geo_selected <- logtpm_GEO[, intersect_genes, drop = FALSE]  # samples x genes

cat("GEO选中基因表达矩阵维度 (samples x genes):", dim(X_geo_selected), "\n")

# 5. 批次校正
# 使用训练数据作为参考批次
X_train_selected <- x_train[, intersect_genes, drop = FALSE]

# 准备批次校正数据
X_train_mat <- t(X_train_selected)  # genes x TCGA样本数
X_geo_mat <- t(X_geo_selected)      # genes x GEO样本数
combined_matrix <- cbind(X_train_mat, X_geo_mat)  # genes x (TCGA + GEO)
batch <- c(rep("TCGA", ncol(X_train_mat)), rep("GEO", ncol(X_geo_mat)))

# 执行ComBat批次校正
combat_corrected <- ComBat(dat = combined_matrix, batch = batch, par.prior = TRUE, prior.plots = FALSE)
X_geo_corrected <- t(combat_corrected[, (ncol(X_train_mat) + 1):ncol(combined_matrix)])

cat("批次校正后GEO矩阵维度 (samples x genes):", dim(X_geo_corrected), "\n")

# 6. 整合临床信息
# 检查clinical_df的结构
cat("clinical_df列名:\n")
print(colnames(clinical_df))

# 根据clinical_df的实际结构创建标签
geo_clinical_processed <- clinical_df %>%
  mutate(
    sample_id = GSM,  # 假设GSM是样本ID列
    # 清理和整理T分期
    pathologic_t_clean = case_when(
      grepl("^T1", t_stage, ignore.case = TRUE) ~ "T1",
      grepl("^T2", t_stage, ignore.case = TRUE) ~ "T2", 
      grepl("^T3", t_stage, ignore.case = TRUE) ~ "T3",
      grepl("^T4", t_stage, ignore.case = TRUE) ~ "T4",
      is.na(t_stage) ~ NA_character_,
      TRUE ~ "Other"
    ),
    t_stage_label = case_when(
      pathologic_t_clean %in% c("T1", "T2") ~ 0,  # 非扩散
      pathologic_t_clean %in% c("T3", "T4") ~ 1,   # 扩散
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(t_stage_label))

cat("GEO临床信息样本数:", nrow(geo_clinical_processed), "\n")
cat("GEO T分期分布:\n")
print(table(geo_clinical_processed$t_stage_label))

# 整合表达数据和临床信息
geo_expression_data <- as.data.frame(X_geo_corrected)
geo_expression_data$sample_id <- rownames(X_geo_corrected)

geo_final_data <- geo_expression_data %>%
  inner_join(geo_clinical_processed, by = "sample_id") %>%
  filter(!is.na(t_stage_label))

cat("整合后的GEO验证数据维度:", dim(geo_final_data), "\n")
cat("最终GEO验证集T分期分布:\n")
print(table(geo_final_data$t_stage_label))
saveRDS(geo_final_data, "seed_125_lasso_results/geo_validation_data.rds")

# 7. 模型预测和评估
X_geo <- as.matrix(geo_final_data[, intersect_genes])
y_geo <- geo_final_data$t_stage_label

# 计算EMT评分
emt_score_geo <- as.numeric(X_geo %*% selected_coef_present)
geo_final_data$emt_score <- emt_score_geo

# 使用逻辑回归模型预测概率
predict_prob_geo <- predict(final_model, 
                            newdata = as.data.frame(X_geo), 
                            type = "response")

# 使用0.5作为分类阈值
predict_class_geo <- ifelse(predict_prob_geo > 0.5, 1, 0)

# 计算性能指标
geo_accuracy <- mean(predict_class_geo == y_geo)
roc_data_geo <- roc(y_geo, predict_prob_geo)
geo_auc <- auc(roc_data_geo)

# 计算其他指标
geo_sensitivity <- sum(predict_class_geo == 1 & y_geo == 1) / sum(y_geo == 1)
geo_specificity <- sum(predict_class_geo == 0 & y_geo == 0) / sum(y_geo == 0)
geo_precision <- sum(predict_class_geo == 1 & y_geo == 1) / sum(predict_class_geo == 1)

cat("\n=== GEO外部验证结果 ===\n")
cat("验证集样本数:", length(y_geo), "\n")
cat("T1/T2 (非扩散):", sum(y_geo == 0), "\n")
cat("T3/T4 (扩散):", sum(y_geo == 1), "\n")
cat("准确率:", round(geo_accuracy, 3), "\n")
cat("AUC:", round(geo_auc, 3), "\n")

# 保存评分结果
geo_scores <- data.frame(
  sample_id = geo_final_data$sample_id,
  emt_score = geo_final_data$emt_score,
  true_label = geo_final_data$t_stage_label,
  predicted_prob = predict_prob_geo,
  predicted_class = predict_class_geo,
  t_stage = geo_final_data$pathologic_t_clean
)

write.csv(geo_scores, "seed_125_lasso_results/geo_validation_scores.csv", row.names = FALSE)

# 可视化
library(pROC)
library(ggplot2)

# 加载训练集和测试集的预测结果
model_data <- readRDS("lasso_model_data_preprocessed.rds")
x_train <- model_data$x_train
x_test <- model_data$x_test
y_train <- model_data$y_train
y_test <- model_data$y_test

# 加载选中的基因和模型
selected_genes_df <- read.csv("seed_125_lasso_results/core_emt_genes.csv", stringsAsFactors = FALSE)
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

final_model <- readRDS("seed_125_lasso_results/final_model.rds")

# 训练集ROC曲线
x_train_selected <- x_train[, selected_genes, drop = FALSE]
train_pred_prob <- predict(final_model, 
                           newdata = as.data.frame(x_train_selected), 
                           type = "response")
roc_train <- roc(y_train, train_pred_prob)
auc_train <- auc(roc_train)
auc_ci_train <- ci.auc(roc_train, conf.level = 0.95)

# 准备训练集ROC数据
roc_coords_train <- coords(roc_train, "all", transpose = FALSE)
roc_data_train <- data.frame(
  Sensitivity = roc_coords_train$sensitivity,
  FPR = 1 - roc_coords_train$specificity,
  Threshold = roc_coords_train$threshold,
  Dataset = "Training Set"
)

# 生成训练集ROC图
roc_plot_train <- ggplot(roc_data_train, aes(x = FPR, y = Sensitivity)) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", 
              color = "gray50", 
              linewidth = 0.6,
              alpha = 0.7) +
  geom_line(color = "#0072B2",  # 蓝色
            linewidth = 1.2,
            alpha = 0.9) +
  annotate("text", 
           x = 0.65, y = 0.15, 
           label = paste("AUC =", round(auc_train, 3), "\n",
                         "95% CI: [", round(auc_ci_train[1], 3), ", ", 
                         round(auc_ci_train[3], 3), "]", sep = ""),
           size = 4,
           hjust = 0,
           vjust = 0,
           family = "sans") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "False Positive Rate") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "Sensitivity") +
  ggtitle("Training Set ROC Curve") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid.major = element_line(color = "gray90", 
                                    linewidth = 0.2,
                                    linetype = "solid"),
    panel.grid.minor = element_blank()
  )

# 保存训练集ROC图
pdf("seed_125_lasso_results/training_set_roc.pdf", 
    width = 5, height = 4.5, useDingbats = FALSE)
print(roc_plot_train)
dev.off()

cat("训练集AUC: ", round(auc_train, 3), "\n")
cat("训练集AUC 95% CI: [", round(auc_ci_train[1], 3), ", ", 
    round(auc_ci_train[3], 3), "]\n\n")

# 测试集ROC曲线
x_test_selected <- x_test[, selected_genes, drop = FALSE]
test_pred_prob <- predict(final_model, 
                          newdata = as.data.frame(x_test_selected), 
                          type = "response")
roc_test <- roc(y_test, test_pred_prob)
auc_test <- auc(roc_test)
auc_ci_test <- ci.auc(roc_test, conf.level = 0.95)

# 准备测试集ROC数据
roc_coords_test <- coords(roc_test, "all", transpose = FALSE)
roc_data_test <- data.frame(
  Sensitivity = roc_coords_test$sensitivity,
  FPR = 1 - roc_coords_test$specificity,
  Threshold = roc_coords_test$threshold,
  Dataset = "Test Set"
)

# 生成测试集ROC图
roc_plot_test <- ggplot(roc_data_test, aes(x = FPR, y = Sensitivity)) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", 
              color = "gray50", 
              linewidth = 0.6,
              alpha = 0.7) +
  geom_line(color = "#D55E00",  # 橙红色
            linewidth = 1.2,
            alpha = 0.9) +
  annotate("text", 
           x = 0.65, y = 0.15, 
           label = paste("AUC =", round(auc_test, 3), "\n",
                         "95% CI: [", round(auc_ci_test[1], 3), ", ", 
                         round(auc_ci_test[3], 3), "]", sep = ""),
           size = 4,
           hjust = 0,
           vjust = 0,
           family = "sans") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "False Positive Rate") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "Sensitivity") +
  ggtitle("Test Set ROC Curve") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid.major = element_line(color = "gray90", 
                                    linewidth = 0.2,
                                    linetype = "solid"),
    panel.grid.minor = element_blank()
  )

# 保存测试集ROC图
pdf("seed_125_lasso_results/test_set_roc.pdf", 
    width = 5, height = 4.5, useDingbats = FALSE)
print(roc_plot_test)
dev.off()

cat("测试集AUC: ", round(auc_test, 3), "\n")
cat("测试集AUC 95% CI: [", round(auc_ci_test[1], 3), ", ", 
    round(auc_ci_test[3], 3), "]\n\n")

# 外部验证集ROC曲线
# 加载GEO验证结果
geo_scores <- read.csv("seed_125_lasso_results/geo_validation_scores.csv")
roc_geo <- roc(geo_scores$true_label, geo_scores$predicted_prob)
auc_geo <- auc(roc_geo)
auc_ci_geo <- ci.auc(roc_geo, conf.level = 0.95)

# 准备外部验证集ROC数据
roc_coords_geo <- coords(roc_geo, "all", transpose = FALSE)
roc_data_geo <- data.frame(
  Sensitivity = roc_coords_geo$sensitivity,
  FPR = 1 - roc_coords_geo$specificity,
  Threshold = roc_coords_geo$threshold,
  Dataset = "External Validation Set"
)

# 生成外部验证集ROC图
roc_plot_geo <- ggplot(roc_data_geo, aes(x = FPR, y = Sensitivity)) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", 
              color = "gray50", 
              linewidth = 0.6,
              alpha = 0.7) +
  geom_line(color = "#009E73",  # 绿色
            linewidth = 1.2,
            alpha = 0.9) +
  annotate("text", 
           x = 0.65, y = 0.15, 
           label = paste("AUC =", round(auc_geo, 3), "\n",
                         "95% CI: [", round(auc_ci_geo[1], 3), ", ", 
                         round(auc_ci_geo[3], 3), "]", sep = ""),
           size = 4,
           hjust = 0,
           vjust = 0,
           family = "sans") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "False Positive Rate") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0.02, 0.02),
                     name = "Sensitivity") +
  ggtitle("External Validation Set ROC Curve") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid.major = element_line(color = "gray90", 
                                    linewidth = 0.2,
                                    linetype = "solid"),
    panel.grid.minor = element_blank()
  )

# 保存外部验证集ROC图
pdf("seed_125_lasso_results/external_validation_roc.pdf", 
    width = 5, height = 4.5, useDingbats = FALSE)
print(roc_plot_geo)
dev.off()