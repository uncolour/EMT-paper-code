# Step 3-3: LASSO建模
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/TCGA GEO-Model')
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)
library(caret)

# 加载预处理的数据
model_data <- readRDS("lasso_model_data_preprocessed.rds")

# 提取数据
x_train <- model_data$x_train
x_test <- model_data$x_test
y_train <- model_data$y_train
y_test <- model_data$y_test
combined_data <- model_data$combined_data
common_genes <- model_data$common_genes
batch_info <- model_data$batch_info

cat("训练集维度:", dim(x_train), "\n")
cat("测试集维度:", dim(x_test), "\n")
cat("训练集标签分布:\n")
print(table(y_train))
cat("测试集标签分布:\n")
print(table(y_test))

# 创建结果目录
result_dir <- "seed_125_lasso_results"
if(!dir.exists(result_dir)) {
  dir.create(result_dir)
}

# 1. LASSO特征选择
set.seed(2)
cv.lasso <- cv.glmnet(x_train, y_train, 
                      alpha = 1,           # LASSO回归
                      family = "binomial", # 二分类
                      nfolds = 10,         # 10折交叉验证
                      type.measure = "deviance",  # 使用deviance而不是class
                      standardize = TRUE,   # 标准化特征
                      nlambda = 100)        # 增加lambda网格密度

png(file.path(result_dir, "lasso_cv_plot.png"), width = 800, height = 600)
plot(cv.lasso, main = "LASSO Cross-Validation")
abline(v = log(cv.lasso$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = log(cv.lasso$lambda.1se), col = "blue", lty = 2, lwd = 2)
legend("topright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, lwd = 2, cex = 0.8)
dev.off()

cat("最优lambda (min):", cv.lasso$lambda.min, "\n")
cat("最优lambda (1se):", cv.lasso$lambda.1se, "\n")

# 使用lambda.min（更宽松的标准）
lasso_model_min <- glmnet(x_train, y_train, 
                          alpha = 1, 
                          family = "binomial",
                          lambda = cv.lasso$lambda.min)

# 使用lambda.1se（更严格的标准）
lasso_model_1se <- glmnet(x_train, y_train, 
                          alpha = 1, 
                          family = "binomial",
                          lambda = cv.lasso$lambda.1se)

# 使用介于min和1se之间的lambda值
lambda_custom <- (cv.lasso$lambda.min + cv.lasso$lambda.1se) / 2
lasso_model_custom <- glmnet(x_train, y_train, 
                             alpha = 1, 
                             family = "binomial",
                             lambda = lambda_custom)

# 比较三种方法选择的基因数量
extract_genes <- function(model, threshold = 0) {
  coef_matrix <- as.matrix(coef(model))
  selected_indices <- which(abs(coef_matrix) > threshold)
  selected_genes <- rownames(coef_matrix)[selected_indices]
  selected_genes <- selected_genes[selected_genes != "(Intercept)"]
  return(selected_genes)
}

genes_min <- extract_genes(lasso_model_min)
genes_1se <- extract_genes(lasso_model_1se)
genes_custom <- extract_genes(lasso_model_custom)

cat("方法1 (lambda.min) 选择基因数:", length(genes_min), "\n")
cat("方法2 (lambda.1se) 选择基因数:", length(genes_1se), "\n")
cat("方法3 (custom lambda) 选择基因数:", length(genes_custom), "\n")

# 选择lambda.1se（更严格的标准）
lasso_model <- lasso_model_1se
selected_genes <- genes_1se

# 提取选中的基因和系数
lasso_coef <- as.matrix(coef(lasso_model))
selected_indices <- which(lasso_coef != 0)
selected_genes <- rownames(lasso_coef)[selected_indices][-1] # 去除截距项
selected_coef <- lasso_coef[selected_indices][-1]

cat("最终选中的基因数量:", length(selected_genes), "\n")
cat("选中的基因及系数:\n")
for(i in 1:length(selected_genes)) {
  cat(selected_genes[i], ":", round(selected_coef[i], 4), "\n")
}

# 保存选中的基因
selected_genes_df <- data.frame(
  gene = selected_genes,
  coefficient = selected_coef,
  absolute_coef = abs(selected_coef),
  direction = ifelse(selected_coef > 0, "Risk", "Protective")
) %>%
  arrange(desc(absolute_coef))

write.csv(selected_genes_df, file.path(result_dir, "core_emt_genes.csv"), row.names = FALSE)

# LASSO特征选择结果可视化
# 创建CV图PDF（简洁无标题版）
pdf(file.path(result_dir, "lasso_cv_plot.pdf"), width = 8, height = 6)
par(mar = c(5, 5, 2, 2))

plot(cv.lasso, 
     xlab = expression(log(lambda)), 
     ylab = "Binomial Deviance",
     cex.lab = 1.3,
     cex.axis = 1.1)

# 添加lambda.1se垂直线
abline(v = log(cv.lasso$lambda.1se), col = "blue", lty = 2, lwd = 2)
dev.off()

# 创建系数路径图PDF
pdf(file.path(result_dir, "lasso_coefficient_path_no_labels.pdf"), width = 14, height = 8)

# 设置图形参数
par(mar = c(5, 5, 2, 10))

# 获取glmnet模型对象
lasso_full <- cv.lasso$glmnet.fit

# 获取所有lambda值和系数
lambda_seq <- lasso_full$lambda
coef_matrix <- lasso_full$beta
gene_names <- rownames(coef_matrix)

# 反转lambda顺序：从大到小 -> 从小到大
lambda_seq_rev <- rev(lambda_seq)
coef_matrix_rev <- coef_matrix[, ncol(coef_matrix):1]

# 创建空图框架
x_vals <- -log(lambda_seq_rev)
x_range <- range(x_vals)
y_range <- range(coef_matrix_rev)

# 扩展y轴范围
y_padding <- diff(y_range) * 0.15
plot_ylim <- c(y_range[1] - y_padding, y_range[2] + y_padding)

# 绘制基础图形
plot(x_range, plot_ylim, 
     type = "n",
     xlab = "-log(lambda)",  
     ylab = "Coefficient",
     cex.lab = 1.3,
     cex.axis = 1.1)

# 绘制所有特征的系数路径（灰色细线）
for(i in 1:nrow(coef_matrix_rev)) {
  lines(x_vals, coef_matrix_rev[i,], 
        col = adjustcolor("gray60", alpha = 0.4), 
        lwd = 0.8)
}

# 找到选中的基因索引
selected_indices <- match(selected_genes, gene_names)

# 根据您观察的线条颜色顺序定义颜色（从上到下：绿、红、蓝、橙、紫）
cat("\n=== 调试信息 ===\n")

# 在λ.1se处获取选中的基因系数值，查看实际顺序
lambda_1se_idx <- which.min(abs(lambda_seq_rev - cv.lasso$lambda.1se))
coef_at_1se <- coef_matrix_rev[selected_indices, lambda_1se_idx]

# 创建一个数据框显示实际系数值和顺序
actual_order <- data.frame(
  gene = selected_genes,
  coef_at_1se = coef_at_1se,
  index = selected_indices,
  stringsAsFactors = FALSE
)

# 按系数值从大到小排序（从上到下）
actual_order <- actual_order[order(actual_order$coef_at_1se, decreasing = TRUE), ]
print(actual_order)

# 根据实际观察到的顺序分配颜色
color_palette <- c("#4DAF4A", "#E41A1C", "#377EB8", "#FF7F00", "#984EA3")
gene_colors <- setNames(color_palette, actual_order$gene)

cat("\n颜色分配:\n")
for(i in 1:nrow(actual_order)) {
  cat(actual_order$gene[i], "->", gene_colors[actual_order$gene[i]], "(位置", i, ")\n")
}

# 突出显示选中的5个基因 - 使用按实际顺序分配的颜色
for(j in 1:length(selected_indices)) {
  i <- selected_indices[j]
  gene_name <- selected_genes[j]
  color_val <- gene_colors[gene_name]
  
  if(!is.na(i)) {
    lines(x_vals, coef_matrix_rev[i,], 
          col = color_val, 
          lwd = 3)
    cat("已绘制线条:", gene_name, "颜色:", color_val, "\n")
  }
}

# 添加lambda.1se垂直线
lambda_1se_x <- -log(cv.lasso$lambda.1se)
abline(v = lambda_1se_x, col = "blue", lty = 2, lwd = 2.5)

# 获取线条尾部（λ最小处）的系数值
lambda_min_idx <- length(lambda_seq_rev)
coef_at_tail <- coef_matrix_rev[selected_indices, lambda_min_idx]

# 创建完整的基因信息（保持颜色分配）
gene_info <- data.frame(
  name = selected_genes,
  coef_1se = coef_at_1se,
  coef_tail = coef_at_tail,
  color = gene_colors[selected_genes],
  stringsAsFactors = FALSE
)

cat("\n基因完整信息:\n")
print(gene_info)

# 在λ.1se处添加点
for(j in 1:nrow(gene_info)) {
  gene_name <- gene_info$name[j]
  color_val <- gene_info$color[j]
  
  points(x = lambda_1se_x,
         y = gene_info$coef_1se[j],
         col = color_val, 
         pch = 19, 
         cex = 1.8,
         lwd = 2)
}
dev.off()

# 2. 构建最终预测模型
# 使用选中的基因构建模型
x_train_selected <- x_train[, selected_genes, drop = FALSE]
x_test_selected <- x_test[, selected_genes, drop = FALSE]

# 训练逻辑回归模型
final_model <- glm(y_train ~ ., 
                   data = as.data.frame(x_train_selected), 
                   family = binomial())
cat("最终模型摘要:\n")
print(summary(final_model))

# 模型预测
train_pred_prob <- predict(final_model, as.data.frame(x_train_selected), type = "response")
test_pred_prob <- predict(final_model, as.data.frame(x_test_selected), type = "response")

# 计算EMT评分
emt_score_train <- as.numeric(x_train_selected %*% selected_coef)
emt_score_test <- as.numeric(x_test_selected %*% selected_coef)
x_all <- rbind(x_train, x_test)
emt_score_all <- as.numeric(x_all[, selected_genes, drop = FALSE] %*% selected_coef)

# 将评分添加到combined_data中
combined_data$emt_score <- emt_score_all

# 3. 模型评估
# 训练集性能
train_pred_class <- ifelse(train_pred_prob > 0.5, 1, 0)
train_accuracy <- mean(train_pred_class == y_train)
train_auc <- auc(roc(y_train, train_pred_prob))

# 测试集性能
test_pred_class <- ifelse(test_pred_prob > 0.5, 1, 0)
test_accuracy <- mean(test_pred_class == y_test)
test_auc <- auc(roc(y_test, test_pred_prob))

cat("\n=== 模型性能评估 ===\n")
cat("训练集 - 准确率:", round(train_accuracy, 3), " AUC:", round(train_auc, 3), "\n")
cat("测试集 - 准确率:", round(test_accuracy, 3), " AUC:", round(test_auc, 3), "\n")

# 4. 可视化结果
# ROC曲线
roc_data <- roc(y_test, test_pred_prob)
png(file.path(result_dir, "roc_curve.png"), width = 600, height = 600)
plot(roc_data, main = paste("ROC Curve (AUC =", round(test_auc, 3), ")"), 
     col = "blue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "red")
dev.off()

# 特征重要性图
importance_df <- selected_genes_df %>%
  arrange(absolute_coef)

p_importance <- ggplot(importance_df, aes(x = reorder(gene, absolute_coef), y = absolute_coef, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Core EMT Genes Importance",
       x = "Gene", y = "Coefficient Magnitude",
       subtitle = paste("Selected by LASSO (lambda.1se)")) +
  theme_minimal() +
  scale_fill_manual(values = c("Risk" = "red", "Protective" = "blue")) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(file.path(result_dir, "feature_importance.png"), p_importance, width = 10, height = 8, dpi = 300)

# 5. 保存模型和结果
# 保存模型
saveRDS(lasso_model, file.path(result_dir, "lasso_model.rds"))
saveRDS(final_model, file.path(result_dir, "final_model.rds"))
saveRDS(combined_data, file.path(result_dir, "final_data_with_scores.rds"))

# 保存评分
emt_scores <- data.frame(
  sample_id = combined_data$sample_id,
  data_source = combined_data$data_source,
  emt_score = combined_data$emt_score,
  t_stage_label = combined_data$t_stage_label,
  gleason_total = combined_data$gleason_total,
  pathologic_n = combined_data$pathologic_n_clean
)

write.csv(emt_scores, file.path(result_dir, "emt_scores.csv"), row.names = FALSE)

# 保存评估结果
evaluation_results <- list(
  selected_genes = selected_genes,
  coefficients = selected_coef,
  train_accuracy = train_accuracy,
  train_auc = train_auc,
  test_accuracy = test_accuracy,
  test_auc = test_auc,
  lambda_min = cv.lasso$lambda.min,
  lambda_1se = cv.lasso$lambda.1se,
  lambda_custom = lambda_custom
)

saveRDS(evaluation_results, file.path(result_dir, "model_evaluation_results.rds"))
