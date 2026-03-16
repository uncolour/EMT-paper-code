# Step 3-5: Model Functional Evaluation
rm(list = ls())
gc()

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# 设置工作目录和输出路径
setwd('/home/datahup/xzp/paper/PCa/data/TCGA GEO-Model')
output_dir <- '/home/datahup/xzp/paper/PCa/fig/TCGA GEO Model'

# 1. 加载所有数据
# 加载模型数据
model_data <- readRDS("lasso_model_data_preprocessed.rds")
combined_data <- model_data$combined_data

selected_genes_df <- read.csv("seed_125_lasso_results/core_emt_genes.csv", stringsAsFactors = FALSE)
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

geo_final_data <- readRDS("seed_125_lasso_results/geo_validation_data.rds")

# 2. 计算各数据集的EMT评分
x_train <- model_data$x_train
x_test <- model_data$x_test
y_train <- model_data$y_train
y_test <- model_data$y_test

# 训练集EMT评分
x_train_selected <- x_train[, selected_genes, drop = FALSE]
emt_score_train <- as.numeric(x_train_selected %*% selected_coef)
train_scores <- data.frame(
  sample_id = model_data$training_samples$sample_id,
  emt_score = emt_score_train,
  label = y_train,
  dataset = "Training"
)

# 测试集EMT评分
x_test_selected <- x_test[, selected_genes, drop = FALSE]
emt_score_test <- as.numeric(x_test_selected %*% selected_coef)
test_scores <- data.frame(
  sample_id = model_data$testing_samples$sample_id,
  emt_score = emt_score_test,
  label = y_test,
  dataset = "Testing"
)

# 验证集EMT评分
X_geo_selected <- as.matrix(geo_final_data[, selected_genes])
emt_score_geo <- as.numeric(X_geo_selected %*% selected_coef)
geo_scores <- data.frame(
  sample_id = geo_final_data$sample_id,
  emt_score = emt_score_geo,
  label = geo_final_data$t_stage_label,
  dataset = "Validation"
)

# 3. 整合临床信息 - 直接使用label作为T分期分组
train_indices <- model_data$train_indices
test_indices <- model_data$test_indices

# 训练集添加临床信息
train_scores_with_clinical <- train_scores %>%
  mutate(
    sample_id_tcga = sample_id,
    t_stage_group = ifelse(label == 0, "T1/T2 (Non-invasive)", "T3/T4 (Invasive)"),
    gleason_total = as.numeric(as.character(combined_data$gleason_total[train_indices])),
    dataset_type = "Training"
  )

# 测试集添加临床信息
test_scores_with_clinical <- test_scores %>%
  mutate(
    sample_id_tcga = sample_id,
    t_stage_group = ifelse(label == 0, "T1/T2 (Non-invasive)", "T3/T4 (Invasive)"),
    gleason_total = as.numeric(as.character(combined_data$gleason_total[test_indices])),
    dataset_type = "Testing"
  )

# 验证集添加临床信息
geo_scores_with_clinical <- geo_scores %>%
  mutate(
    sample_id_tcga = sample_id,
    t_stage_group = ifelse(label == 0, "T1/T2 (Non-invasive)", "T3/T4 (Invasive)"),
    gleason_total = sapply(geo_final_data$gleason_total, function(x) {
      if (is.na(x)) return(NA)
      if (grepl("=", x)) {
        total <- as.numeric(strsplit(x, "=")[[1]][1])
        return(total)
      } else {
        return(as.numeric(x))
      }
    }),
    dataset_type = "Validation"
  )

# 4. 数据清理和标准化
prepare_dataset <- function(data) {
  data_clean <- data %>%
    filter(!is.na(gleason_total)) %>%
    mutate(
      gleason_group = ifelse(gleason_total <= 7, "Gleason ≤7 (Low/Intermediate)", "Gleason ≥8 (High)")
    ) %>%
    filter(gleason_group %in% c("Gleason ≤7 (Low/Intermediate)", "Gleason ≥8 (High)"))
  
  return(data_clean)
}

# 处理各数据集
train_clean <- prepare_dataset(train_scores_with_clinical)
test_clean <- prepare_dataset(test_scores_with_clinical)
validation_clean <- prepare_dataset(geo_scores_with_clinical)

cat("训练集样本数:", nrow(train_clean), "\n")
cat("测试集样本数:", nrow(test_clean), "\n")
cat("验证集样本数:", nrow(validation_clean), "\n")

# 5. 生成箱线图
create_boxplot <- function(data, dataset_name, var_type, fill_colors) {
  if (var_type == "t_stage") {
    x_var <- "t_stage_group"
    title <- paste("EMT Score by T Stage -", dataset_name)
  } else {
    x_var <- "gleason_group"
    title <- paste("EMT Score by Gleason Grade -", dataset_name)
  }
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = emt_score, fill = .data[[x_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.2, aes(color = .data[[x_var]])) +
    labs(
      title = title,
      x = if(var_type == "t_stage") "T Stage Group" else "Gleason Grade",
      y = "EMT Score"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors)
  
  return(p)
}

# 定义颜色
t_stage_colors <- c("#2E8B57", "#CD5C5C")
gleason_colors <- c("#4682B4", "#FF6347")

# 生成各数据集的箱线图
p_train_t <- create_boxplot(train_clean, "Training Set", "t_stage", t_stage_colors)
p_test_t <- create_boxplot(test_clean, "Testing Set", "t_stage", t_stage_colors)
p_validation_t <- create_boxplot(validation_clean, "Validation Set", "t_stage", t_stage_colors)

p_train_g <- create_boxplot(train_clean, "Training Set", "gleason", gleason_colors)
p_test_g <- create_boxplot(test_clean, "Testing Set", "gleason", gleason_colors)
p_validation_g <- create_boxplot(validation_clean, "Validation Set", "gleason", gleason_colors)

# 6. 组合和保存图表
# T分期组合图
t_stage_combined <- grid.arrange(
  p_train_t, p_test_t, p_validation_t,
  nrow = 1,
  top = textGrob("EMT Score Distribution by T Stage Groups", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Gleason评分组合图
gleason_combined <- grid.arrange(
  p_train_g, p_test_g, p_validation_g,
  nrow = 1,
  top = textGrob("EMT Score Distribution by Gleason Grade Groups", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# 保存组合图
ggsave(file.path(output_dir, "EMT_TStage_By_Dataset.pdf"), t_stage_combined, 
       width = 18, height = 6, dpi = 300)
ggsave(file.path(output_dir, "EMT_Gleason_By_Dataset.pdf"), gleason_combined, 
       width = 18, height = 6, dpi = 300)

# 保存单个图表
ggsave(file.path(output_dir, "EMT_TStage_Training.pdf"), p_train_t, width = 6, height = 6)
ggsave(file.path(output_dir, "EMT_TStage_Testing.pdf"), p_test_t, width = 6, height = 6)
ggsave(file.path(output_dir, "EMT_TStage_Validation.pdf"), p_validation_t, width = 6, height = 6)

ggsave(file.path(output_dir, "EMT_Gleason_Training.pdf"), p_train_g, width = 6, height = 6)
ggsave(file.path(output_dir, "EMT_Gleason_Testing.pdf"), p_test_g, width = 6, height = 6)
ggsave(file.path(output_dir, "EMT_Gleason_Validation.pdf"), p_validation_g, width = 6, height = 6)

# 7. 生存分析
library(survival)
library(survminer)
library(timeROC)

# 加载生存数据
load("/home/datahup/xzp/paper/PCa/data/AUcell-Model/survival_df.RData")

# 为验证集数据添加生存信息
validation_survival_data <- validation_clean %>%
  left_join(survival_df, by = c("sample_id_tcga" = "GSM")) %>%
  filter(!is.na(event) & !is.na(time))

# 检查生存数据
cat("验证集生存数据样本数:", nrow(validation_survival_data), "\n")
cat("BCR事件分布:\n")
event_table <- table(validation_survival_data$event)
print(event_table)
cat("事件率:", round(event_table[2]/sum(event_table)*100, 1), "%\n")

# 准备生存分析数据
validation_survival_data <- validation_survival_data %>%
  mutate(
    bcr_event = event,
    survival_time = time,
    
    # 使用中位数分组
    emt_group = ifelse(emt_score > median(emt_score, na.rm = TRUE), 
                       "High_EMT", "Low_EMT"),
    
    # 创建三分位分组
    emt_tertile = cut(emt_score, 
                      breaks = quantile(emt_score, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Low", "Intermediate", "High"),
                      include.lowest = TRUE)
  )

# 生存数据汇总
cat("\n=== 生存分析数据集汇总 ===\n")
cat("总样本数:", nrow(validation_survival_data), "\n")
cat("EMT分组分布:\n")
print(table(validation_survival_data$emt_group))
cat("EMT三分位分布:\n")
print(table(validation_survival_data$emt_tertile))

cat("随访时间统计:\n")
cat("中位随访时间:", round(median(validation_survival_data$survival_time), 1), "月\n")
cat("随访时间范围:", round(range(validation_survival_data$survival_time), 1), "月\n")

# 8. Kaplan-Meier生存分析
# 按中位数分组
km_fit_median <- survfit(Surv(survival_time, bcr_event) ~ emt_group, 
                         data = validation_survival_data)

# 按三分位分组
km_fit_tertile <- survfit(Surv(survival_time, bcr_event) ~ emt_tertile, 
                          data = validation_survival_data)

# 绘制中位数分组KM曲线
p_km_median <- ggsurvplot(
  km_fit_median,
  data = validation_survival_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  surv.median.line = "hv",
  palette = c("#E41A1C", "#377EB8"),
  legend.labs = c("High EMT", "Low EMT"),
  legend.title = "EMT Group",
  xlab = "Time to BCR (Months)",
  ylab = "BCR-Free Survival Probability",
  ggtheme = theme_classic()
)

# 绘制三分位分组KM曲线
p_km_tertile <- ggsurvplot(
  km_fit_tertile,
  data = validation_survival_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
  legend.labs = c("Low EMT", "Intermediate EMT", "High EMT"),
  legend.title = "EMT Tertile",
  xlab = "Time to BCR (Months)",
  ylab = "BCR-Free Survival Probability",
  ggtheme = theme_classic()
)

# 保存KM曲线
ggsave(file.path(output_dir, "BCR_KM_Curve_EMT_Median.pdf"), 
       p_km_median$plot, width = 8, height = 8)
ggsave(file.path(output_dir, "BCR_KM_Curve_EMT_Tertile.pdf"), 
       p_km_tertile$plot, width = 8, height = 8)
cat("Kaplan-Meier曲线已保存\n")

# 9. 核心基因的Cox回归分析
selected_genes_df <- read.csv("seed_125_lasso_results/core_emt_genes.csv", stringsAsFactors = FALSE)

# 获取五个核心基因
core_genes <- selected_genes_df$gene
cat("分析的五个核心基因:\n")
print(core_genes)

# 准备基因表达数据
# 从geo_final_data中提取基因表达数据
geo_final_data <- readRDS("seed_125_lasso_results/geo_validation_data.rds")

# 确保基因表达数据存在
gene_expression_data <- geo_final_data %>%
  select(sample_id, all_of(core_genes)) %>%
  # 标准化基因表达数据（z-score）
  mutate(across(all_of(core_genes), ~scale(.)[,1]))

# 合并生存数据
gene_survival_data <- validation_survival_data %>%
  left_join(gene_expression_data, by = c("sample_id_tcga" = "sample_id")) %>%
  # 清理数据
  filter(!is.na(survival_time) & !is.na(bcr_event)) %>%
  # 确保所有基因表达数据都不缺失
  filter(if_all(all_of(core_genes), ~!is.na(.)))

cat(sprintf("用于基因分析的样本数: %d\n", nrow(gene_survival_data)))
cat(sprintf("BCR事件数: %d (%.1f%%)\n", 
            sum(gene_survival_data$bcr_event),
            mean(gene_survival_data$bcr_event) * 100))

# 检查基因表达数据的分布
cat("\n基因表达数据统计:\n")
for(gene in core_genes) {
  cat(sprintf("%s: 均值 = %.3f, SD = %.3f, 范围 = [%.3f, %.3f]\n",
              gene,
              mean(gene_survival_data[[gene]], na.rm = TRUE),
              sd(gene_survival_data[[gene]], na.rm = TRUE),
              min(gene_survival_data[[gene]], na.rm = TRUE),
              max(gene_survival_data[[gene]], na.rm = TRUE)))
}

# 单因素Cox回归分析（每个基因单独分析）
univ_gene_results <- data.frame()

for (gene in core_genes) {
  # 单因素Cox回归
  formula <- as.formula(paste("Surv(survival_time, bcr_event) ~", gene))
  cox_model <- coxph(formula, data = gene_survival_data)
  cox_summary <- summary(cox_model)
  
  univ_gene_results <- rbind(univ_gene_results, data.frame(
    Gene = gene,
    coef = cox_summary$coefficients[1, "coef"],
    HR = cox_summary$conf.int[1, "exp(coef)"],
    HR_lower = cox_summary$conf.int[1, "lower .95"],
    HR_upper = cox_summary$conf.int[1, "upper .95"],
    p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
    Analysis = "Univariate",
    stringsAsFactors = FALSE
  ))
}

# 多因素Cox回归分析（所有基因同时分析）
# 创建公式
multiv_formula <- as.formula(paste(
  "Surv(survival_time, bcr_event) ~", 
  paste(core_genes, collapse = " + ")
))

# 拟合多因素模型
multiv_gene_cox <- coxph(multiv_formula, data = gene_survival_data)
multiv_summary <- summary(multiv_gene_cox)

# 提取多因素结果
multiv_gene_results <- data.frame(
  Gene = core_genes,
  coef = multiv_summary$coefficients[, "coef"],
  HR = multiv_summary$conf.int[, "exp(coef)"],
  HR_lower = multiv_summary$conf.int[, "lower .95"],
  HR_upper = multiv_summary$conf.int[, "upper .95"],
  p_value = multiv_summary$coefficients[, "Pr(>|z|)"],
  Analysis = "Multivariate",
  stringsAsFactors = FALSE
)

# 合并结果
gene_combined_results <- rbind(univ_gene_results, multiv_gene_results)

# 准备绘图数据 - 增加HR范围格式化
gene_plot_data <- gene_combined_results %>%
  mutate(
    # 格式化P值
    P_Value_Formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01 ~ sprintf("%.3f", p_value),
      p_value < 0.05 ~ sprintf("%.3f", p_value),
      TRUE ~ sprintf("%.3f", p_value)
    ),
    # 格式化HR值 - 保留2位小数
    HR_Formatted = sprintf("%.2f", HR),
    # 格式化95% CI - 保留2位小数
    CI_Formatted = sprintf("%.2f-%.2f", HR_lower, HR_upper),
    # 组合HR和95% CI
    HR_with_CI = sprintf("%s (%s)", HR_Formatted, CI_Formatted),
    # 显著性标记
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # 添加显著性标记到HR_with_CI
    HR_with_CI_Sig = paste0(HR_with_CI, Significance),
    # 获取基因的系数（来自LASSO模型）
    LASSO_Coefficient = selected_genes_df$coefficient[match(Gene, selected_genes_df$gene)],
    # 按LASSO系数绝对值排序
    Gene = factor(Gene, levels = core_genes[order(abs(selected_genes_df$coefficient), decreasing = TRUE)]),
    # 分析类型因子
    Analysis = factor(Analysis, levels = c("Univariate", "Multivariate"))
  ) %>%
  arrange(Gene, Analysis) %>%
  group_by(Analysis) %>%
  mutate(
    row_num = row_number(),
    rev_row_num = max(row_num) - row_num + 1  # 反向编号，从上到下显示
  ) %>%
  ungroup()

# 输出结果
cat("\n单因素Cox回归分析结果（核心基因）:\n")
print(univ_gene_results)

cat("\n多因素Cox回归分析结果（核心基因）:\n")
print(multiv_gene_results)

# 模型评估
cat("\n多因素模型统计量:\n")
cat(sprintf("C-index: %.3f\n", multiv_summary$concordance[1]))
cat(sprintf("似然比检验 P值: %.3f\n", multiv_summary$logtest[3]))

# 比例风险假设检验
ph_test <- cox.zph(multiv_gene_cox)
cat("\n比例风险假设检验:\n")
print(ph_test)

# 创建核心基因森林图函数 - 改进版本（HR列显示范围和HR值）
create_gene_forest_plot_with_range <- function(data, title, color, filename) {
  
  # 计算合适的HR显示范围（用于右侧森林图的坐标轴）
  hr_values <- c(data$HR_lower, data$HR, data$HR_upper)
  hr_min <- max(0.01, min(hr_values, na.rm = TRUE))
  hr_max <- min(100, max(hr_values, na.rm = TRUE))
  
  # 根据HR范围确定坐标轴刻度
  if (hr_max > 10) {
    hr_breaks <- c(0.01, 0.1, 1, 10, 100)
    hr_labels <- c("0.01", "0.1", "1", "10", "100")
    x_limits <- c(0.01, 100)
  } else if (hr_max > 5) {
    hr_breaks <- c(0.1, 0.5, 1, 2, 5, 10)
    hr_labels <- c("0.1", "0.5", "1", "2", "5", "10")
    x_limits <- c(0.1, 10)
  } else if (hr_max > 2) {
    hr_breaks <- c(0.2, 0.5, 1, 2, 5)
    hr_labels <- c("0.2", "0.5", "1", "2", "5")
    x_limits <- c(0.2, 5)
  } else {
    hr_breaks <- c(0.3, 0.5, 0.8, 1, 1.5, 2)
    hr_labels <- c("0.3", "0.5", "0.8", "1", "1.5", "2")
    x_limits <- c(0.3, 2)
  }
  
  # 调整坐标轴范围以确保所有点都可见
  buffer_factor <- 0.9
  if (hr_min < x_limits[1]) {
    x_limits[1] <- hr_min * buffer_factor
  }
  if (hr_max > x_limits[2]) {
    x_limits[2] <- hr_max / buffer_factor
  }
  
  cat(sprintf("森林图HR范围: %.3f - %.3f\n", x_limits[1], x_limits[2]))
  
  # 创建表格部分（包含列标题） - 调整列宽适应HR范围显示
  table_plot <- ggplot(data, aes(x = 0, y = rev_row_num)) +
    # 添加列标题背景
    annotate("rect", xmin = 0, xmax = 0.85, ymin = max(data$rev_row_num) + 0.7, 
             ymax = max(data$rev_row_num) + 1.3, fill = "grey90", alpha = 0.8) +
    # 列标题文字 - 调整列宽度
    annotate("text", x = 0.12, y = max(data$rev_row_num) + 1, 
             label = "Gene", hjust = 0, size = 4.0, fontface = "bold") +
    annotate("text", x = 0.35, y = max(data$rev_row_num) + 1, 
             label = "P Value", hjust = 0.5, size = 4.0, fontface = "bold") +
    annotate("text", x = 0.60, y = max(data$rev_row_num) + 1, 
             label = "HR (95% CI)", hjust = 0.5, size = 4.0, fontface = "bold") +
    annotate("text", x = 0.78, y = max(data$rev_row_num) + 1, 
             label = "LASSO Coef", hjust = 0.5, size = 4.0, fontface = "bold") +
    # 数据内容 - 使用HR_with_CI_Sig包含范围和显著性
    geom_text(aes(label = Gene), x = 0.12, hjust = 0, size = 3.8, fontface = "bold") +
    geom_text(aes(label = P_Value_Formatted), x = 0.35, hjust = 0.5, size = 3.8) +
    geom_text(aes(label = HR_with_CI_Sig), x = 0.60, hjust = 0.5, size = 3.5) + # 减小字体以适应范围
    geom_text(aes(label = sprintf("%.3f", LASSO_Coefficient)), 
              x = 0.78, hjust = 0.5, size = 3.8, 
              color = ifelse(data$LASSO_Coefficient > 0, "red", "blue")) +
    # 坐标轴设置
    scale_y_continuous(
      limits = c(0.5, max(data$rev_row_num) + 1.5)
    ) +
    scale_x_continuous(
      limits = c(0, 0.85),
      expand = c(0, 0)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 5, 0, 10),
      plot.title = element_text(size = 14, face = "bold", color = color, hjust = 0.5, 
                                margin = margin(b = 10))
    ) +
    labs(title = title)
  
  # 创建森林图部分
  forest_plot <- ggplot(data, aes(x = HR, y = rev_row_num)) +
    # 参考线（HR=1）
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    # 误差线
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), 
                   height = 0.2, size = 0.9, color = color) +
    # 点
    geom_point(size = 3.5, color = color, shape = 16) +
    # 对数坐标轴（自适应范围）
    scale_x_log10(
      limits = x_limits,
      breaks = hr_breaks,
      labels = hr_labels
    ) +
    # 坐标轴标题
    labs(x = "Hazard Ratio (95% CI)", y = "") +
    # 主题设置
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black"),
      panel.grid.major.x = element_line(color = "grey90", size = 0.4),
      plot.margin = margin(0, 15, 0, 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    # 添加刻度标记
    scale_y_continuous(
      limits = c(0.5, max(data$rev_row_num) + 0.5)
    )
  
  # 合并表格和森林图 - 调整宽度比例
  combined_plot <- table_plot + forest_plot + 
    plot_layout(widths = c(1.4, 1))  # 增加表格部分宽度
  
  # 保存为单独的PDF
  ggsave(file.path(output_dir, filename), 
         combined_plot, width = 15, height = 6, dpi = 300)  # 增加宽度以适应更长的HR范围文本
  
  cat(sprintf("已保存: %s\n", filename))
  
  return(combined_plot)
}

# 创建并保存单因素森林图
gene_univ_plot <- create_gene_forest_plot_with_range(
  gene_plot_data %>% filter(Analysis == "Univariate"),
  "Univariate Cox Analysis of Core EMT Genes",
  "#1B9E77",  # 绿色
  "Core_Genes_Univariate_ForestPlot.pdf"
)

# 创建并保存多因素森林图  
gene_multiv_plot <- create_gene_forest_plot_with_range(
  gene_plot_data %>% filter(Analysis == "Multivariate"),
  "Multivariate Cox Analysis of Core EMT Genes",
  "#D95F02",  # 橙色
  "Core_Genes_Multivariate_ForestPlot.pdf"
)

# 10. 时间依赖性ROC分析
time_points <- c(12, 36, 60)  # 1年、3年、5年

roc_analysis <- timeROC(
  T = validation_survival_data$survival_time,
  delta = validation_survival_data$bcr_event,
  marker = validation_survival_data$emt_score,
  cause = 1,
  times = time_points,
  ROC = TRUE
)

# 准备ROC曲线数据
roc_data <- data.frame(
  Time = rep(c("1-Year", "3-Year", "5-Year"), each = length(roc_analysis$TP[,1])),
  FPR = c(roc_analysis$FP[,1], roc_analysis$FP[,2], roc_analysis$FP[,3]),
  TPR = c(roc_analysis$TP[,1], roc_analysis$TP[,2], roc_analysis$TP[,3])
)

# 创建AUC标注数据 - 在对角线下方的右侧区域
annotation_x <- 0.75  # 统一在x轴0.75位置标注
annotation_y_start <- 0.25  # 第一个标注的y位置
annotation_spacing <- 0.08  # 标注之间的垂直间距

auc_annotations <- data.frame(
  Time = c("1-Year", "3-Year", "5-Year"),
  AUC = roc_analysis$AUC,
  x = annotation_x,
  y = annotation_y_start - (0:2) * annotation_spacing,  # 从上到下排列
  label = sprintf("%s: AUC = %.3f", 
                  c("1-Year", "3-Year", "5-Year"), 
                  roc_analysis$AUC)
)

# 绘制ROC曲线
p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 0.8) +
  # 添加AUC标注
  geom_label(
    data = auc_annotations,
    aes(x = x, y = y, label = label, fill = Time),
    color = "white",
    size = 3.5,
    fontface = "bold",
    hjust = 0.5,
    vjust = 0.5,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  labs(
    title = "Time-Dependent ROC Analysis",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Time Point"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# 保存图形
ggsave(file.path(output_dir, "Time_Dependent_ROC_Validation.pdf"), 
       p_roc, width = 8, height = 6)

# 输出AUC值到控制台
cat("时间依赖性ROC AUC值:\n")
for(i in 1:length(time_points)) {
  cat(sprintf("%d-Year AUC: %.3f\n", 
              time_points[i]/12,
              roc_analysis$AUC[i]))
}

# 生成生存分析汇总表格
survival_summary <- validation_survival_data %>%
  group_by(emt_group) %>%
  summarise(
    n = n(),
    events = sum(bcr_event),
    event_rate = round(events/n * 100, 1),
    median_survival = round(median(survival_time), 1),
    mean_emt = round(mean(emt_score), 3),
    .groups = 'drop'
  )
print(survival_summary)

# 三分位分组汇总
tertile_summary <- validation_survival_data %>%
  group_by(emt_tertile) %>%
  summarise(
    n = n(),
    events = sum(bcr_event),
    event_rate = round(events/n * 100, 1),
    median_survival = round(median(survival_time), 1),
    mean_emt = round(mean(emt_score), 3),
    .groups = 'drop'
  )
print(tertile_summary)

# 11. 多因素和单因素Cox回归分析与组合森林图
multiv_cox_data <- validation_survival_data %>%
  # 选择需要的临床变量
  select(sample_id_tcga, survival_time, bcr_event, emt_score,
         t_stage_group, gleason_total) %>%
  # 添加其他临床变量
  left_join(
    geo_final_data %>% 
      select(GSM, n_stage, psa_level) %>%
      rename(sample_id_tcga = GSM),
    by = "sample_id_tcga"
  ) %>%
  # 数据清理和变量转换
  filter(!is.na(gleason_total) & !is.na(n_stage) & !is.na(psa_level)) %>%
  mutate(
    # T分期：T1/T2 = 0, T3/T4 = 1
    T_Stage_numeric = ifelse(t_stage_group == "T1/T2 (Non-invasive)", 0, 1),
    # Gleason评分：≤7 = 0, ≥8 = 1
    Gleason_numeric = ifelse(gleason_total <= 7, 0, 1),
    # N分期：N0 = 0, N+ = 1
    N_Stage_numeric = case_when(
      n_stage == "N0" ~ 0,
      n_stage %in% c("N1", "N2") ~ 1,
      TRUE ~ NA_real_
    ),
    # PSA水平：连续变量，转换为数值型
    PSA_Numeric = as.numeric(psa_level),
    # EMT评分：连续变量
    EMT_score = emt_score
  ) %>%
  filter(!is.na(T_Stage_numeric) & !is.na(Gleason_numeric) & 
           !is.na(N_Stage_numeric) & !is.na(PSA_Numeric) &
           !is.na(EMT_score))

cat("多因素Cox分析有效样本数:", nrow(multiv_cox_data), "\n")

# 多因素Cox回归分析
multiv_cox <- coxph(Surv(survival_time, bcr_event) ~ 
                      T_Stage_numeric + Gleason_numeric + N_Stage_numeric + 
                      PSA_Numeric + EMT_score, 
                    data = multiv_cox_data)

multiv_summary <- summary(multiv_cox)

# 提取多因素结果 - 简化变量名
multiv_results <- data.frame(
  Variable = c("T_Stage", "Gleason", "N_Stage", "PSA", "EMT"),
  Variable_name = c("T Stage",
                    "Gleason Score", 
                    "N Stage",
                    "PSA Level",
                    "EMT Score"),
  Variable_abbr = c("T Stage (T3-T4 vs T1-T2)",
                    "Gleason Score (≥8 vs ≤7)", 
                    "N Stage (N+ vs N0)",
                    "PSA Level",
                    "EMT Score"),
  coef = multiv_summary$coefficients[, "coef"],
  HR = multiv_summary$conf.int[, "exp(coef)"],
  HR_lower = multiv_summary$conf.int[, "lower .95"],
  HR_upper = multiv_summary$conf.int[, "upper .95"],
  p_value = multiv_summary$coefficients[, "Pr(>|z|)"],
  Analysis = "Multivariate"
)

# 单因素Cox回归分析
univ_variables <- c("T_Stage_numeric", "Gleason_numeric", "N_Stage_numeric", 
                    "PSA_Numeric", "EMT_score")
univ_variable_names <- c("T Stage",
                         "Gleason Score", 
                         "N Stage",
                         "PSA Level",
                         "EMT Score")
univ_variable_abbr <- c("T Stage (T3-T4 vs T1-T2)",
                        "Gleason Score (≥8 vs ≤7)", 
                        "N Stage (N+ vs N0)",
                        "PSA Level",
                        "EMT Score")

univ_results <- data.frame()

for (i in 1:length(univ_variables)) {
  var <- univ_variables[i]
  var_name <- univ_variable_names[i]
  var_abbr <- univ_variable_abbr[i]
  
  formula <- as.formula(paste("Surv(survival_time, bcr_event) ~", var))
  cox_model <- coxph(formula, data = multiv_cox_data)
  
  cox_summary <- summary(cox_model)
  
  univ_results <- rbind(univ_results, data.frame(
    Variable = var,
    Variable_name = var_name,
    Variable_abbr = var_abbr,
    coef = cox_summary$coefficients[1, "coef"],
    HR = cox_summary$conf.int[1, "exp(coef)"],
    HR_lower = cox_summary$conf.int[1, "lower .95"],
    HR_upper = cox_summary$conf.int[1, "upper .95"],
    p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
    Analysis = "Univariate",
    stringsAsFactors = FALSE
  ))
}

# 合并单因素和多因素结果
combined_results <- rbind(univ_results, multiv_results)

# 准备绘图数据 - 增加HR范围格式化
plot_data <- combined_results %>%
  mutate(
    # 格式化P值
    P_Value_Formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    # 格式化HR值 - 保留2位小数
    HR_Formatted = sprintf("%.2f", HR),
    # 格式化95% CI - 保留2位小数
    CI_Formatted = sprintf("%.2f-%.2f", HR_lower, HR_upper),
    # 组合HR和95% CI
    HR_with_CI = sprintf("%s (%s)", HR_Formatted, CI_Formatted),
    # 显著性标记
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # 添加显著性标记到HR_with_CI
    HR_with_CI_Sig = paste0(HR_with_CI, Significance),
    # 固定变量顺序 - 使用简化变量名
    Variable_name = factor(Variable_name, 
                           levels = c("T Stage",
                                      "Gleason Score", 
                                      "N Stage",
                                      "PSA Level",
                                      "EMT Score")),
    # 分析类型因子
    Analysis = factor(Analysis, levels = c("Univariate", "Multivariate"))
  ) %>%
  arrange(Variable_name, Analysis) %>%
  group_by(Analysis) %>%
  mutate(
    row_num = row_number(),
    rev_row_num = max(row_num) - row_num + 1  # 反向编号，从上到下显示
  ) %>%
  ungroup()

# 输出格式化结果到控制台（包含完整信息）
cat("\n=== 格式化的Cox回归分析结果 ===\n")
cat("\n单因素分析:\n")
for(i in 1:nrow(univ_results)) {
  var_name <- univ_results$Variable_name[i]
  var_abbr <- univ_results$Variable_abbr[i]
  hr <- sprintf("%.2f", univ_results$HR[i])
  ci_lower <- sprintf("%.2f", univ_results$HR_lower[i])
  ci_upper <- sprintf("%.2f", univ_results$HR_upper[i])
  p_value <- ifelse(univ_results$p_value[i] < 0.001, "<0.001", 
                    sprintf("%.3f", univ_results$p_value[i]))
  
  cat(sprintf("%s: HR = %s (95%% CI: %s-%s), P = %s\n", 
              var_abbr, hr, ci_lower, ci_upper, p_value))
}

cat("\n多因素分析:\n")
for(i in 1:nrow(multiv_results)) {
  var_name <- multiv_results$Variable_name[i]
  var_abbr <- multiv_results$Variable_abbr[i]
  hr <- sprintf("%.2f", multiv_results$HR[i])
  ci_lower <- sprintf("%.2f", multiv_results$HR_lower[i])
  ci_upper <- sprintf("%.2f", multiv_results$HR_upper[i])
  p_value <- ifelse(multiv_results$p_value[i] < 0.001, "<0.001", 
                    sprintf("%.3f", multiv_results$p_value[i]))
  
  cat(sprintf("%s: HR = %s (95%% CI: %s-%s), P = %s\n", 
              var_abbr, hr, ci_lower, ci_upper, p_value))
}

# 创建带有HR范围的森林图函数（使用简化变量名）
create_clinical_forest_plot_with_range <- function(data, title, color, filename) {
  
  # 计算合适的HR显示范围（用于右侧森林图的坐标轴）
  hr_values <- c(data$HR_lower, data$HR, data$HR_upper)
  hr_min <- max(0.01, min(hr_values, na.rm = TRUE))
  hr_max <- min(100, max(hr_values, na.rm = TRUE))
  
  # 根据HR范围确定坐标轴刻度
  if (hr_max > 10) {
    hr_breaks <- c(0.01, 0.1, 1, 10, 100)
    hr_labels <- c("0.01", "0.1", "1", "10", "100")
    x_limits <- c(0.01, 100)
  } else if (hr_max > 5) {
    hr_breaks <- c(0.1, 0.5, 1, 2, 5, 10)
    hr_labels <- c("0.1", "0.5", "1", "2", "5", "10")
    x_limits <- c(0.1, 10)
  } else if (hr_max > 2) {
    hr_breaks <- c(0.2, 0.5, 1, 2, 5)
    hr_labels <- c("0.2", "0.5", "1", "2", "5")
    x_limits <- c(0.2, 5)
  } else {
    hr_breaks <- c(0.3, 0.5, 0.8, 1, 1.5, 2)
    hr_labels <- c("0.3", "0.5", "0.8", "1", "1.5", "2")
    x_limits <- c(0.3, 2)
  }
  
  # 调整坐标轴范围以确保所有点都可见
  buffer_factor <- 0.9
  if (hr_min < x_limits[1]) {
    x_limits[1] <- hr_min * buffer_factor
  }
  if (hr_max > x_limits[2]) {
    x_limits[2] <- hr_max / buffer_factor
  }
  
  cat(sprintf("森林图HR范围: %.3f - %.3f\n", x_limits[1], x_limits[2]))
  
  # 创建表格部分 - 优化列宽和布局
  table_plot <- ggplot(data, aes(x = 0, y = rev_row_num)) +
    # 添加列标题背景
    annotate("rect", xmin = 0, xmax = 0.65, ymin = max(data$rev_row_num) + 0.7, 
             ymax = max(data$rev_row_num) + 1.3, fill = "grey90", alpha = 0.8) +
    # 列标题文字 - 重新分配列宽度（简化变量名后可以恢复原宽度）
    annotate("text", x = 0.15, y = max(data$rev_row_num) + 1, 
             label = "Variable", hjust = 0, size = 4.0, fontface = "bold") +
    annotate("text", x = 0.38, y = max(data$rev_row_num) + 1, 
             label = "P Value", hjust = 0.5, size = 4.0, fontface = "bold") +
    annotate("text", x = 0.55, y = max(data$rev_row_num) + 1, 
             label = "HR (95% CI)", hjust = 0.5, size = 4.0, fontface = "bold") +
    # 数据内容 - 使用简化变量名
    geom_text(aes(label = Variable_name), x = 0.15, hjust = 0, size = 3.8) +
    geom_text(aes(label = P_Value_Formatted), x = 0.38, hjust = 0.5, size = 3.8) +
    geom_text(aes(label = HR_with_CI_Sig), x = 0.55, hjust = 0.5, size = 3.5) +
    # 坐标轴设置
    scale_y_continuous(
      limits = c(0.5, max(data$rev_row_num) + 1.5)
    ) +
    scale_x_continuous(
      limits = c(0, 0.65),
      expand = c(0, 0)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 5, 0, 10),
      plot.title = element_text(size = 12, face = "bold", color = color, hjust = 0)
    ) +
    labs(title = title)
  
  # 创建森林图部分
  forest_plot <- ggplot(data, aes(x = HR, y = rev_row_num)) +
    # 参考线
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    # 误差线
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), 
                   height = 0.2, size = 0.9, color = color) +
    # 点
    geom_point(size = 3.5, color = color, shape = 16) +
    # 对数坐标轴（自适应范围）
    scale_x_log10(
      limits = x_limits,
      breaks = hr_breaks,
      labels = hr_labels
    ) +
    # 坐标轴标题
    labs(x = "Hazard Ratio (95% CI)", y = "") +
    # 主题设置
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black"),
      panel.grid.major.x = element_line(color = "grey90", size = 0.4),
      plot.margin = margin(0, 15, 0, 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    # 添加刻度标记
    scale_y_continuous(
      limits = c(0.5, max(data$rev_row_num) + 0.5)
    )
  
  # 合并表格和森林图
  combined <- table_plot + forest_plot + 
    plot_layout(widths = c(1, 1.5))
  
  # 保存为单独的PDF
  ggsave(file.path(output_dir, filename), 
         combined, width = 14, height = 6, dpi = 300)
  
  cat(sprintf("森林图已保存: %s\n", filename))
  
  return(combined)
}

# 创建并保存单因素森林图
univ_plot <- create_clinical_forest_plot_with_range(
  plot_data %>% filter(Analysis == "Univariate"),
  "Univariate Cox Regression Analysis",
  "#377EB8",  # 蓝色
  "Univariate_Cox_ForestPlot.pdf"
)

# 创建并保存多因素森林图
multiv_plot <- create_clinical_forest_plot_with_range(
  plot_data %>% filter(Analysis == "Multivariate"),
  "Multivariate Cox Regression Analysis",
  "#D55E00",  # 橙色
  "Multivariate_Cox_ForestPlot.pdf"
)
# 12. 列线图和校准曲线分析
library(rms)
library(pec)

nomogram_data <- validation_survival_data %>%
  select(sample_id_tcga, survival_time, bcr_event, emt_score,
         t_stage_group, gleason_total) %>%
  left_join(
    geo_final_data %>% 
      select(GSM, n_stage) %>%
      rename(sample_id_tcga = GSM),
    by = "sample_id_tcga"
  ) %>%
  filter(!is.na(gleason_total) & !is.na(n_stage)) %>%
  mutate(
    T_Stage_numeric = ifelse(t_stage_group == "T1/T2 (Non-invasive)", 0, 1),
    Gleason_numeric = ifelse(gleason_total <= 7, 0, 1),
    N_Stage_numeric = case_when(
      n_stage == "N0" ~ 0,
      n_stage %in% c("N1", "N2") ~ 1,
      TRUE ~ NA_real_
    ),
    EMT_score = emt_score
  ) %>%
  filter(!is.na(T_Stage_numeric) & !is.na(Gleason_numeric) & 
           !is.na(N_Stage_numeric) & !is.na(EMT_score))

cat(sprintf("样本数: %d\n", nrow(nomogram_data)))
cat(sprintf("事件数: %d (%.1f%%)\n", 
            sum(nomogram_data$bcr_event),
            mean(nomogram_data$bcr_event) * 100))

# 设置rms包的工作环境
dd <- datadist(nomogram_data)
options(datadist = 'dd')

# 创建列线图模型
cph_for_nomogram <- cph(Surv(survival_time, bcr_event) ~ 
                          T_Stage_numeric + Gleason_numeric + N_Stage_numeric + 
                          EMT_score,
                        data = nomogram_data,
                        surv = TRUE, x = TRUE, y = TRUE)

cat(sprintf("C-index: %.3f\n", cph_for_nomogram$stats["C"]))

# 创建列线图
tryCatch({
  surv_prob_func <- function(x, time) {
    base_surv <- survest(cph_for_nomogram, 
                         newdata = data.frame(T_Stage_numeric = 0,
                                              Gleason_numeric = 0,
                                              N_Stage_numeric = 0,
                                              EMT_score = 0),
                         times = time)$surv
    return(1 - base_surv^exp(x))
  }
  
  nomogram_plot <- nomogram(cph_for_nomogram,
                            fun = list('1-Year' = function(x) surv_prob_func(x, 12),
                                       '3-Year' = function(x) surv_prob_func(x, 36),
                                       '5-Year' = function(x) surv_prob_func(x, 60)),
                            fun.at = list(c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                                          c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                                          c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)),
                            funlabel = c("1-Year BCR Risk",
                                         "3-Year BCR Risk", 
                                         "5-Year BCR Risk"),
                            lp = FALSE)
  
  pdf(file.path(output_dir, "EMT_Nomogram_NoPSA.pdf"), width = 14, height = 10)
  plot(nomogram_plot, lmgp = 0.2, cex.axis = 0.7, cex.var = 0.9)
  dev.off()
  cat("列线图已保存: EMT_Nomogram_NoPSA.pdf\n")
  
}, error = function(e) {
  cat("列线图创建失败，尝试简化版本...\n")
  
  nomogram_plot_simple <- nomogram(cph_for_nomogram,
                                   fun = function(x) 1 - (0.85)^exp(x),
                                   fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                                   funlabel = "3-Year BCR Risk",
                                   lp = FALSE)
  
  pdf(file.path(output_dir, "EMT_Nomogram_NoPSA_Simple.pdf"), width = 12, height = 8)
  plot(nomogram_plot_simple, lmgp = 0.3, cex.axis = 0.8, cex.var = 0.9)
  dev.off()
  cat("简化版列线图已保存: EMT_Nomogram_NoPSA_Simple.pdf\n")
})

# 校准曲线
set.seed(1234)

# 创建校准曲线
cal_1y <- calibrate(
  cph_for_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 12,          # 1年
  m = 50,
  B = 1000
)

cal_3y <- calibrate(
  cph_for_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 36,          # 3年
  m = 50,
  B = 1000
)

cal_5y <- calibrate(
  cph_for_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 60,          # 5年
  m = 50,
  B = 1000
)

# 准备校准数据
prepare_cal_data <- function(cal_object, time_label) {
  data.frame(
    pred_prob = cal_object[, "mean.predicted"],      # 预测风险
    obs_prob = cal_object[, "KM"],                   # 观察风险
    obs_corrected = cal_object[, "KM.corrected"],    # 校正后观察风险
    std_err = cal_object[, "std.err"],               # 标准误
    time = time_label
  )
}

# 合并数据
cal_data <- rbind(
  prepare_cal_data(cal_1y, "1-Year"),
  prepare_cal_data(cal_3y, "3-Year"),
  prepare_cal_data(cal_5y, "5-Year")
)

# 设置专业配色
colors <- c("1-Year" = "#E41A1C",    # 红色
            "3-Year" = "#377EB8",    # 蓝色
            "5-Year" = "#4DAF4A")    # 绿色

# 创建图形
p_cal <- ggplot(cal_data, aes(x = pred_prob, y = obs_corrected, 
                              color = time, linetype = time, shape = time)) +
  
  # 对角线参考线（理想校准线）
  geom_abline(intercept = 0, slope = 1, 
              color = "gray60", 
              linetype = "dashed", 
              size = 0.8) +
  
  # 校准曲线
  geom_line(size = 1.2) +
  
  # 数据点 - 改为实心
  geom_point(size = 3.5) +  # 移除了fill和stroke参数，使用实心点
  
  # 误差线
  geom_errorbar(aes(ymin = obs_corrected - 1.96 * std_err,
                    ymax = obs_corrected + 1.96 * std_err),
                width = 0.02, 
                size = 0.7) +
  
  # 坐标轴标签 - 移除了标题
  labs(
    x = "Predicted BCR Risk",
    y = "Observed BCR Risk (Bias-Corrected)",
    color = "Follow-up Time",
    linetype = "Follow-up Time",
    shape = "Follow-up Time"
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_blank(),  # 移除了标题
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", size = 0.7),
    axis.ticks = element_line(color = "black", size = 0.7),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),
    panel.grid.major = element_line(color = "gray92", size = 0.2),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # 颜色和样式映射
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("1-Year" = "solid", 
                                   "3-Year" = "dashed", 
                                   "5-Year" = "dotted")) +
  scale_shape_manual(values = c("1-Year" = 16,    # 实心圆形
                                "3-Year" = 17,    # 实心三角形
                                "5-Year" = 15)) + # 实心方形
  
  # 坐标轴范围（0-100%）
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  
  # 百分比格式
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.02)))

ggsave(file.path(output_dir, "Calibration_Curves_Final.pdf"), 
       p_cal, width = 8, height = 6, dpi = 300)

# 13. 三图组合可视化
# 准备数据
validation_risk_data <- validation_survival_data %>%
  arrange(emt_score) %>%
  mutate(
    patient_id = 1:n(),  # 重新编号为1,2,3...
    risk_group = ifelse(emt_score > median(emt_score), "High", "Low"),
    bcr_status = ifelse(bcr_event == 1, "Recurrence", "Non-recurrence"),
    risk_group = factor(risk_group, levels = c("Low", "High")),
    bcr_status = factor(bcr_status, levels = c("Non-recurrence", "Recurrence")),
    # 为散点图添加随机抖动
    jitter_y = runif(n(), 0.8, 1.2)  # Y坐标随机抖动
  )

# 计算高低风险组分界点
low_risk_count <- sum(validation_risk_data$risk_group == "Low")
cutoff_point <- low_risk_count + 0.5

# 颜色方案（保持原脚本的蓝-红配色）
risk_colors <- c("Low" = "#377EB8", "High" = "#E41A1C")
bcr_colors <- c("Non-recurrence" = "#377EB8", "Recurrence" = "#E41A1C")

# 图1：EMT风险评分分布图
p1 <- ggplot(validation_risk_data, aes(x = patient_id)) +
  # 使用geom_ribbon创建填充区域
  geom_ribbon(aes(ymin = median(emt_score), 
                  ymax = emt_score, 
                  fill = risk_group), 
              alpha = 0.7) +
  
  # 中位线
  geom_hline(yintercept = median(validation_risk_data$emt_score), 
             colour = "black", 
             linetype = "dashed", 
             size = 0.8, 
             alpha = 0.7) +
  
  # 高低风险组分隔线
  geom_vline(xintercept = cutoff_point,
             colour = "black", 
             linetype = "solid", 
             size = 1) +
  
  # 添加分组标签
  annotate("text", 
           x = low_risk_count/2, 
           y = max(validation_risk_data$emt_score) * 1.05,
           label = "Low Risk", 
           color = "#377EB8", 
           fontface = "bold", 
           size = 4.5) +
  
  annotate("text", 
           x = low_risk_count + (nrow(validation_risk_data) - low_risk_count)/2,
           y = max(validation_risk_data$emt_score) * 1.05,
           label = "High Risk", 
           color = "#E41A1C", 
           fontface = "bold", 
           size = 4.5) +
  
  # 颜色填充
  scale_fill_manual(values = risk_colors) +
  
  # 坐标轴标签
  labs(x = "Patients (ordered by EMT risk score)",
       y = "EMT Risk Score") +
  
  # 主题设置 - 添加X轴显示
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey92", size = 0.3),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  
  # 调整X轴刻度
  scale_x_continuous(
    breaks = c(1, round(nrow(validation_risk_data)/4), 
               round(nrow(validation_risk_data)/2),
               round(nrow(validation_risk_data)*3/4), 
               nrow(validation_risk_data)),
    labels = c("1", "25%", "50%", "75%", "100%")
  )

# 保存图1
ggsave(file.path(output_dir, "Panel1_EMT_Risk_Score_Distribution.pdf"),
       p1, width = 12, height = 5, dpi = 300)

# 图2：BCR状态分布图
p2 <- ggplot(validation_risk_data, aes(x = patient_id, y = jitter_y)) +
  # 数据点（使用散点，加入随机抖动）
  geom_point(aes(color = bcr_status, shape = bcr_status), 
             size = 2.5, 
             alpha = 0.7) +
  
  # 高低风险组分隔线
  geom_vline(xintercept = cutoff_point,
             colour = "black", 
             linetype = "solid", 
             size = 1) +
  
  # 添加Y轴标签
  labs(x = "Patients (ordered by EMT risk score)",
       y = "BCR Status",
       color = "BCR Status",
       shape = "BCR Status") +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none",  # 移除图例
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  
  # 调整X轴刻度（与图1保持一致）
  scale_x_continuous(
    breaks = c(1, round(nrow(validation_risk_data)/4), 
               round(nrow(validation_risk_data)/2),
               round(nrow(validation_risk_data)*3/4), 
               nrow(validation_risk_data)),
    labels = c("1", "25%", "50%", "75%", "100%")
  ) +
  
  # 调整Y轴范围和标签
  scale_y_continuous(
    limits = c(0.7, 1.3),
    breaks = c(0.8, 1.0, 1.2),
    labels = c("", "", "")  # 隐藏Y轴刻度标签
  ) +
  
  # 颜色和形状设置（但图例已移除）
  scale_color_manual(values = bcr_colors) +
  scale_shape_manual(values = c("Non-recurrence" = 16,  # 实心圆
                                "Recurrence" = 17))     # 三角形

# 保存图2
ggsave(file.path(output_dir, "Panel2_BCR_Status_Distribution.pdf"),
       p2, width = 12, height = 4, dpi = 300)

# 图3：基因表达热图
p3_complex <- Heatmap(
  heatmap_data_scaled,
  name = "Z-score",  # 右侧只显示Z-score
  
  # 颜色映射
  col = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")),
  
  # 完全移除顶部注释
  top_annotation = NULL,  # 移除所有顶部注释
  
  # 基因行设置
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  show_row_names = TRUE,
  row_title = "EMT signature genes",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 样本列设置
  show_column_names = FALSE,
  column_split = column_anno$Risk_Group,  # 仍按风险组分隔（保持垂直线分隔）
  column_gap = unit(2, "mm"),
  column_title = NULL,  # 移除列标题
  
  # 聚类设置
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # 图例设置（只保留Z-score图例）
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10, face = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(3, "cm"),
    direction = "vertical"
  ),
  
  # 隐藏其他图例
  show_heatmap_legend = TRUE
)

# 保存图3
pdf(file.path(output_dir, "Panel3_EMT_Gene_Expression_Heatmap.pdf"),
    width = 12, height = 6)
draw(p3_complex, 
     heatmap_legend_side = "right",  # Z-score图例在右侧
     merge_legend = TRUE)
dev.off()

# 图片拼接
pdf(NULL)  # 创建空PDF设备
heatmap_grob <- grid.grabExpr(
  draw(p3_complex, 
       heatmap_legend_side = "right",
       merge_legend = TRUE),
  wrap = TRUE
)
dev.off()

# 现在使用grid.arrange正确拼接
pdf(file.path(output_dir, "Three_Panels_Vertical_Combined.pdf"), 
    width = 12, height = 16)

grid.arrange(
  p1,            # 图1：EMT风险评分分布图
  p2,            # 图2：BCR状态分布图  
  heatmap_grob,  # 图3：基因表达热图（已转换为grob）
  ncol = 1,
  heights = c(1.2, 1, 2)  # 调整高度比例
)

dev.off()



# TCGA生存分析
#tcga_logtpm <- readRDS("/home/datahup/xzp/paper/PCa/data/AUcell-Model/TCGA_PRAD_tpm.rds")
#tcga_clinical <- read.delim("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/clinical_data.tsv")
#cat("TCGA表达数据维度:", dim(tcga_logtpm), "\n")
#cat("TCGA临床数据维度:", dim(tcga_clinical), "\n")

# 提取TCGA BCR生存数据
#tcga_bcr_data <- tcga_clinical %>%
#  mutate(
#    sample_id = Sample.ID,
#    sample_id_clean = gsub("-", ".", gsub("-01$", "", sample_id)),
#    bcr_time = as.numeric(Disease.Free..Months.),
#    bcr_status = case_when(
#      grepl("^0|DiseaseFree", Disease.Free.Status, ignore.case = TRUE) ~ 0,
#      grepl("^1|Recurred", Disease.Free.Status, ignore.case = TRUE) ~ 1,
#      TRUE ~ NA_real_
#    )
#  ) %>%
# filter(!is.na(bcr_status) & !is.na(bcr_time) & bcr_time > 0) %>%
#  dplyr::select(sample_id, sample_id_clean, bcr_time, bcr_status)
#cat("TCGA有效BCR数据样本数:", nrow(tcga_bcr_data), "\n")
#print(table(tcga_bcr_data$bcr_status))
#cat("事件率:", round(mean(tcga_bcr_data$bcr_status) * 100, 1), "%\n")

# TCGA样本匹配和EMT评分计算
# 创建TCGA表达数据映射表
#tcga_expression_mapping <- data.frame(
#  original_id = colnames(tcga_logtpm),
#  dot_id = gsub("-", ".", colnames(tcga_logtpm))
#)

# 检查重复
#tcga_duplicates <- tcga_expression_mapping$dot_id[duplicated(tcga_expression_mapping$dot_id)]
#cat("TCGA重复的dot_id数量:", length(tcga_duplicates), "\n")

# 找到真正匹配的样本
#actually_matched_tcga <- tcga_bcr_data$sample_id_clean[tcga_bcr_data$sample_id_clean %in% tcga_expression_mapping$dot_id]
#cat("实际匹配的TCGA样本数:", length(actually_matched_tcga), "\n")

# 过滤BCR数据
#tcga_bcr_filtered <- tcga_bcr_data %>%
#  filter(sample_id_clean %in% actually_matched_tcga) %>%
#  arrange(sample_id_clean)

# 获取对应的表达数据列名（去重）
#tcga_matched_mapping <- tcga_expression_mapping %>%
#  filter(dot_id %in% tcga_bcr_filtered$sample_id_clean) %>%
#  distinct(dot_id, .keep_all = TRUE) %>%
#  arrange(match(dot_id, tcga_bcr_filtered$sample_id_clean))
#tcga_matched_original_ids <- tcga_matched_mapping$original_id

# 提取表达数据
#tcga_expression_final <- tcga_logtpm[, tcga_matched_original_ids, drop = FALSE]

# 检查模型基因在TCGA中的可用性
#available_genes_tcga <- selected_genes[selected_genes %in% rownames(tcga_expression_final)]
#cat("TCGA中可用的EMT基因:", length(available_genes_tcga), "/", length(selected_genes), "\n")

# 计算TCGA EMT评分
#X_tcga_final <- t(tcga_expression_final[available_genes_tcga, , drop = FALSE])
#emt_scores_tcga_final <- as.numeric(X_tcga_final %*% selected_coef[available_genes_tcga])

#cat("最终TCGA BCR数据样本数:", nrow(tcga_bcr_filtered), "\n")
#cat("最终TCGA表达数据样本数:", ncol(tcga_expression_final), "\n")
#cat("样本数量一致:", nrow(tcga_bcr_filtered) == length(emt_scores_tcga_final), "\n")

# 合并TCGA生存数据（两分组）
#cat("\n=== 合并TCGA生存数据 ===\n")
#tcga_survival_data <- tcga_bcr_filtered %>%
#  mutate(
#    emt_score = emt_scores_tcga_final,
#    emt_group = ifelse(emt_score > median(emt_score), "High_EMT", "Low_EMT")
#  )

#cat("TCGA生存分析最终样本数:", nrow(tcga_survival_data), "\n")
#cat("中位数:", round(median(tcga_survival_data$emt_score), 3), "\n")
#cat("范围:", round(min(tcga_survival_data$emt_score), 3), "-", 
#    round(max(tcga_survival_data$emt_score), 3), "\n")
#print(table(tcga_survival_data$emt_group))

# TCGA Kaplan-Meier生存分析（两分组）
#km_fit_tcga <- survfit(Surv(bcr_time, bcr_status) ~ emt_group, data = tcga_survival_data)

#p_km_tcga <- ggsurvplot(
#  km_fit_tcga,
#  data = tcga_survival_data,
#  pval = TRUE,
#  pval.method = TRUE,
#  conf.int = TRUE,
#  risk.table = TRUE,
#  surv.median.line = "hv",
#  palette = c("#E41A1C", "#377EB8"),
#  legend.labs = c("High EMT", "Low EMT"),
#  xlab = "Time to BCR (Months)",
#  ylab = "BCR-Free Survival Probability",
#  ggtheme = theme_classic()
#)

# 保存TCGA KM曲线
#ggsave(file.path(output_dir, "TCGA_BCR_KM_Curve_EMT.pdf"), 
#       p_km_tcga$plot, width = 8, height = 8)

# TCGA Cox回归分析
# 单变量Cox模型（连续EMT评分）
#cox_tcga_continuous <- coxph(Surv(bcr_time, bcr_status) ~ emt_score, data = tcga_survival_data)
#cox_tcga_summary <- summary(cox_tcga_continuous)

#cat("TCGA单变量Cox回归分析 (连续EMT评分):\n")
#cat("Hazard Ratio:", round(exp(cox_tcga_continuous$coefficients), 3), "\n")
#cat("95% CI:", round(exp(confint(cox_tcga_continuous)), 3), "\n")
#cat("p-value:", signif(cox_tcga_summary$coefficients[,"Pr(>|z|)"], 3), "\n")

# 单变量Cox模型（分类EMT评分）
#cox_tcga_categorical <- coxph(Surv(bcr_time, bcr_status) ~ emt_group, data = tcga_survival_data)
#cox_tcga_cat_summary <- summary(cox_tcga_categorical)

#cat("\nTCGA单变量Cox回归分析 (分类EMT评分):\n")
#cat("Hazard Ratio:", round(exp(cox_tcga_categorical$coefficients), 3), "\n")
#cat("95% CI:", round(exp(confint(cox_tcga_categorical)), 3), "\n")
#cat("p-value:", signif(cox_tcga_cat_summary$coefficients[,"Pr(>|z|)"], 3), "\n")

# TCGA时间依赖性ROC分析
#roc_tcga <- timeROC(
#  T = tcga_survival_data$bcr_time,
#  delta = tcga_survival_data$bcr_status,
#  marker = tcga_survival_data$emt_score,
#  cause = 1,
#  times = time_points,
#  ROC = TRUE
#)
#  cat(sprintf("%d-Year AUC: %.3f\n", time_points[i]/12, roc_tcga$AUC[i]))

# 绘制TCGA ROC曲线
#roc_data_tcga <- data.frame(
#  Time = rep(c("1-Year", "3-Year", "5-Year"), each = length(roc_tcga$TP[,1])),
#  FPR = c(roc_tcga$FP[,1], roc_tcga$FP[,2], roc_tcga$FP[,3]),
#  TPR = c(roc_tcga$TP[,1], roc_tcga$TP[,2], roc_tcga$TP[,3])
#)

#p_roc_tcga <- ggplot(roc_data_tcga, aes(x = FPR, y = TPR, color = Time)) +
#  geom_line(size = 1.2) +
#  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
#  labs(
#    x = "False Positive Rate",
#    y = "True Positive Rate",
#    color = "Time Point"
#  ) +
#  theme_classic() +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold"),
#    legend.position = "bottom"
#  ) +
#  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))

#ggsave(file.path(output_dir, "TCGA_Time_Dependent_ROC.pdf"), 
#       p_roc_tcga, width = 8, height = 6)

# TCGA生存分析结果汇总
# 生成TCGA生存分析汇总表格
#survival_summary_tcga <- tcga_survival_data %>%
#  group_by(emt_group) %>%
#  summarise(
#    n = n(),
#    events = sum(bcr_status),
#    event_rate = round(events/n * 100, 1),
#    median_survival = round(median(bcr_time), 1),
#    mean_emt = round(mean(emt_score), 3),
#    median_emt = round(median(emt_score), 3),
#    .groups = 'drop'
#  )
#print(survival_summary_tcga)
#saveRDS(tcga_survival_data, "TCGA_EMT_Survival_Data.rds")