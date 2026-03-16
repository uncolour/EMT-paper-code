# Step 4-2: Mutation Analysis
rm(list = ls())
gc()

library(dplyr)
library(maftools)
library(stringr)
library(ggplot2)
library(ggpubr)

# 设置工作目录和输出路径
setwd('/home/datahup/xzp/paper/PCa/data/Mutation Analysis')
output_dir <- '/home/datahup/xzp/paper/PCa/fig/Mutation Analysis'

# 1. 加载突变数据
cat("=== 1. 加载突变数据 ===\n")
maf_gz_files <- list.files("extracted_maf", pattern = "\\.maf\\.gz$", recursive = TRUE, full.names = TRUE)
cat("找到", length(maf_gz_files), "个maf.gz文件\n")

# 读取MAF文件
maf_list <- list()
for (i in 1:length(maf_gz_files)) {
  if (i %% 50 == 0) cat("正在处理第", i, "/", length(maf_gz_files), "个文件\n")
  
  tryCatch({
    maf_single <- read.maf(maf = maf_gz_files[i])
    maf_list[[i]] <- maf_single
  }, error = function(e) {
    cat("  警告: 文件", basename(maf_gz_files[i]), "读取失败:", e$message, "\n")
  })
}

# 移除失败的文件
maf_list <- maf_list[!sapply(maf_list, is.null)]
cat("成功读取", length(maf_list), "个MAF文件\n")

# 合并MAF文件
maf_combined <- merge_mafs(maf_list)
mutation_samples <- as.character(maf_combined@data$Tumor_Sample_Barcode)
unique_mutation_samples <- unique(mutation_samples)
mutation_samples_short <- substr(unique_mutation_samples, 1, 12)

# 去重处理
mutation_unique_indices <- !duplicated(mutation_samples_short)
unique_mutation_samples_final <- unique_mutation_samples[mutation_unique_indices]
mutation_samples_short_final <- mutation_samples_short[mutation_unique_indices]

cat("合并后突变样本数:", length(unique_mutation_samples), "\n")
cat("去重后突变样本数:", length(unique_mutation_samples_final), "\n")

# 保存合并后的MAF
saveRDS(maf_combined, "TCGA_combined_maf.rds")
cat("合并MAF已保存: TCGA_combined_maf.rds\n")

# 2. 加载表达数据并去重
cat("\n=== 2. 加载表达数据 ===\n")
expression_dir <- '/home/datahup/xzp/paper/PCa/data/AUcell-Model'
logtpm_matrix <- readRDS(file.path(expression_dir, "TCGA_PRAD_tpm.rds"))

expression_samples <- colnames(logtpm_matrix)
expression_samples_short <- substr(expression_samples, 1, 12)

# 表达数据去重
expression_unique_indices <- !duplicated(expression_samples_short)
logtpm_unique <- logtpm_matrix[, expression_unique_indices]
expression_samples_short_unique <- expression_samples_short[expression_unique_indices]

cat("原始表达样本数:", length(expression_samples), "\n")
cat("去重后表达样本数:", length(expression_samples_short_unique), "\n")

# 3. 样本匹配
cat("\n=== 3. 样本匹配 ===\n")
overlapping_samples <- intersect(expression_samples_short_unique, mutation_samples_short_final)
cat("重叠样本数:", length(overlapping_samples), "\n")

# 4. EMT评分计算
cat("\n=== 4. 计算EMT评分 ===\n")
emt_model_path <- "../TCGA GEO-Model/seed_125_lasso_results/core_emt_genes.csv"
selected_genes_df <- read.csv(emt_model_path)
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

cat("EMT模型基因数:", length(selected_genes), "\n")

# 提取重叠样本的表达数据
overlapping_expression_indices <- expression_samples_short_unique %in% overlapping_samples
logtpm_overlapping <- logtpm_unique[, overlapping_expression_indices]
expression_samples_overlapping <- expression_samples_short_unique[overlapping_expression_indices]

cat("用于EMT评分的样本数:", ncol(logtpm_overlapping), "\n")

# 提取EMT基因表达数据
tcga_selected_genes <- logtpm_overlapping[rownames(logtpm_overlapping) %in% selected_genes, ]
cat("在TCGA中找到的EMT基因数:", nrow(tcga_selected_genes), "\n")

# 检查缺失基因
missing_genes <- setdiff(selected_genes, rownames(logtpm_overlapping))
if(length(missing_genes) > 0) {
  cat("警告: 以下基因在TCGA数据中缺失:", paste(missing_genes, collapse = ", "), "\n")
}

# 准备预测数据
tcga_X <- t(tcga_selected_genes)
tcga_X <- as.data.frame(tcga_X)
tcga_X <- tcga_X[, selected_genes, drop = FALSE]

# 计算EMT评分
tcga_emt_scores <- as.numeric(as.matrix(tcga_X) %*% selected_coef)

# 创建EMT分组
emt_median <- median(tcga_emt_scores, na.rm = TRUE)
tcga_emt_data <- data.frame(
  sample_id = colnames(logtpm_overlapping),
  sample_id_short = expression_samples_overlapping,
  emt_score = tcga_emt_scores,
  EMTGroup = ifelse(tcga_emt_scores > emt_median, "High_EMT", "Low_EMT")
)

cat("EMT评分统计:\n")
print(summary(tcga_emt_scores))
cat("高EMT组:", sum(tcga_emt_data$EMTGroup == "High_EMT"), "\n")
cat("低EMT组:", sum(tcga_emt_data$EMTGroup == "Low_EMT"), "\n")

# 创建样本映射表
cat("\n=== 5. 创建样本映射表 ===\n")
sample_mapping <- data.frame()
for(short_id in overlapping_samples) {
  # 在表达数据中找到对应的样本
  expr_match <- tcga_emt_data[tcga_emt_data$sample_id_short == short_id, ]
  # 在突变数据中找到对应的样本
  mut_match <- unique_mutation_samples_final[mutation_samples_short_final == short_id]
  
  sample_mapping <- rbind(sample_mapping, data.frame(
    expression_id = expr_match$sample_id,
    mutation_id = mut_match,
    short_id = short_id,
    emt_score = expr_match$emt_score,
    EMTGroup = expr_match$EMTGroup
  ))
}

cat("最终匹配样本数:", nrow(sample_mapping), "\n")

# 保存EMT评分数据
write.csv(tcga_emt_data, "TCGA_EMT_Scores.csv", row.names = FALSE)
cat("EMT评分数据已保存: TCGA_EMT_Scores.csv\n")

# 创建MAF子集
matched_mutation_samples <- sample_mapping$mutation_id
maf_subset_data <- maf_combined@data[maf_combined@data$Tumor_Sample_Barcode %in% matched_mutation_samples, ]
maf_subset_clin <- maf_combined@clinical.data[maf_combined@clinical.data$Tumor_Sample_Barcode %in% matched_mutation_samples, ]

# 创建新的MAF对象
maf_subset <- read.maf(maf = maf_subset_data,
                       clinicalData = maf_subset_clin,
                       isTCGA = TRUE)

cat("子集MAF样本数:", length(unique(maf_subset@data$Tumor_Sample_Barcode)), "\n")
cat("子集MAF突变数:", nrow(maf_subset@data), "\n")

# 添加EMT分组信息
maf_subset@clinical.data$EMTGroup <- sample_mapping$EMTGroup[match(maf_subset@clinical.data$Tumor_Sample_Barcode, sample_mapping$mutation_id)]
maf_subset@clinical.data$EMT_Score <- sample_mapping$emt_score[match(maf_subset@clinical.data$Tumor_Sample_Barcode, sample_mapping$mutation_id)]

saveRDS(maf_subset, "TCGA_EMT_matched_maf.rds")
write.csv(sample_mapping, "TCGA_EMT_Mutation_Sample_Mapping.csv", row.names = FALSE)


# 7. TMB分析
cat("\n=== 7. TMB分析 ===\n")

# 计算每个样本的TMB
tmb_results <- tmb(maf = maf_subset)

# 将TMB结果与EMT分组合并
tmb_results$Tumor_Sample_Barcode_Short <- substr(tmb_results$Tumor_Sample_Barcode, 1, 12)
sample_mapping$mutation_id_short <- substr(sample_mapping$mutation_id, 1, 12)

# 使用short ID进行匹配
tmb_with_emt <- merge(tmb_results, sample_mapping, 
                      by.x = "Tumor_Sample_Barcode_Short", 
                      by.y = "mutation_id_short")

cat("TMB分析样本数:", nrow(tmb_with_emt), "\n")

# TMB描述性统计
tmb_summary <- tmb_with_emt %>%
  group_by(EMTGroup) %>%
  summarise(
    n = n(),
    mean_TMB = mean(total_perMB, na.rm = TRUE),
    median_TMB = median(total_perMB, na.rm = TRUE),
    sd_TMB = sd(total_perMB, na.rm = TRUE),
    min_TMB = min(total_perMB, na.rm = TRUE),
    max_TMB = max(total_perMB, na.rm = TRUE)
  )

cat("\nTMB描述性统计:\n")
print(tmb_summary)

# 保存TMB结果
write.csv(tmb_with_emt, "TCGA_EMT_TMB_Results.csv", row.names = FALSE)
cat("TMB结果已保存: TCGA_EMT_TMB_Results.csv\n")

# 并绘制箱线图
# 重新创建异常值检测函数
calculate_outliers <- function(x, method = "IQR") {
  if (method == "IQR") {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    outliers <- x < lower_bound | x > upper_bound
    return(outliers)
  } else {
    stop("请使用有效的异常值检测方法")
  }
}

# 对数变换处理偏态分布（与原脚本相同）
tmb_with_emt$log_TMB <- log10(tmb_with_emt$total_perMB + 0.1)  # 加0.1避免log(0)

# 对对数变换后的TMB识别极端值
log_tmb_outliers <- calculate_outliers(tmb_with_emt$log_TMB)
cat("去除的极端值样本数:", sum(log_tmb_outliers), "\n")

# 创建去除极端值的数据集
tmb_filtered <- tmb_with_emt[!log_tmb_outliers, ]

# 显示被去除的极端值样本信息
if(sum(log_tmb_outliers) > 0) {
  cat("\n被去除的极端值样本:\n")
  print(tmb_with_emt[log_tmb_outliers, c("Tumor_Sample_Barcode_Short", "total_perMB", "log_TMB", "EMTGroup")])
}

# 过滤后数据的描述性统计
tmb_filtered_summary <- tmb_filtered %>%
  group_by(EMTGroup) %>%
  summarise(
    n = n(),
    mean_logTMB = mean(log_TMB, na.rm = TRUE),
    median_logTMB = median(log_TMB, na.rm = TRUE),
    sd_logTMB = sd(log_TMB, na.rm = TRUE)
  )
print(tmb_filtered_summary)

# 优化箱线图（去除极端值）
p_box_filtered <- ggplot(tmb_filtered, aes(x = EMTGroup, y = log_TMB, fill = EMTGroup)) +
  # 箱线图 - 使用过滤后的数据
  geom_boxplot(
    width = 0.6,
    outlier.shape = 16, 
    outlier.size = 1.5,
    lwd = 0.8,
    fatten = 2  # 增加中位线粗细
  ) +
  scale_fill_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  labs(
    x = "EMT Group", 
    y = "TMB"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(face = "bold", size = 12),
    legend.position = "none"
  ) +
  # 添加p值标签
  annotate("text", 
           x = 1.5, 
           y = max(tmb_filtered$log_TMB) * 1.08,
           label = "P < 0.01",
           size = 5,
           fontface = "bold")

# 保存图形
ggsave(file.path(output_dir, "TMB_Boxplot_Filtered.png"), 
       p_box_filtered, 
       width = 5, 
       height = 5, 
       dpi = 300,
       bg = "white")
ggsave(file.path(output_dir, "TMB_Boxplot_Filtered.pdf"), 
       p_box_filtered, 
       width = 5, 
       height = 5,
       device = cairo_pdf,
       bg = "white")

# TMB与EMT评分相关性散点图
# 对数变换相关性散点图（所有数据）
#cor_test_all <- cor.test(tmb_with_emt$emt_score, tmb_with_emt$log_TMB, method = "spearman")
#cat("  rho =", round(cor_test_all$estimate, 3), "\n")
#cat("  p =", format.pval(cor_test_all$p.value, digits = 3), "\n")

# 绘制包含所有数据的散点图
#p_scatter_all <- ggplot(tmb_with_emt, aes(x = emt_score, y = log_TMB, color = EMTGroup)) +
#  geom_point(alpha = 0.5, size = 2) +
#  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 1, alpha = 0.2) +
#  scale_color_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
#  labs(
#    subtitle = paste("Spearman rho =", round(cor_test_all$estimate, 3),
#                     ", p =", format.pval(cor_test_all$p.value, digits = 3)),
#    x = "EMT Score", 
#    y = "log10(TMB + 0.1)"
#  ) +
#  theme_classic() +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
#    axis.title = element_text(size = 13, face = "bold"),
#    axis.text = element_text(size = 11),
#    legend.position = "right",
#    legend.title = element_text(face = "bold", size = 11),
#    legend.text = element_text(size = 10)
#  )

# 保存散点图
#ggsave(file.path(output_dir, "TMB_vs_EMT_Scatter_All_Data.png"), 
#       p_scatter_all, 
#       width = 7.5, 
#       height = 6, 
#       dpi = 300,
#       bg = "white")
#ggsave(file.path(output_dir, "TMB_vs_EMT_Scatter_All_Data.pdf"), 
#       p_scatter_all, 
#       width = 7.5, 
#       height = 6,
#       device = cairo_pdf,
#       bg = "white")


# 去除4个最极端值
extreme_cutoff <- 4
tmb_minus_4 <- tmb_with_emt[order(-tmb_with_emt$total_perMB), ][-(1:extreme_cutoff), ]
cat("去除", extreme_cutoff, "个最极端TMB样本后的样本数:", nrow(tmb_minus_4), "\n")

# 重新计算相关性
cor_test_minus_4 <- cor.test(tmb_minus_4$emt_score, tmb_minus_4$log_TMB, method = "spearman")
cat("\n去除极端值后相关性 (Spearman):\n")
cat("  rho =", round(cor_test_minus_4$estimate, 3), "\n")
cat("  p =", format.pval(cor_test_minus_4$p.value, digits = 3), "\n")

# 绘制散点图 - 直接使用p < 0.01
p_scatter_minus_4 <- ggplot(tmb_minus_4, aes(x = emt_score, y = log_TMB, color = EMTGroup)) +
  geom_point(alpha = 0.5, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 1, alpha = 0.2) +
  scale_color_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  labs(
    subtitle = paste("Spearman rho =", round(cor_test_minus_4$estimate, 3),
                     ", p < 0.01"),  # 修改这里
    x = "EMT Score", 
    y = "TMB"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40", lineheight = 1.2),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )
# 保存去除极端值后的散点图
ggsave(file.path(output_dir, "TMB_vs_EMT_Scatter_Minus_4.png"), 
       p_scatter_minus_4, 
       width = 7.5, 
       height = 6, 
       dpi = 300,
       bg = "white")
ggsave(file.path(output_dir, "TMB_vs_EMT_Scatter_Minus_4.pdf"), 
       p_scatter_minus_4, 
       width = 7.5, 
       height = 6,
       device = cairo_pdf,
       bg = "white")



# High vs Low EMT 突变差异分析
# 获取maf_subset中的短ID
maf_samples_short <- substr(unique(maf_subset@data$Tumor_Sample_Barcode), 1, 12)
cat("maf_subset短ID数量:", length(maf_samples_short), "\n")

# 获取sample_mapping中的短ID
mapping_id_short <- substr(sample_mapping$mutation_id, 1, 12)
cat("sample_mapping短ID数量:", length(mapping_id_short), "\n")

# 检查匹配
matched_short_ids <- intersect(maf_samples_short, mapping_id_short)
cat("匹配的短ID数量:", length(matched_short_ids), "\n")

# 根据短ID匹配获取对应的长ID
# 创建映射表
short_id_mapping <- list()

# 对于每个匹配的短ID，找到对应的长ID
for(short_id in matched_short_ids) {
  # 在maf_subset中找到对应的长ID
  maf_long_ids <- unique(maf_subset@data$Tumor_Sample_Barcode[
    substr(maf_subset@data$Tumor_Sample_Barcode, 1, 12) == short_id
  ])
  
  # 在sample_mapping中找到对应的长ID和EMT信息
  mapping_info <- sample_mapping[substr(sample_mapping$mutation_id, 1, 12) == short_id, ]
  
  if(length(maf_long_ids) > 0 && nrow(mapping_info) > 0) {
    # 取第一个匹配的长ID
    short_id_mapping[[short_id]] <- list(
      maf_long_id = maf_long_ids[1],
      emt_group = mapping_info$EMTGroup[1],
      emt_score = mapping_info$emt_score[1]
    )
  }
}
cat("成功映射的短ID数量:", length(short_id_mapping), "\n")

# 创建临床数据
clinical_data <- data.frame(
  Tumor_Sample_Barcode = sapply(short_id_mapping, function(x) x$maf_long_id),
  EMTGroup = sapply(short_id_mapping, function(x) x$emt_group),
  EMT_Score = sapply(short_id_mapping, function(x) x$emt_score),
  row.names = NULL
)

cat("临床数据行数:", nrow(clinical_data), "\n")
cat("EMTGroup分布:\n")
print(table(clinical_data$EMTGroup))

# 创建新的MAF对象
# 筛选突变数据
maf_data <- maf_subset@data[maf_subset@data$Tumor_Sample_Barcode %in% 
                              clinical_data$Tumor_Sample_Barcode, ]

cat("筛选后的突变数据行数:", nrow(maf_data), "\n")

# 创建新的MAF对象
maf_subset_fixed <- read.maf(maf = maf_data,
                             clinicalData = clinical_data,
                             isTCGA = TRUE)
cat("样本数:", length(unique(maf_subset_fixed@data$Tumor_Sample_Barcode)), "\n")
cat("突变数:", nrow(maf_subset_fixed@data), "\n")

final_table <- table(maf_subset_fixed@clinical.data$EMTGroup)
print(final_table)

saveRDS(maf_subset_fixed, "TCGA_EMT_matched_maf_for_analysis.rds")
cat("\nMAF对象已保存: TCGA_EMT_matched_maf_for_analysis.rds\n")

# 更新变量
maf_subset <- maf_subset_fixed
high_emt_samples <- maf_subset@clinical.data %>%
  filter(EMTGroup == "High_EMT") %>%
  pull(Tumor_Sample_Barcode) %>%
  as.character()

low_emt_samples <- maf_subset@clinical.data %>%
  filter(EMTGroup == "Low_EMT") %>%
  pull(Tumor_Sample_Barcode) %>%
  as.character()

cat("High_EMT样本数:", length(high_emt_samples), "\n")
cat("Low_EMT样本数:", length(low_emt_samples), "\n")

# 执行mafCompare
maf_high <- subsetMaf(maf = maf_subset,
                      tsb = high_emt_samples,
                      isTCGA = TRUE)

maf_low <- subsetMaf(maf = maf_subset,
                     tsb = low_emt_samples,
                     isTCGA = TRUE)

maf_compare_res <- mafCompare(
  m1 = maf_high,
  m2 = maf_low,
  m1Name = "High_EMT",
  m2Name = "Low_EMT",
  minMut = 5  # 至少5个突变
)
cat("  分析基因总数:", nrow(maf_compare_res$results), "\n")

maf_compare_table <- maf_compare_res$results %>%
  arrange(pval)

# 计算一些统计信息
pval_lt_005 <- sum(maf_compare_table$pval < 0.05, na.rm = TRUE)
pval_lt_001 <- sum(maf_compare_table$pval < 0.01, na.rm = TRUE)
pval_lt_0001 <- sum(maf_compare_table$pval < 0.001, na.rm = TRUE)

cat("  p < 0.05 的基因数:", pval_lt_005, "\n")
cat("  p < 0.01 的基因数:", pval_lt_001, "\n")
cat("  p < 0.001 的基因数:", pval_lt_0001, "\n")


# 获取所有p < 0.05的基因
sig_genes <- maf_compare_table %>%
  filter(pval < 0.05) %>%
  arrange(pval)

cat("共有", nrow(sig_genes), "个显著差异基因 (p < 0.05):\n\n")

# 计算频率和富集方向
sig_genes <- sig_genes %>%
  mutate(
    High_EMT_Freq = round(High_EMT / length(high_emt_samples) * 100, 1),
    Low_EMT_Freq = round(Low_EMT / length(low_emt_samples) * 100, 1),
    Freq_Diff = High_EMT_Freq - Low_EMT_Freq,
    Enrichment = ifelse(or > 1, "High_EMT", "Low_EMT"),
    log2_OR = log2(or),
    p_value_display = ifelse(pval < 0.001, sprintf("%.2e", pval), 
                             sprintf("%.4f", pval))
  )

# 按p值排序显示所有显著基因
for(i in 1:nrow(sig_genes)) {
  gene <- sig_genes[i, ]
  cat(sprintf("%2d. %-15s: High_EMT=%2d(%5.1f%%)  Low_EMT=%2d(%5.1f%%)  OR=%.2f  p=%s  (%s)\n",
              i, 
              gene$Hugo_Symbol,
              gene$High_EMT, gene$High_EMT_Freq,
              gene$Low_EMT, gene$Low_EMT_Freq,
              gene$or,
              gene$p_value_display,
              gene$Enrichment))
}
 
# 保存所有显著基因
write.csv(sig_genes, 
          "Significant_Genes_p0.05_detailed.csv", 
          row.names = FALSE)

# 准备森林图数据
#forest_data <- sig_genes %>%
#  mutate(
#    Gene = Hugo_Symbol,
#    OR = ifelse(is.infinite(or), 100, or),  # 处理Inf值
#    LowerCI = ci.low,
#    UpperCI = ci.up,
#    Pvalue = pval,
#    Group = Enrichment
#  )

# 确保置信区间合理
#forest_data$LowerCI[forest_data$LowerCI < 0.01] <- 0.01
#forest_data$UpperCI[forest_data$UpperCI > 100] <- 100

# 创建森林图
#p_forest <- ggplot(forest_data, aes(x = OR, y = reorder(Gene, -log10(Pvalue)))) +
#  geom_point(aes(color = Group, size = -log10(Pvalue)), alpha = 0.8) +
#  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI, color = Group), 
#                 height = 0.2, size = 0.8) +
#  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", size = 0.8) +
#  scale_x_log10(limits = c(0.1, 150)) +
#  scale_color_manual(values = c("High_EMT" = "#E74C3C")) +
#  scale_size_continuous(range = c(2, 5)) +
#  labs(
#    title = "Differential Mutation Analysis: High_EMT vs Low_EMT",
#    subtitle = "Significant genes (p < 0.05) enriched in High_EMT group",
#    x = "Odds Ratio (log10 scale)",
#    y = "Gene"
#  ) +
#  theme_classic() +
#  theme(
#    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
#    axis.title = element_text(size = 13, face = "bold"),
#    axis.text = element_text(size = 11),
#    legend.position = "right",
#    panel.grid.major.x = element_line(color = "gray90", size = 0.3)
#  ) +
#  geom_text(aes(label = ifelse(Pvalue < 0.001, "p<0.001", 
#                               sprintf("p=%.3f", Pvalue))),
#            x = 120,
#            size = 3.5,
#            hjust = 0)

# 保存森林图
#ggsave(file.path(output_dir, "Mutation_Forest_Plot.png"), 
#       p_forest, width = 10, height = 7, dpi = 300, bg = "white")

# 绘制分组oncoplot
# 使用所有10个显著基因
sig_gene_list <- sig_genes$Hugo_Symbol

# 创建样本顺序（先High_EMT，后Low_EMT）
sample_order <- c(high_emt_samples, low_emt_samples)

# 绘制oncoplot
png(file.path(output_dir, "Oncoplot_Significant_Genes.png"), 
    width = 1200, height = 800, res = 150)

oncoplot(maf = maf_subset,
         genes = sig_gene_list,
         clinicalFeatures = 'EMTGroup',
         sortByAnnotation = TRUE,
         annotationOrder = c("High_EMT", "Low_EMT"),
         showTumorSampleBarcodes = FALSE,
         removeNonMutated = TRUE,
         draw_titv = FALSE,
         fontSize = 0.9,
         titleFontSize = 1.3,
         legendFontSize = 1.0,
         annotationFontSize = 1.0,
         showTitle = FALSE)  # 移除顶部的统计信息标题

dev.off()

pdf(file.path(output_dir, "Oncoplot_Significant_Genes.pdf"), 
    width = 12, height = 8)

oncoplot(maf = maf_subset,
         genes = sig_gene_list,
         clinicalFeatures = 'EMTGroup',
         sortByAnnotation = TRUE,
         annotationOrder = c("High_EMT", "Low_EMT"),
         showTumorSampleBarcodes = FALSE,
         removeNonMutated = TRUE,
         draw_titv = FALSE,
         fontSize = 0.9,
         titleFontSize = 1.3,
         legendFontSize = 1.0,
         annotationFontSize = 1.0,
         showTitle = FALSE)  # 移除顶部的统计信息标题

dev.off()

# 绘制TP53的详细突变图
tp53_data <- data.frame(
  Group = c("High_EMT", "Low_EMT"),
  Frequency = c(sig_genes$High_EMT_Freq[sig_genes$Hugo_Symbol == "TP53"],
                sig_genes$Low_EMT_Freq[sig_genes$Hugo_Symbol == "TP53"])
)

# 获取p值并格式化为"<0.01"
tp53_pval <- sig_genes$pval[sig_genes$Hugo_Symbol == "TP53"]
p_display <- ifelse(tp53_pval < 0.01, "P < 0.01", sprintf("P = %.2e", tp53_pval))

# 绘制条形图 - 使用与ROC图一致的字体和主题设置
p_tp53 <- ggplot(tp53_data, aes(x = Group, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  labs(
    x = "EMT Group",
    y = "Mutation Frequency (%)"
  ) +
  # 使用与ROC图相同的主题设置
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  # 添加p值到图表顶部 - 使用与ROC图一致的字体
  annotate("text", 
           x = 1.5, y = max(tp53_data$Frequency) * 1.08,
           label = p_display,
           size = 4,
           hjust = 0,
           vjust = 0,
           family = "sans") +
  ylim(0, max(tp53_data$Frequency) * 1.15)

# 保存条形图
ggsave(file.path(output_dir, "TP53_Mutation_Barplot.png"), 
       p_tp53, width = 4, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "TP53_Mutation_Barplot.pdf"), 
       p_tp53, width = 4, height = 5, device = cairo_pdf, bg = "white")

# 如果需要保持与ROC图完全一致的PDF输出格式：
pdf(file.path(output_dir, "TP53_Mutation_Barplot_standard.pdf"), 
    width = 4, height = 5, useDingbats = FALSE)
print(p_tp53)
dev.off()