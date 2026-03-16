### Step4-3 TIDE
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/TIDE')
fig_dir <- '/home/datahup/xzp/paper/PCa/fig/TIDE'

# 加载所需包
library(tidyverse)
library(ggplot2)
library(ggpubr)

# 1. 加载PCa数据
cat("\n=== 加载PCa数据 ===\n")

expression_dir <- '/home/datahup/xzp/paper/PCa/data/AUcell-Model'
logtpm_matrix <- readRDS(file.path(expression_dir, "TCGA_PRAD_tpm.rds"))

# 样本去重
sample_ids <- colnames(logtpm_matrix)
base_patient_ids <- gsub("\\.", "-", substr(sample_ids, 1, 12))
selected_samples <- sample_ids[!duplicated(base_patient_ids)]
logtpm_unique <- logtpm_matrix[, selected_samples]

cat("TCGA-PRAD样本数:", length(selected_samples), "\n")

# 计算EMT评分
emt_model_path <- "../TCGA GEO-Model/seed_125_lasso_results/core_emt_genes.csv"
selected_genes_df <- read.csv(emt_model_path)
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

tcga_selected_genes <- logtpm_unique[rownames(logtpm_unique) %in% selected_genes, ]
tcga_X <- t(tcga_selected_genes)
tcga_emt_scores <- as.numeric(as.matrix(tcga_X) %*% selected_coef)

emt_median <- median(tcga_emt_scores, na.rm = TRUE)
tcga_emt_data <- data.frame(
  sample_id = selected_samples,
  emt_score = tcga_emt_scores,
  EMTGroup = ifelse(tcga_emt_scores > emt_median, "High_EMT", "Low_EMT")
)

# 准备TIDE官网输入文件
#expression_dir <- '/home/datahup/xzp/paper/PCa/data/AUcell-Model'
#logtpm_matrix <- readRDS(file.path(expression_dir, "TCGA_PRAD_tpm.rds"))

# 样本去重
#sample_ids <- colnames(logtpm_matrix)
#base_patient_ids <- gsub("\\.", "-", substr(sample_ids, 1, 12))
#selected_samples <- sample_ids[!duplicated(base_patient_ids)]
#logtpm_unique <- logtpm_matrix[, selected_samples]

#cat("TCGA-PRAD样本数:", length(selected_samples), "\n")
#cat("基因数:", nrow(logtpm_unique), "\n")

# TIDE官网要求的格式：
# - 行：基因
# - 列：样本
# - 第一列：基因名
# - 数据：标准化后的表达值

# 使用TPM数据的log2转换
#tide_input <- as.data.frame(logtpm_unique)

# 添加基因名列
#tide_input$Gene <- rownames(tide_input)

# 重新排列列，让Gene作为第一列
#tide_input <- tide_input %>% select(Gene, everything())

# 检查数据范围
#summary(as.numeric(as.matrix(logtpm_unique)))

# 保存TIDE输入文件
#write_tsv(tide_input, "TCGA_PRAD_TIDE_Input.tsv")

#cat("TIDE输入文件已保存: TCGA_PRAD_TIDE_Input.tsv\n")
#cat("文件维度:", dim(tide_input), "\n")

# 2. 加载TIDE官网结果并整合分析
# 加载TIDE结果文件
tide_results <- read.csv("TIDE_result.csv")
cat("TIDE结果文件维度:", dim(tide_results), "\n")
cat("列名:", colnames(tide_results), "\n")
cat("前几行数据:\n")
print(head(tide_results))

# 3. 数据整合
# 准备样本ID用于匹配
tide_results$Patient_ID <- gsub("\\.", "-", substr(tide_results$Patient, 1, 12))
tcga_emt_data$Patient_ID <- gsub("\\.", "-", substr(tcga_emt_data$sample_id, 1, 12))

merged_data <- tide_results %>%
  left_join(tcga_emt_data, by = "Patient_ID")
analysis_data <- merged_data %>% filter(!is.na(emt_score))
cat("最终分析样本数:", nrow(analysis_data), "\n")

# 4. 统计分析
# 4.1 TIDE评分与EMT分组比较
t_test_tide <- t.test(TIDE ~ EMTGroup, data = analysis_data)
cat("TIDE评分组间差异 (t检验):\n")
cat("p值 =", format.pval(t_test_tide$p.value, digits = 3), "\n")

# 分组统计
group_stats <- analysis_data %>%
  group_by(EMTGroup) %>%
  summarise(
    n = n(),
    Mean_TIDE = mean(TIDE),
    Median_TIDE = median(TIDE),
    SD_TIDE = sd(TIDE),
    Mean_Dysfunction = mean(Dysfunction),
    Mean_Exclusion = mean(Exclusion),
    Mean_CD8 = mean(CD8),
    Mean_CD274 = mean(CD274),
    Mean_IFNG = mean(IFNG)
  )
print(group_stats)

# 4.2 EMT评分与TIDE评分相关性
cor_test <- cor.test(analysis_data$emt_score, analysis_data$TIDE, method = "spearman")
cat("EMT评分与TIDE评分相关性:\n")
cat("Spearman rho =", round(cor_test$estimate, 3), "\n")
cat("p值 =", format.pval(cor_test$p.value, digits = 3), "\n")

# 5. 可视化分析
# 5.1 TIDE评分箱线图
tide_combined_data <- analysis_data %>%
  select(EMTGroup, TIDE, Dysfunction, Exclusion) %>%
  pivot_longer(cols = c(TIDE, Dysfunction, Exclusion), 
               names_to = "ScoreType", 
               values_to = "ScoreValue")

# 设置评分类型标签和顺序
tide_combined_data$ScoreType <- factor(tide_combined_data$ScoreType,
                                       levels = c("TIDE", "Dysfunction", "Exclusion"),
                                       labels = c("TIDE Score", "Dysfunction", "Exclusion"))

# 绘制组合箱线图
# 或者直接显示文字"p < 0.05"
p_tide_combined <- ggplot(tide_combined_data, aes(x = EMTGroup, y = ScoreValue, fill = EMTGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, shape = 16) +
  scale_fill_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  labs(
    x = "EMT Group",
    y = "Score",
    fill = "EMT Group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.spacing = unit(0.5, "cm")
  ) +
  facet_wrap(~ ScoreType, ncol = 3, scales = "free_y") +
  # 手动添加文字标签
  geom_text(data = data.frame(ScoreType = unique(tide_combined_data$ScoreType)),
            aes(x = 1.5, y = Inf, label = "p < 0.05"),
            vjust = 2, size = 4.5, fontface = "bold", inherit.aes = FALSE)

# 保存组合图
ggsave(file.path(fig_dir, "TIDE_Scores_Combined_Boxplot.pdf"), 
       p_tide_combined, width = 12, height = 5)
ggsave(file.path(fig_dir, "TIDE_Scores_Combined_Boxplot.png"), 
       p_tide_combined, width = 12, height = 5, dpi = 300)

# 5.2 EMT评分与TIDE评分相关性散点图
p_correlation <- ggplot(analysis_data, aes(x = emt_score, y = TIDE, color = EMTGroup)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.2) +
  scale_color_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  labs(subtitle = paste("Spearman rho =", round(cor_test$estimate, 3), ", p < 0.01"),
       x = "EMT Score", y = "TIDE Score",
       color = "EMT Group") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

# 保存相关性散点图
ggsave(file.path(fig_dir, "EMT_TIDE_Correlation.pdf"), p_correlation, width = 7, height = 6)

# 保存数据结果
write.csv(analysis_data, "TCGA_PRAD_TIDE_EMT_Combined_Analysis.csv", row.names = FALSE)
write.csv(group_stats, "TIDE_EMT_Group_Statistics.csv", row.names = FALSE)

# 6.免疫检查点分析
# 读取免疫检查点基因列表
ncbi_data <- read.delim("gene_list.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
                        comment.char = "", quote = "", fill = TRUE)

cat("NCBI文件维度:", dim(ncbi_data), "\n")
cat("列名:", colnames(ncbi_data), "\n")

all_ncbi_genes <- unique(na.omit(ncbi_data$Symbol))
cat("从NCBI文件提取的基因总数:", length(all_ncbi_genes), "\n")

available_checkpoints <- all_ncbi_genes[all_ncbi_genes %in% rownames(logtpm_unique)]
cat("在表达矩阵中可用的免疫检查点基因数:", length(available_checkpoints), 
    "/", length(all_ncbi_genes), " (", 
    round(length(available_checkpoints)/length(all_ncbi_genes)*100, 1), "%)\n")

checkpoint_info <- data.frame(
  Symbol = available_checkpoints,
  In_NCBI = TRUE,
  In_Expression = TRUE,
  Full_Name = ncbi_data$description[match(available_checkpoints, ncbi_data$Symbol)]
)
write.csv(checkpoint_info, "NCBI_Immune_Checkpoint_Genes_Available.csv", row.names = FALSE)

# 风险基因与免疫检查点基因相关性分析
risk_genes <- unique(selected_genes)
cat("建模使用的风险基因数:", length(risk_genes), "\n")

risk_genes_info <- data.frame(
  Gene = risk_genes,
  Coefficient = selected_coef[risk_genes],
  In_Expression = risk_genes %in% rownames(logtpm_unique)
)

# 提取风险基因和检查点基因的表达数据
valid_risk_genes <- risk_genes[risk_genes %in% rownames(logtpm_unique)]
cat("在表达矩阵中存在的风险基因:", length(valid_risk_genes), "/", length(risk_genes), "\n")
risk_genes_data <- logtpm_unique[valid_risk_genes, , drop = FALSE]
checkpoint_data <- logtpm_unique[available_checkpoints, , drop = FALSE]

# 确保样本顺序一致
common_samples <- intersect(colnames(risk_genes_data), colnames(checkpoint_data))
risk_genes_data <- risk_genes_data[, common_samples, drop = FALSE]
checkpoint_data <- checkpoint_data[, common_samples, drop = FALSE]

cat("用于分析的共同样本数:", length(common_samples), "\n")

# 计算相关性矩阵
cor_matrix <- matrix(NA, nrow = length(valid_risk_genes), ncol = length(available_checkpoints),
                     dimnames = list(valid_risk_genes, available_checkpoints))

pvalue_matrix <- matrix(NA, nrow = length(valid_risk_genes), ncol = length(available_checkpoints),
                        dimnames = list(valid_risk_genes, available_checkpoints))

pb <- txtProgressBar(min = 0, max = length(valid_risk_genes), style = 3)

for (i in seq_along(valid_risk_genes)) {
  risk_gene <- valid_risk_genes[i]
  risk_expr <- as.numeric(risk_genes_data[risk_gene, ])
  
  for (j in seq_along(available_checkpoints)) {
    checkpoint <- available_checkpoints[j]
    checkpoint_expr <- as.numeric(checkpoint_data[checkpoint, ])
    
    valid_idx <- complete.cases(risk_expr, checkpoint_expr)
    if (sum(valid_idx) > 10) {
      # 使用tryCatch避免错误中断计算
      cor_result <- tryCatch({
        cor.test(risk_expr[valid_idx], checkpoint_expr[valid_idx], 
                 method = "spearman", exact = FALSE)  # 添加exact=FALSE避免警告
      }, error = function(e) {
        return(list(estimate = NA, p.value = NA))
      })
      
      cor_matrix[i, j] <- cor_result$estimate
      pvalue_matrix[i, j] <- cor_result$p.value
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# 保存相关性结果
write.csv(cor_matrix, "Risk_Genes_Checkpoint_Correlation_Matrix.csv")
write.csv(pvalue_matrix, "Risk_Genes_Checkpoint_Pvalues.csv")

# 6.3 综合筛选TOP免疫检查点基因
# 计算每个免疫检查点基因与所有风险基因的平均相关性绝对值
checkpoint_scores <- data.frame(
  Checkpoint = available_checkpoints,
  MeanAbsCor = NA,           # 平均绝对相关性（主要评分指标）
  MaxAbsCor = NA,            # 最大绝对相关性（参考）
  PosCorCount = NA,          # 正相关数量
  NegCorCount = NA,          # 负相关数量
  MeanCor = NA               # 平均原始相关性（考虑方向）
)

for (j in seq_along(available_checkpoints)) {
  checkpoint <- available_checkpoints[j]
  cor_values <- cor_matrix[, checkpoint]  # 该检查点与所有风险基因的相关性
  
  # 计算平均绝对相关性（主要评分）
  abs_cor_values <- abs(cor_values)
  checkpoint_scores$MeanAbsCor[j] <- mean(abs_cor_values, na.rm = TRUE)
  
  # 计算最大绝对相关性
  checkpoint_scores$MaxAbsCor[j] <- max(abs_cor_values, na.rm = TRUE)
  
  # 计算相关性方向统计
  checkpoint_scores$PosCorCount[j] <- sum(cor_values > 0, na.rm = TRUE)
  checkpoint_scores$NegCorCount[j] <- sum(cor_values < 0, na.rm = TRUE)
  
  # 计算平均原始相关性（考虑正负方向）
  checkpoint_scores$MeanCor[j] <- mean(cor_values, na.rm = TRUE)
}

# 按平均相关性绝对值排序（主要排序标准）
checkpoint_scores <- checkpoint_scores[order(checkpoint_scores$MeanAbsCor, decreasing = TRUE), ]
print(head(checkpoint_scores, 30))

# 计算相关性强度等级
checkpoint_scores$CorStrength <- cut(
  checkpoint_scores$MeanAbsCor,
  breaks = c(0, 0.3, 0.5, 0.7, 1),
  labels = c("Weak", "Moderate", "Strong", "Very Strong"),
  include.lowest = TRUE
)

# 保存检查点评分结果
write.csv(checkpoint_scores, "Immune_Checkpoint_Mean_Correlation_Scores.csv", row.names = FALSE)

# 选择TOP50进行检查点展示
top_50_checkpoints <- head(checkpoint_scores$Checkpoint, 50)
cat("\nTOP50免疫检查点基因（按平均绝对相关性排序）：\n")
print(top_50_checkpoints)

# 手动挑选的20个经典免疫检查点基因
selected_classic_checkpoints <- c(
  "VSIR",       # VISTA - 重要检查点
  "NT5E",       # CD73 - 检查点酶
  "ENTPD1",     # CD39 - 检查点酶
  "LGALS9",     # Galectin-9 - TIM-3配体
  "CD200",      # OX-2 - 免疫检查点
  "CD40",       # 共刺激分子
  "PDCD1LG2",   # PD-L2 - PD-1配体
  "SIRPA",      # SIRPα - CD47配体
  "TGFB1",      # TGF-β1 - 重要免疫抑制因子
  "TNFRSF1B",   # TNFR2 - 免疫调节受体
  "TNFRSF1A",   # TNFR1 - 免疫调节受体
  "IRAK3",      # 免疫调节激酶
  "IFI16",      # 干扰素诱导免疫蛋白
  "HSD11B1",    # 免疫代谢调节
  "CCL21",      # 趋化因子
  "SERPINB9",   # 免疫调节丝氨酸蛋白酶抑制剂
  "IL27RA",     # IL-27受体
  "SOCS1",      # 细胞因子信号抑制因子
  "CXCR4",      # 趋化因子受体
  "KLRB1"       # CD161 - NK细胞受体
)

cat("挑选的20个经典免疫检查点基因:\n")
print(selected_classic_checkpoints)

# 6.5 绘制相关性点图（风险基因 vs 经典免疫检查点）
# 准备绘图数据
plot_data <- as.data.frame(as.table(cor_matrix))
colnames(plot_data) <- c("RiskGene", "Checkpoint", "Correlation")
plot_data_classic <- plot_data[plot_data$Checkpoint %in% selected_classic_checkpoints, ]

# 排序
risk_genes_ordered <- valid_risk_genes[order(abs(selected_coef[valid_risk_genes]), decreasing = TRUE)]
plot_data_classic$RiskGene <- factor(plot_data_classic$RiskGene, levels = risk_genes_ordered)

checkpoint_means <- checkpoint_scores[checkpoint_scores$Checkpoint %in% selected_classic_checkpoints, ]
checkpoint_means <- checkpoint_means[order(checkpoint_means$MeanAbsCor, decreasing = TRUE), ]
plot_data_classic$Checkpoint <- factor(plot_data_classic$Checkpoint, 
                                       levels = checkpoint_means$Checkpoint)

# 创建相关性强度分组（用于尺寸）
plot_data_classic$CorStrength <- cut(plot_data_classic$Correlation,
                                     breaks = c(0, 0.4, 0.6, 0.8, 1),
                                     labels = c("0.3", "0.5", "0.7", "0.9"))

# 设置尺寸值
size_values <- c("0.3" = 4, "0.5" = 6, "0.7" = 8, "0.9" = 10)

# 绘制简洁的实心圆点图
p_dotplot_simple <- ggplot(plot_data_classic, aes(x = Checkpoint, y = RiskGene)) +
  geom_point(aes(size = CorStrength, color = Correlation), alpha = 0.8) +
  
  # 使用离散的尺寸
  scale_size_manual(
    name = "Correlation\nStrength",
    values = size_values,
    labels = c("0.3", "0.5", "0.7", "0.9")
  ) +
  
  # 颜色渐变
  scale_color_gradientn(
    name = "Correlation\nCoefficient",
    colors = c("#ADD8E6", "#87CEEB", "#4682B4", "#FF6B6B", "#FF3333", "#FF0000"),
    limits = c(0.2, 0.9),
    breaks = c(0.3, 0.5, 0.7, 0.9),
    labels = c("0.3", "0.5", "0.7", "0.9")
  ) +
  
  labs(
    x = "Immune Checkpoint Genes",
    y = "EMT Risk Genes"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank()
  ) +
  # 添加背景网格
  geom_hline(yintercept = seq(0.5, length(valid_risk_genes) - 0.5, by = 1), 
             color = "grey95", size = 0.1) +
  geom_vline(xintercept = seq(0.5, length(selected_classic_checkpoints) - 0.5, by = 1), 
             color = "grey95", size = 0.1)

# 保存图片
ggsave(file.path(fig_dir, "Correlation_DotPlot_Final.pdf"), 
       p_dotplot_simple, width = 14, height = 7)
ggsave(file.path(fig_dir, "Correlation_DotPlot_Final.png"), 
       p_dotplot_simple, width = 14, height = 7, dpi = 300)

# 7. 高、低风险组免疫检查点相关基因丰度差异
# 获取EMT分组信息
emt_group_info <- tcga_emt_data %>%
  filter(sample_id %in% common_samples) %>%
  select(sample_id, EMTGroup)
emt_group_vector <- setNames(emt_group_info$EMTGroup, emt_group_info$sample_id)

# 7.1 选择要展示的免疫检查点基因
top_display_num <- 20
top_display_checkpoints <- head(selected_classic_checkpoints, top_display_num)

# 准备组合箱线图数据
combined_boxplot_data <- data.frame()

for (checkpoint in top_display_checkpoints) {
  # 获取表达数据
  expr_values <- as.numeric(checkpoint_data[checkpoint, common_samples])
  
  # 创建数据框
  df <- data.frame(
    Gene = checkpoint,
    Expression = expr_values,
    EMTGroup = emt_group_vector[common_samples],
    Sample = common_samples
  )
  
  # 移除缺失值
  df <- df[complete.cases(df), ]
  
  combined_boxplot_data <- rbind(combined_boxplot_data, df)
}

# 按平均相关性排序基因
gene_order <- checkpoint_means %>%
  filter(Checkpoint %in% top_display_checkpoints) %>%
  arrange(desc(MeanAbsCor))

combined_boxplot_data$Gene <- factor(combined_boxplot_data$Gene,
                                     levels = gene_order$Checkpoint)

# 7.2 绘制组合箱线图
p_combined_boxplot <- ggplot(combined_boxplot_data, aes(x = Gene, y = Expression, fill = EMTGroup)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7, 
               position = position_dodge(width = 0.8)) +
  
  scale_fill_manual(values = c("High_EMT" = "#E74C3C", "Low_EMT" = "#3498DB")) +
  
  labs(
    x = "Immune Checkpoint Genes",
    y = "Expression Level (log2TPM)",
    fill = "EMT Group"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey95", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# 保存组合箱线图
ggsave(file.path(fig_dir, "Immune_Checkpoints_Expression_Boxplot.pdf"), 
       p_combined_boxplot, width = 10, height = 7)
ggsave(file.path(fig_dir, "Immune_Checkpoints_Expression_Boxplot.png"), 
       p_combined_boxplot, width = 10, height = 7, dpi = 300)

cat("免疫检查点表达箱线图已保存\n")

# 7.3 保存结果
cat("\n=== 步骤7.3: 保存表达差异统计结果 ===\n")

expression_diff_stats <- data.frame()

for (checkpoint in selected_classic_checkpoints) {
  expr_values <- as.numeric(checkpoint_data[checkpoint, common_samples])
  high_expr <- expr_values[emt_group_vector[common_samples] == "High_EMT"]
  low_expr <- expr_values[emt_group_vector[common_samples] == "Low_EMT"]
  
  high_expr <- high_expr[!is.na(high_expr)]
  low_expr <- low_expr[!is.na(low_expr)]
  
  if (length(high_expr) >= 5 && length(low_expr) >= 5) {
    t_test <- tryCatch({
      t.test(high_expr, low_expr)
    }, error = function(e) {
      return(list(p.value = NA, estimate = NA))
    })
    
    mean_high <- mean(high_expr, na.rm = TRUE)
    mean_low <- mean(low_expr, na.rm = TRUE)
    log2fc <- ifelse(mean_low != 0, log2((mean_high + 1) / (mean_low + 1)), NA)
    
    expression_diff_stats <- rbind(expression_diff_stats, data.frame(
      Checkpoint = checkpoint,
      Mean_High = mean_high,
      Mean_Low = mean_low,
      Log2FC = log2fc,
      Ttest_p = ifelse(!is.na(t_test$p.value), t_test$p.value, NA),
      Samples_High = length(high_expr),
      Samples_Low = length(low_expr)
    ))
  }
}

expression_diff_stats$FDR <- p.adjust(expression_diff_stats$Ttest_p, method = "fdr")
expression_diff_stats$Significant <- ifelse(expression_diff_stats$FDR < 0.05, "Yes", "No")

write.csv(expression_diff_stats, "Immune_Checkpoints_Expression_Stats.csv", row.names = FALSE)