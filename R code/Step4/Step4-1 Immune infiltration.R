# Step 4-1: immune infiltration 
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/immune infiltration')
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(tidyr)
library(corrplot)

# 1. 加载数据和模型 - 使用原脚本的数据源
logtpm_matrix <- readRDS("../AUcell-Model/TCGA_PRAD_tpm.rds")
cat("TCGA表达数据维度:", dim(logtpm_matrix), "\n")
cat("TCGA样本数:", ncol(logtpm_matrix), "\n")

# 加载当前模型的基因
selected_genes_df <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/core_emt_genes.csv")
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

cat("选中的EMT基因:", paste(selected_genes, collapse = ", "), "\n")

# 样本去重
sample_ids <- colnames(logtpm_matrix)
base_patient_ids <- gsub("\\.", "-", substr(sample_ids, 1, 12))

cat("原始样本数:", length(sample_ids), "\n")
cat("唯一患者数:", length(unique(base_patient_ids)), "\n")

# 去重：每个患者保留第一个样本
selected_samples <- sample_ids[!duplicated(base_patient_ids)]
logtpm_unique <- logtpm_matrix[, selected_samples]

cat("去重后样本数:", length(selected_samples), "\n")
cat("去重后数据维度:", dim(logtpm_unique), "\n")

# 计算EMT评分
tcga_selected_genes <- logtpm_unique[rownames(logtpm_unique) %in% selected_genes, ]
cat("在TCGA中找到的EMT基因数:", nrow(tcga_selected_genes), "\n")

# 检查缺失基因
missing_genes <- setdiff(selected_genes, rownames(logtpm_unique))
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
  sample_id_tcga = selected_samples,
  sample_id_clean = gsub("\\.", "-", substr(selected_samples, 1, 12)),
  emt_score = tcga_emt_scores,
  EMTGroup = ifelse(tcga_emt_scores > emt_median, "High_EMT", "Low_EMT")
)

cat("EMT评分统计:\n")
print(summary(tcga_emt_scores))
cat("高EMT组:", sum(tcga_emt_data$EMTGroup == "High_EMT"), "\n")
cat("低EMT组:", sum(tcga_emt_data$EMTGroup == "Low_EMT"), "\n")

# 2. 免疫浸润分析
# 加载免疫特征基因集
immune_signatures <- read.csv("CellReports.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 解析免疫特征
cell_types <- immune_signatures[, 1]
available_immune_signatures <- list()

for(i in 1:nrow(immune_signatures)) {
  cell_type <- as.character(cell_types[i])
  genes <- as.character(immune_signatures[i, -1])
  genes <- genes[genes != "" & !is.na(genes)]
  genes <- unique(genes)
  
  available_genes <- genes[genes %in% rownames(logtpm_unique)]
  
  if(length(available_genes) >= 3) {
    available_immune_signatures[[cell_type]] <- available_genes
  }
}

cat("可用的免疫细胞类型数:", length(available_immune_signatures), "\n")

# 执行ssGSEA
ssgsea_results <- gsva(as.matrix(logtpm_unique), 
                       available_immune_signatures,
                       method = "ssgsea",
                       kcdf = "Gaussian",
                       parallel.sz = 4,
                       ssgsea.norm = TRUE)

cat("免疫细胞类型数:", nrow(ssgsea_results), "\n")
cat("样本数:", ncol(ssgsea_results), "\n")

# 保存结果
saveRDS(ssgsea_results, "ssgsea_immune_results_current_model.rds")

# 数据整合
immune_df <- as.data.frame(t(ssgsea_results))
immune_df$id <- rownames(immune_df)
immune_df$id_clean <- gsub("\\.", "-", substr(immune_df$id, 1, 12))

# 合并数据
common_samples <- intersect(immune_df$id_clean, tcga_emt_data$sample_id_clean)
cat("共同样本数:", length(common_samples), "\n")

final_immune <- immune_df[immune_df$id_clean %in% common_samples, ]
final_emt <- tcga_emt_data[tcga_emt_data$sample_id_clean %in% common_samples, ]

# 确保顺序一致
final_immune <- final_immune[match(common_samples, final_immune$id_clean), ]
final_emt <- final_emt[match(common_samples, final_emt$sample_id_clean), ]

analysis_data <- cbind(final_immune, EMTGroup = final_emt$EMTGroup)

cat("最终分析数据样本数:", nrow(analysis_data), "\n")

# 3. 免疫浸润差异分析
immune_cell_types <- colnames(analysis_data)[!colnames(analysis_data) %in% c("id", "id_clean", "EMTGroup")]

# 执行t检验
immune_results <- data.frame(
  cell_type = immune_cell_types,
  p_value = NA,
  mean_high = NA,
  mean_low = NA,
  fold_change = NA,
  t_statistic = NA
)

for(i in 1:length(immune_cell_types)) {
  cell_type <- immune_cell_types[i]
  high_group <- analysis_data[analysis_data$EMTGroup == "High_EMT", cell_type]
  low_group <- analysis_data[analysis_data$EMTGroup == "Low_EMT", cell_type]
  
  t_test <- t.test(high_group, low_group)
  immune_results$p_value[i] <- t_test$p.value
  immune_results$mean_high[i] <- mean(high_group)
  immune_results$mean_low[i] <- mean(low_group)
  immune_results$fold_change[i] <- mean(high_group) / mean(low_group)
  immune_results$t_statistic[i] <- t_test$statistic
}

# 多重检验校正
immune_results$fdr <- p.adjust(immune_results$p_value, method = "fdr")
immune_results$log2_fc <- log2(immune_results$fold_change)

# 筛选显著结果
significant_immune <- immune_results[immune_results$fdr < 0.05, ]
significant_immune <- significant_immune[order(abs(significant_immune$log2_fc), decreasing = TRUE), ]

cat("显著差异的免疫细胞类型数 (FDR < 0.05):", nrow(significant_immune), "\n")

# 保存结果
write.csv(immune_results, "immune_cell_differences_current_model.csv", row.names = FALSE)
write.csv(significant_immune, "significant_immune_cells_current_model.csv", row.names = FALSE)

# 4.可视化
figure_dir <- '/home/datahup/xzp/paper/PCa/fig/immune infiltration'

# 准备绘图数据（包含所有28个免疫细胞）
plot_data <- analysis_data %>%
  select(-id, -id_clean) %>%
  pivot_longer(cols = -EMTGroup, 
               names_to = "Immune_Cell_Type", 
               values_to = "Estimating_Score") %>%
  rename(RiskGroup = EMTGroup) %>%
  mutate(RiskGroup = ifelse(RiskGroup == "High_EMT", "High", "Low"))

# 获取所有免疫细胞类型（不筛选，显示全部）
all_cells <- colnames(analysis_data)[!colnames(analysis_data) %in% c("id", "id_clean", "EMTGroup")]
cat("总共免疫细胞类型数:", length(all_cells), "\n")
cat("免疫细胞类型:", paste(all_cells, collapse = ", "), "\n")

# 使用所有细胞进行可视化
plot_data_all <- plot_data[plot_data$Immune_Cell_Type %in% all_cells, ]

# 生成小提琴图（包含所有28个免疫细胞）
p_all <- ggplot(plot_data_all, aes(x = reorder(Immune_Cell_Type, -Estimating_Score), 
                                   y = Estimating_Score, fill = RiskGroup)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.7,
              width = 1.2, trim = TRUE,
              color = NA) +
  geom_boxplot(width = 0.25, show.legend = FALSE,
               position = position_dodge(0.9),
               color = 'black', alpha = 0.8,
               outlier.shape = 21, outlier.size = 0.8,
               outlier.fill = "gray", outlier.alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, 
                                   color = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 10),
        legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = c("High" = "#FF6B6B", "Low" = "#4ECDC4"),
                    name = "EMT Group",
                    labels = c("High (High EMT)", "Low (Low EMT)")) +
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup, label = ..p.signif..), 
                     method = "wilcox.test", 
                     label.y = max(plot_data_all$Estimating_Score) * 1.05,
                     size = 2.8,
                     bracket.size = 0.3,
                     tip.length = 0.02)

# 保存包含所有免疫细胞的图
ggsave(filename = file.path(figure_dir, "ssGSEA_all_immune_cells_current_model.pdf"), 
       plot = p_all, width = 18, height = 10, dpi = 300)

# EMT评分与免疫细胞相关性分析
cor_data <- analysis_data %>%
  select(-id, -id_clean, -EMTGroup) %>%
  cbind(emt_score = final_emt$emt_score)

# 计算完整的Spearman相关性矩阵（包含所有免疫细胞和EMT评分）
Cor_imm <- cor(cor_data, method = "spearman")
cat("相关性矩阵维度:", dim(Cor_imm), "\n")

# 计算显著性
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# 计算p值矩阵
p.mat <- cor.mtest(cor_data)

# 提取EMT评分与其他免疫细胞的相关性结果
emt_correlations <- Cor_imm["emt_score", -which(colnames(Cor_imm) == "emt_score")]
emt_pvalues <- p.mat["emt_score", -which(colnames(p.mat) == "emt_score")]

# 创建相关性结果数据框
cor_results <- data.frame(
  immune_cell = names(emt_correlations),
  correlation = emt_correlations,
  p_value = emt_pvalues
)

# 计算FDR校正
cor_results$fdr <- p.adjust(cor_results$p_value, method = "fdr")

# 按相关性绝对值排序
cor_results <- cor_results[order(abs(cor_results$correlation), decreasing = TRUE), ]

cat("EMT评分与免疫细胞相关性分析完成\n")
cat("显著相关的免疫细胞数 (p < 0.05):", sum(cor_results$p_value < 0.05), "\n")
cat("显著相关的免疫细胞数 (FDR < 0.05):", sum(cor_results$fdr < 0.05), "\n")

# 保存相关性结果
write.csv(cor_results, "emt_immune_correlation_results_current_model.csv", row.names = FALSE)

# 生成完整的相关性热图
library(corrplot)

# 设置颜色
col <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                          "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))

# 绘制完整的相关性热图
pdf(file = file.path(figure_dir, "emt_immune_correlation_current_model.pdf"), width = 14, height = 12)

corrplot(Cor_imm, 
         type = "upper",
         title = 'Spearman Correlation between EMT Score and Immune Cells', 
         mar = c(0, 0, 1.2, 0),
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.8,
         cl.cex = 0.8,
         addCoef.col = "black", 
         number.cex = 0.6,
         method = "ellipse",
         col = col(100),
         diag = FALSE,
         order = "original")

dev.off()

# 生成综合的EMT基因相关性热图
# 准备综合相关性数据
comprehensive_cor_data <- gene_immune_data[, c(selected_genes, immune_cell_types)]

# 计算综合相关性矩阵
cor_comprehensive <- cor(comprehensive_cor_data, method = "spearman")

# 计算显著性
p.mat_comprehensive <- cor.mtest(comprehensive_cor_data)

# 提取EMT基因与免疫细胞的相关性子矩阵
emt_immune_cor <- cor_comprehensive[selected_genes, immune_cell_types]
emt_immune_p <- p.mat_comprehensive[selected_genes, immune_cell_types]

# 方法1: 使用pheatmap但修复保存问题
tryCatch({
  # 设置颜色
  col_heatmap <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)
  
  # 创建显著性标记矩阵
  significance_matrix <- matrix("", nrow = nrow(emt_immune_cor), ncol = ncol(emt_immune_cor))
  significance_matrix[emt_immune_p < 0.001] <- "***"
  significance_matrix[emt_immune_p < 0.01 & emt_immune_p >= 0.001] <- "**"
  significance_matrix[emt_immune_p < 0.05 & emt_immune_p >= 0.01] <- "*"
  
  # 使用pheatmap保存
  p <- pheatmap(emt_immune_cor,
                scale = "none",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_colnames = TRUE,
                show_rownames = TRUE,
                color = col_heatmap,
                display_numbers = significance_matrix,
                number_color = "black",
                fontsize_row = 10,
                fontsize_col = 8,
                main = "Spearman Correlation: EMT Genes vs Immune Cells")
  
  # 手动保存图形
  pdf(file = file.path(figure_dir, "emt_genes_immune_correlation_comprehensive_current_model.pdf"), 
      width = 16, height = 12)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()})

