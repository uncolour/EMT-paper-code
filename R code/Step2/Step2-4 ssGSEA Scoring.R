# Step 2-4: ssGSEA Score
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/AUcell')
library(Seurat)
library(ggplot2)
library(dplyr)
library(GSVA)
library(msigdbr)

# 1. 加载数据
CancerCell <- readRDS("CancerCell_with_EMT_Analysis.rds")
de_genes <- read.csv("DEGs_High_vs_Low_EMT.csv")
cat("细胞数量:", ncol(CancerCell), "\n")
cat("差异基因数量:", nrow(de_genes), "\n")

emt_geneset <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)
cat("EMT基因集大小:", length(emt_geneset), "\n")

# 2. ssGSEA评分
expr_matrix <- as.matrix(GetAssayData(CancerCell, assay = "RNA", slot = "data"))
gene_sets <- list(EMT_Signature = emt_geneset)

ssgsea_scores <- gsva(expr_matrix, 
                      gene_sets,
                      method = "ssgsea",
                      kcdf = "Gaussian",
                      parallel.sz = 8,
                      verbose = TRUE)

ssgsea_df <- as.data.frame(t(ssgsea_scores))
rownames(ssgsea_df) <- colnames(CancerCell)
CancerCell <- AddMetaData(CancerCell, metadata = ssgsea_df)

# 3. 筛选与EMT显著相关的差异基因
emt_scores <- CancerCell$EMT_Signature

# 只考虑显著差异基因 (adj.p < 0.05)
sig_genes <- de_genes %>% 
  filter(p_val_adj < 0.05) %>% 
  pull(gene) %>% 
  intersect(rownames(expr_matrix))

# 计算每个基因与EMT评分的相关性
cor_results <- data.frame(
  gene = sig_genes,
  cor_coef = NA,
  p_value = NA,
  fdr = NA
)

for (i in 1:length(sig_genes)) {
  gene_expr <- as.numeric(expr_matrix[sig_genes[i], ])
  if (sd(gene_expr) > 0) {
    cor_test <- cor.test(gene_expr, emt_scores, method = "spearman")
    cor_results$cor_coef[i] <- cor_test$estimate
    cor_results$p_value[i] <- cor_test$p.value
  }
}
cor_results$fdr <- p.adjust(cor_results$p_value, method = "fdr")

# 筛选显著相关的基因 (|cor| > 0.3 & FDR < 0.05)
emt_related_genes <- cor_results %>%
  filter(abs(cor_coef) > 0.3 & fdr < 0.05) %>%
  arrange(desc(abs(cor_coef)))
cat("与EMT显著相关的基因数量:", nrow(emt_related_genes), "\n")

write.csv(emt_related_genes, "EMT_Related_Genes.csv", row.names = FALSE)

# 4. 可视化结果
all_genes_scores <- cor_results %>%
  mutate(
    # 计算绝对相关性作为评分
    score = abs(cor_coef),
    # 标记是否显著相关
    significant = ifelse(abs(cor_coef) > 0.3 & fdr < 0.05, "Significant", "Not Significant")
  ) %>%
  arrange(desc(score))  # 按评分降序排列

# 选择TOP20基因
top20_genes <- all_genes_scores %>%
  head(20) %>%
  mutate(
    # 标记是否在显著基因中（前16个）
    bar_color = ifelse(row_number() <= nrow(emt_related_genes), "Significant", "Top20_Not_Sig")
  )

# 创建颜色映射
color_mapping <- c(
  "Significant" = "#E41A1C",  # 红色：显著基因
  "Top20_Not_Sig" = "#B3B3B3"  # 灰色：非显著基因
)

# 绘制简洁的基因评分条形图（无标题，无图例）
p_barplot <- ggplot(top20_genes, aes(x = reorder(gene, score), y = score, fill = bar_color)) +
  # 绘制条形图
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  # 添加显著性阈值线
  geom_hline(yintercept = 0.3, color = "black", linetype = "dashed", linewidth = 0.8) +
  # 在条形上添加数值标签
  geom_text(aes(label = sprintf("%.2f", score)), 
            hjust = -0.1, size = 3.2, fontface = "bold") +
  # 设置颜色
  scale_fill_manual(values = color_mapping) +
  # 坐标轴翻转
  coord_flip() +
  # 设置Y轴范围
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, max(top20_genes$score) * 1.1)
  ) +
  # 坐标轴标签
  labs(
    x = "Gene",
    y = "|ρ|"
  ) +
  # 简洁的主题设置
  theme_classic(base_size = 12) +
  theme(
    # 去除标题和图例
    plot.title = element_blank(),
    legend.position = "none",
    # 坐标轴设置
    axis.text.y = element_text(face = "italic", size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    # 去除网格线
    panel.grid = element_blank(),
    # 调整边距
    plot.margin = margin(10, 15, 10, 10, "pt")
  )

# 显示图形
print(p_barplot)

# 保存为PDF格式
pdf("/home/datahup/xzp/paper/PCa/fig/AUcell/Top20_Genes_Barplot.pdf", 
    width = 8, height = 6)
print(p_barplot)
dev.off()