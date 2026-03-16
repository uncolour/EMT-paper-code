# Step 2-3: highEMT vs lowEMT DEGs
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/AUcell')
library(Seurat)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

CancerCell <- readRDS("EMTscore.rds")
# 1. 根据EMT评分对细胞进行分群
# 使用中位数作为分界点
emt_median <- median(CancerCell$EMT_AUCell_Score)
CancerCell$EMT_Group <- ifelse(CancerCell$EMT_AUCell_Score > emt_median, 
                               "High_EMT", "Low_EMT")
print(table(CancerCell$EMT_Group))

# 2. 可视化分组结果
# 提取上皮细胞
DefaultAssay(CancerCell) <- "RNA"

# 标准Seurat流程
CancerCell <- NormalizeData(CancerCell)
CancerCell <- FindVariableFeatures(CancerCell, nfeatures = 2000)
CancerCell <- ScaleData(CancerCell)
CancerCell <- RunPCA(CancerCell, npcs = 30)

# tSNE降维
CancerCell <- FindNeighbors(CancerCell, dims = 1:25)
CancerCell <- FindClusters(CancerCell, resolution = 0.5)
CancerCell <- RunTSNE(CancerCell, dims = 1:25, check_duplicates = FALSE)

# 保存包含EMT分组的数据
saveRDS(CancerCell, "CancerCell_EMTgroup.rds")

# 3. 生成tSNE图 - 按EMT分组
p1 <- DimPlot(
  CancerCell,
  reduction = "tsne",
  label = TRUE,
  label.size = 4,
  pt.size = 0.8,
  repel = TRUE
) +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(fig_dir, "CancerCell_tSNE_unannotated.pdf"),
  plot = p1,
  width = 8,
  height = 7,
  device = "pdf"
)

# 4. 生成第二张图：高EMT vs 低EMT注释图
p2 <- DimPlot(
  CancerCell,
  reduction = "tsne",
  group.by = "EMT_Group",
  cols = c("High_EMT" = "#e31a1c", "Low_EMT" = "#1f78b4"),
  label = TRUE,
  label.box = TRUE,
  pt.size = 0.8,
  repel = TRUE
) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = file.path(fig_dir, "CancerCell_EMT_Group_tSNE.pdf"),
  plot = p2,
  width = 9,
  height = 7,
  device = "pdf"
)

# 绘制小提琴图
p_vln_group <- VlnPlot(CancerCell, 
                       features = "EMT_AUCell_Score",
                       group.by = "EMT_Group",
                       pt.size = 0,
                       cols = c("#377eb8", "#e41a1c")) +
  ggtitle("EMT Score by Group") +
  theme_minimal()

print(p_vln_group)

# 保存为PDF格式
ggsave("/home/datahup/xzp/paper/PCa/fig/AUcell/EMT_Score_by_Group.pdf", p_vln_group, width = 6, height = 5, dpi = 300)

# 3. 寻找差异基因
Idents(CancerCell) <- CancerCell$EMT_Group

# 寻找高EMT vs 低EMT的差异基因
de_genes <- FindMarkers(CancerCell,
                        ident.1 = "High_EMT",
                        ident.2 = "Low_EMT",
                        test.use = "wilcox",  # 使用Wilcoxon检验
                        logfc.threshold = 0.25, # 降低logFC阈值以捕获更多基因
                        min.pct = 0.1,        # 在至少10%的细胞中表达
                        only.pos = FALSE)     # 包括上调和下调基因

cat("总差异基因数:", nrow(de_genes), "\n")
cat("上调基因数:", sum(de_genes$avg_log2FC > 0 & de_genes$p_val_adj < 0.05), "\n")
cat("下调基因数:", sum(de_genes$avg_log2FC < 0 & de_genes$p_val_adj < 0.05), "\n")

de_genes$gene <- rownames(de_genes)
de_genes$direction <- ifelse(de_genes$avg_log2FC > 0, "Up", "Down")
de_genes$significance <- ifelse(de_genes$p_val_adj < 0.05, "Significant", "Not significant")

# 4. 绘制火山图
library(ggplot2)
library(ggrepel)

for_volcano <- data.frame(
  'avg_log2FC' = de_genes$avg_log2FC,
  'p_val_adj' = de_genes$p_val_adj,
  'State' = 'Not significant',
  'gene' = de_genes$gene,
  row.names = de_genes$gene
)

# 定义阈值（使用您FindMarkers中设置的参数）
fc_threshold <- 0.5  # 与logfc.threshold一致
padj_threshold <- 0.01

# 标记显著差异基因
for_volcano$State[for_volcano$p_val_adj < padj_threshold & 
                    for_volcano$avg_log2FC > fc_threshold] <- 'Up'
for_volcano$State[for_volcano$p_val_adj < padj_threshold & 
                    for_volcano$avg_log2FC < -fc_threshold] <- 'Down'

# 转换p值为-log10
for_volcano$neg_log10_padj <- -log10(for_volcano$p_val_adj)

# 统计各类型基因数量
up_count <- sum(for_volcano$State == 'Up')
down_count <- sum(for_volcano$State == 'Down')
ns_count <- sum(for_volcano$State == 'Not significant')

# 创建标题
plot_title <- paste0('High_EMT vs Low_EMT Differential Expression\n',
                     'Up: ', up_count, ' | Down: ', down_count, ' | NS: ', ns_count)

# 选择要标记的top基因（上下调各10个最显著的）
#top_up_genes <- for_volcano %>%
#  filter(State == "Up") %>%
#  arrange(p_val_adj) %>%
#  head(10) %>%
#  pull(gene)

#top_down_genes <- for_volcano %>%
#  filter(State == "Down") %>%
#  arrange(p_val_adj) %>%
#  head(10) %>%
#  pull(gene)

#top_genes <- c(top_up_genes, top_down_genes)

# 保存火山图
library(scales)

# 自定义非线性变换：极大增强零附近区域的压缩，并继续压缩显著性线
p2 <- ggplot(for_volcano, aes(x = avg_log2FC, y = neg_log10_padj)) +
  # 散点图使用渐变色
  geom_point(aes(size = neg_log10_padj, color = neg_log10_padj), alpha = 0.7) +
  # 添加阈值线（显著性线）
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "#999999") +
  geom_hline(yintercept = -log10(padj_threshold), 
             linetype = "dashed", color = "#999999") +
  # 颜色渐变
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  # 大小渐变
  scale_size_continuous(range = c(1, 1.5)) +
  # 自定义x轴：极大增强零附近的压缩，压缩显著性线
  scale_x_continuous(
    limits = c(-2, 2),  # 设置x轴的范围
    breaks = NULL,       # 不显示刻度线
    labels = NULL,       # 不显示标签
    trans = scales::trans_new(
      "smooth_shift", 
      transform = function(x) {
        # 极大增强压缩，特别是零附近
        ifelse(abs(x) < 0.1, sign(x) * log1p(abs(x)) * 15, x)  # 强烈压缩零附近
      }, 
      inverse = function(x) {
        # 逆变换：恢复原始尺度
        sign(x) * (exp(abs(x)) - 1) / 15  # 恢复到正常尺度，除以15
      }
    )
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),  # 不显示标题
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_blank(),  # 去掉x轴的数值
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = 'log2(Fold Change)',
    y = '-log10(Adjusted p-value)',
    color = '-log10(padj)',
    size = '-log10(padj)'
  )

# 检查图形效果
print(p2)

# 保存为PDF
ggsave("/home/datahup/xzp/paper/PCa/fig/AUcell/EMT_DEG_Volcano_Enhanced.pdf", 
       p2, width = 8, height = 6, dpi = 300)
ggsave("/home/datahup/xzp/paper/PCa/fig/AUcell/EMT_DEG_Volcano_Basic.png", 
       p1, width = 8, height = 6, dpi = 300)
