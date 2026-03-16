# Step 2-2: AUcell EMT Scoring for Cancer Cells
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/AUcell')
library(Seurat)
library(AUCell)
library(msigdbr)
library(dplyr)
library(ggplot2)

load("copykat_CancerCell.Rdata")
# 1. 从MSigDB获取EMT相关基因集
msigdb_h <- msigdbr(species = "Homo sapiens", category = "H")
emt_geneset <- msigdb_h %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)

print(paste("Number of genes in EMT gene set:", length(emt_geneset)))
head(emt_geneset)

# 2. 准备AUcell输入数据
exprMatrix <- as.matrix(GetAssayData(CancerCell, assay = "RNA", slot = "data"))
exprMatrix <- as.matrix(exprMatrix)   # 确保矩阵是数值矩阵

# 3. 计算AUCell评分
# 构建基因排名
cells_rankings <- AUCell_buildRankings(exprMatrix, 
                                       plotStats = TRUE) # 查看基因分布情况

emt_auc <- AUCell_calcAUC(geneSets = list(EMT = emt_geneset), 
                          rankings = cells_rankings,
                          aucMaxRank = ceiling(0.05 * nrow(cells_rankings))) # 默认参数

emt_scores <- getAUC(emt_auc)["EMT", ]
CancerCell <- AddMetaData(CancerCell, 
                          metadata = emt_scores, 
                          col.name = "EMT_AUCell_Score")

# 4. 可视化结果
library(ggplot2)
library(patchwork)

# 直方图
p_hist <- ggplot(data.frame(EMT_Score = CancerCell$EMT_AUCell_Score), 
                 aes(x = EMT_Score)) +
  geom_histogram(fill = "#2c7bb6", color = "white", alpha = 0.8, bins = 30) +
  labs(title = "EMT Score Distribution",
       x = "AUCell Score", 
       y = "Cell Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

print(p_hist)
ggsave("EMT_Score_Histogram.png", p_hist, width = 6, height = 4, dpi = 300,path = "/home/datahup/xzp/paper/PCa/fig/AUcell")

# UMAP图
p_umap <- FeaturePlot(CancerCell, features = "EMT_AUCell_Score", 
                      reduction = "umap",
                      pt.size = 0.3) +
  scale_colour_gradientn(colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59")) +
  labs(title = "EMT Score on UMAP",
       color = "EMT Score") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_umap)
ggsave("EMT_Score_UMAP.png", p_umap, width = 7, height = 5, dpi = 300,path = "/home/datahup/xzp/paper/PCa/fig/AUcell")

saveRDS(CancerCell, file = "EMTscore.rds")
write.csv(data.frame(Cell = colnames(CancerCell), 
                     EMT_Score = CancerCell$EMT_AUCell_Score),
          file = "AUCell_Scores.csv", row.names = FALSE)
