#Step2-1  CopyKAT
rm(list = ls())
gc()
getwd()

setwd('/home/datahup/xzp/paper/PCa/data/AUcell')
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(copykat)


# CopyKAT
epithelial_types <- c("Luminal Epithelial Cells", "Epithelial Cells")
load(file = '../step1/res0.3_PCa_harmony.Rdata')
levels(PCa_harmony) <- as.character(0:19)
new.cluster.ids <- c(
  "0"  = "T Cells",
  "1"  = "Smooth Muscle Cells",
  "2"  = "Luminal Epithelial Cells",
  "3"  = "Epithelial Cells",
  "4"  = "Macrophages",
  "5"  = "Epithelial Cells",
  "6"  = "Endothelial Cells",      
  "7"  = "Smooth Muscle Cells",
  "8"  = "T Cells",
  "9"  = "Mast Cells",             
  "10" = "B Cells",
  "11" = "T Cells",
  "12" = "Luminal Epithelial Cells",
  "13" = "Epithelial Cells",
  "14" = "Epithelial Cells",
  "15" = "Luminal Epithelial Cells",
  "16" = "Epithelial Cells",
  "17" = "Proliferating Cells",
  "18" = "Pericytes",
  "19" = "Pericytes"
)
PCa_harmony <- RenameIdents(PCa_harmony, new.cluster.ids)
PCa_harmony@meta.data[["CellType"]] <- Idents(PCa_harmony)

epithelial_T_cells <- rownames(PCa_harmony@meta.data)[PCa_harmony@meta.data$CellType %in% epithelial_types]
PCa_Epi <- subset(PCa_harmony, cells = epithelial_T_cells)
print(levels(PCa_Epi))
table(PCa_Epi@meta.data$CellType)
save(PCa_Epi,file = "/home/datahup/xzp/paper/PCa/data/AUcell/PCa_Epi.Rdata")

##运行时间较长  Time difference of 9.183234 days
load(file = "/home/datahup/xzp/paper/PCa/data/AUcell/PCa_Epi.Rdata")
counts <- as.matrix(PCa_Epi@assays$RNA@counts)
cnv <- copykat(rawmat=counts, ngene.chr=5, sam.name="PCa", n.cores=8)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
saveRDS(cnv, "copykat_PCaEpi.rds")

CNV <- readRDS(file = "copykat_PCaEpi.rds")

mallignant <- read.delim("PCa_copykat_prediction.txt")
mall <- mallignant
rownames(mall) <- mall[,1]
length(unique(mall$cell.names)) # 有重复
mall <- mall[!duplicated(mall$cell.names),] # 去除重复
rownames(mall) <- mall[,1]
names(mall)
b <- data.frame('copykat.pred' = mall[,-1])
rownames(b) <- rownames(mall)
table(mall$copykat.pred)

# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(PCa_Epi, metadata = b)
table(scRNA@meta.data[["copykat.pred"]])

library(ggplot2)
Epi <- scRNA
p1 <- DimPlot(Epi, group.by = "RNA_snn_res.0.2", label = T)
p2 <-DimPlot(Epi, group.by = "copykat.pred")
pc <- p1 + p2
ggsave("pred_mallignant.png", pc, width = 14, height = 8,path = "/home/datahup/xzp/paper/PCa/data/AUcell")

save(Epi,file = "copykat_Epi.Rdata")

table(Epi@meta.data[["copykat.pred"]])
CancerCell  = Epi[,Epi@meta.data$copykat.pred %in% "aneuploid"] 
Normal_Epi  = Epi[,Epi@meta.data$copykat.pred %in% "diploid"] 

save(CancerCell,file = "copykat_CancerCell.Rdata")
save(Normal_Epi,file = "copykat_Normal_Epi.Rdata")

# 加载数据
load("copykat_Epi.Rdata")
load("copykat_CancerCell.Rdata")
load("copykat_Normal_Epi.Rdata")

# 加载包
library(ggplot2)
library(Seurat)

# 设置输出目录
fig_dir <- "/home/datahup/xzp/paper/PCa/fig/CopyKAT"
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# 1. 对上皮细胞(Epi)进行tSNE降维
DefaultAssay(Epi) <- "RNA"

# 标准Seurat流程
Epi <- NormalizeData(Epi)
Epi <- FindVariableFeatures(Epi, nfeatures = 2000)
Epi <- ScaleData(Epi)
Epi <- RunPCA(Epi, npcs = 30)

# tSNE降维
Epi <- FindNeighbors(Epi, dims = 1:25)
Epi <- FindClusters(Epi, resolution = 0.5)
Epi <- RunTSNE(Epi, dims = 1:25, check_duplicates = FALSE)

# 保存tSNE降维后的上皮细胞PDF
p1 <- DimPlot(
  Epi,
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
  filename = file.path(fig_dir, "Epithelial_Cells_tSNE.pdf"),
  plot = p1,
  width = 8,
  height = 7,
  device = "pdf"
)

# 2. 参考原脚本简化恶性VS良性的标注
Epi$malignant.simple <- ifelse(
  Epi$copykat.pred == "diploid",
  "Epithelial Cells",
  "Cancer cell"
)

# 统计细胞类型
table(Epi$malignant.simple)

# 生成恶性vs良性细胞的tSNE图
p2 <- DimPlot(
  Epi,
  reduction = "tsne",
  group.by = "malignant.simple",
  label = TRUE,
  label.box = TRUE,
  pt.size = 0.8,
  repel = TRUE
) +
  scale_color_manual(
    values = c("Epithelial Cells" = "#1f78b4", "Cancer cell" = "#e31a1c"),
    name = "Cell Type"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "Malignant_vs_Normal_tSNE.pdf"),
  plot = p2,
  width = 9,
  height = 7,
  device = "pdf"
)