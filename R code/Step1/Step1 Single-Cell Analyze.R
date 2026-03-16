#Step1-1  Single-Cell Analyze
rm(list = ls())
gc()
getwd()

library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

# GSE141445
setwd("/home/datahup/xzp/paper/PCa/data/raw_data_scRNA/GSE141445/")

PCa13_counts <- read.table("data.raw.matrix.txt", 
                           header = TRUE, 
                           row.names = 1, 
                           sep = "\t", 
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

cat("GSE141445数据维度:", dim(PCa13_counts), "\n")
cat("细胞数量:", ncol(PCa13_counts), "\n")
cat("基因数量:", nrow(PCa13_counts), "\n")

PCa13 <- CreateSeuratObject(
  counts = as.matrix(PCa13_counts),
  project = "GSE141445",
  min.cells = 3,
  min.features = 200
)

# 添加样本信息
PCa13$sample <- "GSE141445"
PCa13$dataset <- "GSE141445"
PCa13$tissue <- "Tumor"

# 添加线粒体基因比例
PCa13[["percent.mt"]] <- PercentageFeatureSet(PCa13, pattern = "^MT-")

# GSE181294
setwd("/home/datahup/xzp/paper/PCa/data/raw_data_scRNA/GSE181294/")

# 处理GSM61392X批次样本
sample_prefixes <- c(
  "GSM6133921_S5",
  "GSM6133922_S6", 
  "GSM6133923_S7",
  "GSM6133924_S8",
  "GSM6133925_S9",
  "GSM6133926_S10",
  "GSM6133927_S11",
  "GSM6133928_S12"
)

seurat_list <- list()

for (prefix in sample_prefixes) {
  cat("Processing:", prefix, "...\n")
  
  mtx_path <- paste0(prefix, ".counts.mtx.gz")
  genes_path <- paste0(prefix, ".genes.csv.gz")
  barcodes_path <- paste0(prefix, ".barcode.csv.gz")
  
  # 方法1：使用ReadMtx并跳过标题行
  counts <- ReadMtx(
    mtx = mtx_path,
    features = genes_path,
    cells = barcodes_path,
    feature.column = 1,
    skip.cell = 1,  # 跳过barcodes的标题行
    skip.feature = 1  # 跳过genes的标题行
  )
  
  # 创建Seurat对象
  seu <- CreateSeuratObject(
    counts = counts, 
    project = prefix, 
    min.cells = 3, 
    min.features = 200
  )
  
  # 重命名细胞ID
  seu <- RenameCells(seu, add.cell.id = prefix)
  
  # 添加样本信息
  seu$sample <- prefix
  seu$group <- ifelse(grepl("S[5-8]", prefix), "Benign", "Tumor")
  
  # 添加线粒体基因比例
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  seurat_list[[prefix]] <- seu
  
  cat("Completed:", prefix, "- Cells:", ncol(seu), "Genes:", nrow(seu), "\n\n")
}

# 合并所有样本
PCa8 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
#save(PCa8, file = "PCa8.Rdata")

# 处理GSM549437X批次的样本
files <- list.files(pattern = "\\.count\\.csv\\.gz$", full.names = TRUE)
tumor_files <- files[grepl("-T-", files)]

cat("找到", length(tumor_files), "个癌症样本:\n")
print(basename(tumor_files))

count_list <- list()

for (file in tumor_files) {
  # 提取样本名
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- sub("\\.count$", "", sample_name)
  
  cat("Processing:", sample_name, "...\n")
  
  # 读取count数据
  df <- read.csv(gzfile(file), row.names = 1, check.names = FALSE)
  
  # 重命名列名以包含样本信息
  colnames(df) <- paste(sample_name, colnames(df), sep = "_")
  
  # 存储到列表
  count_list[[sample_name]] <- df
  
  cat("Completed:", sample_name, "- Cells:", ncol(df), "Genes:", nrow(df), "\n\n")
}

# 合并所有基因名
all_genes <- unique(unlist(lapply(count_list, rownames)))
cat("总基因数:", length(all_genes), "\n")

# 填充缺失的基因（用0填充）
count_list_filled <- lapply(count_list, function(mat) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if(length(missing_genes) > 0) {
    add_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                      dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, add_mat)
  }
  mat[all_genes, , drop = FALSE]
})

# 合并所有样本的表达矩阵
merged_counts <- do.call(cbind, count_list_filled)

# 创建Seurat对象
PCa_T <- CreateSeuratObject(counts = merged_counts, project = "PCa_T")

# 添加样本信息到元数据
sample_names <- sapply(strsplit(colnames(PCa_T), "_"), function(x) paste(x[1:3], collapse = "_"))
PCa_T$sample <- sample_names

# 添加组织类型信息
PCa_T$tissue <- ifelse(grepl("-T-", PCa_T$sample), "Tumor", 
                       ifelse(grepl("-N-", PCa_T$sample), "Normal", "Unknown"))

# 添加分级信息
PCa_T$grade <- ifelse(grepl("-HG", PCa_T$sample), "HighGrade",
                      ifelse(grepl("-LG", PCa_T$sample), "LowGrade", "Unknown"))

# 添加线粒体基因比例
PCa_T[["percent.mt"]] <- PercentageFeatureSet(PCa_T, pattern = "^MT-")
#save(PCa_T, file = "PCa18.Rdata")


# GSE193337
setwd("/home/datahup/xzp/paper/PCa/data/raw_data_scRNA/GSE193337/")

samples <- data.frame(
  gsm_id = c("GSM5793828", "GSM5793829", "GSM5793831", "GSM5793832"),
  sample_code = c("P1t", "P2t", "P3t", "P4t"),
  patient = c("Patient_1", "Patient_2", "Patient_3", "Patient_4")
)

seurat_list <- list()

for(i in 1:nrow(samples)) {
  gsm_id <- samples$gsm_id[i]
  sample_code <- samples$sample_code[i]
  
  cat("Processing:", gsm_id, "(", sample_code, ")...\n")
  
  # 构建文件路径
  barcodes_path <- paste0(gsm_id, "_", sample_code, "_barcodes.tsv.gz")
  features_path <- paste0(gsm_id, "_", sample_code, "_features.tsv.gz")
  matrix_path <- paste0(gsm_id, "_", sample_code, "_matrix.mtx.gz")
  
  # 检查文件是否存在
  if(!all(file.exists(c(barcodes_path, features_path, matrix_path)))) {
    cat("Missing files for", gsm_id, "\n")
    next
  }
  
  # 方法1：手动读取并创建矩阵（推荐）
  cat("Reading files manually...\n")
  
  # 读取barcodes（细胞）
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)[, 1]
  
  # 读取features（基因）- 取第一列ENSEMBL ID
  features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)[, 1]
  
  # 读取矩阵
  counts_matrix <- readMM(matrix_path)
  
  # 设置行列名 - 维度已经匹配
  rownames(counts_matrix) <- features
  colnames(counts_matrix) <- barcodes
  
  # 创建Seurat对象
  seu <- CreateSeuratObject(
    counts = counts_matrix, 
    project = gsm_id, 
    min.cells = 3, 
    min.features = 200
  )
  
  # 重命名细胞ID
  seu <- RenameCells(seu, add.cell.id = sample_code)
  
  # 添加样本信息
  seu$sample <- gsm_id
  seu$sample_code <- sample_code
  seu$patient <- samples$patient[i]
  seu$tissue <- "Tumor"
  seu$dataset <- "GSE193337"
  
  # 添加线粒体基因比例（使用ENSEMBL ID模式）
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  seurat_list[[gsm_id]] <- seu
  
  cat("Completed:", gsm_id, "- Cells:", ncol(seu), "Genes:", nrow(seu), "\n\n")
}

# 检查是否有成功读取的样本
if(length(seurat_list) == 0) {
  stop("No samples were successfully loaded.")
}

# 合并所有样本
if(length(seurat_list) == 1) {
  PCa4_T <- seurat_list[[1]]
} else {
  PCa4_T <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
}
save(PCa4_T, file = "PCa4.Rdata")

#合并数据集
PCa_list <- list(PCa1 = PCa13, PCa2 = PCa8, PCa3 = PCa4_T,  PCa4 = PCa_T)
setwd('/home/datahup/xzp/paper/PCa/data/step1')
PCa.Integrate <- merge(
  x = PCa_list[[1]],
  y = PCa_list[2:4]
)
#save(PCa.Integrate, file = "PCa.Integrate.Rdata")
rm(list = ls())
gc()


# Step2 QC
setwd("/home/datahup/xzp/paper/PCa/data/step1")
load('/home/datahup/xzp/paper/PCa/data/step1/PCa.Integrate.Rdata')
library(Seurat)
library(ggplot2)
library(patchwork)

Idents(PCa.Integrate) = PCa.Integrate@meta.data[["orig.ident"]]
PCa.Integrate[["percent.mt"]] <- PercentageFeatureSet(PCa.Integrate, assay = "RNA",pattern = "^MT-")
head(PCa.Integrate@meta.data)
VlnPlot(PCa.Integrate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "QCbefore.png",width = 30, height = 12, path = "../../fig/Step1-1/")

plot1 <- FeatureScatter(PCa.Integrate, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(PCa.Integrate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
combined_plot_before <- plot1 + plot2
print(combined_plot_before)
ggsave(filename = "QCbefore2.png", width = 30, height = 12, path = "../../fig/Step1-1/", plot = combined_plot_before)
ggsave(filename = "QCbefore2.pdf", width = 30, height = 12, path = "../../fig/Step1-1/", plot = combined_plot_before)

#根据上图进行质控
PCa.Integrate <- subset(PCa.Integrate, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15
                          & nCount_RNA > 200 & nCount_RNA < 40000)
PCa.Integrate
VlnPlot(PCa.Integrate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "QCafter.png",width = 30, height = 12, path = "../../fig/Step1-1/")

plot1 <- FeatureScatter(PCa.Integrate, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(PCa.Integrate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
combined_plot_after <- plot1 + plot2
print(combined_plot_after)
ggsave(filename = "QCafter2.png", width = 30, height = 12, path = "../../fig/Step1-1/", plot = combined_plot_after)
ggsave(filename = "QCafter2.pdf", width = 30, height = 12, path = "../../fig/Step1-1/", plot = combined_plot_after)

# Step 3 标准化数据+找高变基因
DefaultAssay(PCa.Integrate) <- "RNA"
PCa.Integrate <- NormalizeData(PCa.Integrate)
PCa.Integrate <- FindVariableFeatures(PCa.Integrate, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(PCa.Integrate), 10)
plot1 <- VariableFeaturePlot(PCa.Integrate)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(filename = "Top10-VarGene.png", plot = plot1 + plot2,width = 20,height = 10,path = "../../fig/Step1-1/")
save(PCa.Integrate,file="PCa.Integrate.Rdata")

# Step 4 数据归一化
all.genes <- rownames(PCa.Integrate)
PCa.Integrate <- ScaleData(PCa.Integrate, features = all.genes)

# Step 5 线性降维 
PCa.Integrate <- RunPCA(PCa.Integrate, features = VariableFeatures(object = PCa.Integrate))
save(PCa.Integrate,file="PCa.Integrate2.Rdata")
load(file = 'PCa.Integrate2.Rdata')
# 查看PCA结果
print(PCa.Integrate[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(PCa.Integrate, dims = 1:2, reduction = "pca")
ggsave(filename = "PCA.png",width = 16,height = 10,path = "../../fig/Step1-1/")

DimPlot(PCa.Integrate, reduction = "pca", raster=FALSE)
ggsave(filename = "PCA2.png",width = 16,height = 10,path = "../../fig/Step1-1/") 

png("../../fig/Step1-1/PC1_HeatmapPlot.png", width = 1600, height = 1000, res = 150)
DimHeatmap(PCa.Integrate, dims = 1, cells = 500, balanced = TRUE)
dev.off()

png("../../fig/Step1-1/PC15_HeatmapPlot.png", width = 1600, height = 1000, res = 150)
DimHeatmap(PCa.Integrate, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

PCa.Integrate <<- JackStraw(PCa.Integrate, num.replicate = 100) # 02h 30m 29s
PCa.Integrate <<- ScoreJackStraw(PCa.Integrate, dims = 1:20)
save(PCa.Integrate,file="PCa.Integrate3.Rdata")
JackStrawPlot(PCa.Integrate, dims = 1:20)
ggsave(filename = "JackPlot.png",width = 16,height = 10,path = "../../fig/Step1-1/")
ElbowPlot(PCa.Integrate)#肘部图 
# 综合以上方法，选择10个主成成分作为参数用于后续分析。
ggsave(filename = "ElbowPlot.png",width = 12,height = 10,path = "../../fig/Step1-1/")

# Step 6 Harmony整合
library(harmony)
load("PCa.Integrate2.Rdata")
PCa.harmony <- RunHarmony(PCa.Integrate, group.by.vars = "orig.ident")
save(PCa.harmony,file="PCa.harmony.Rdata")
#降维聚类
library(Seurat)
PCa.harmony = RunUMAP(PCa.harmony, reduction = "harmony", dims = 1:15)
PCa.harmony = FindNeighbors(PCa.harmony, reduction = "harmony", dims = 1:15) 
PCa.harmony = FindClusters(object = PCa.harmony, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
library(clustree)
clustree(PCa.harmony@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "PCa.harmony_resolution(0-1).pdf",width = 20,height = 14,path = "../../fig/Step1-1/")
save(PCa.harmony,file="PCa.harmony2.Rdata")

Idents(object = PCa.harmony) <- "RNA_snn_res.0.3"
PCa.harmony@meta.data$seurat_clusters = PCa.harmony@meta.data$RNA_snn_res.0.3
head(Idents(PCa.harmony), 5)#查看前5个细胞的分类ID

plot1 = DimPlot(PCa.harmony, reduction = "umap", label=T) 
plot1
ggsave(filename = "PCa.harmony_res0.3.pdf", plot = plot1, width = 8,height = 6,path = "../../fig/Step1-1/")

#调色板
library(ggsci)
pal = pal_ucscgb(alpha = 0.8)(22)
plot2 = DimPlot(PCa.harmony, reduction = "umap", label=F, group.by = 'orig.ident') 
plot2
# ggsave(filename = "scRNA_harmony_patient.pdf", plot = plot2, width = 10,height = 6,path = "./Fig")
ggsave(filename = "PCa.harmony_res0.3_orig.ident.pdf", plot = plot2, width = 10,height = 6,path = "../../fig/Step1-1/")
save(PCa.harmony,file = 'PCa.harmony3.Rdata')

# Step 7  找每个簇的差异基因
load(file = 'PCa.harmony3.Rdata')
library(dplyr)
PCa.markers <- FindAllMarkers(PCa.harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-PCa.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(PCa.harmony, features = top10$gene) + NoLegend()
ggsave(filename = "Top10-MarkerGene.png",width = 32,height = 20,path ="../../fig/Step1-1/")
save(top10,PCa.markers,file = "res0.3_Markers.Rdata")
save(PCa.harmony, file = "res0.3_PCa_harmony.Rdata")

rm(list = ls())
gc()

# Step 8 自动化注释  SingleR
load(file = "res0.3_Markers.Rdata")
load(file = "res0.3_PCa_harmony.Rdata")
sce = PCa.harmony
library(celldex)
library(SingleR)
library(Seurat)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.3
clusters=sce@meta.data$seurat_clusters

load(file="/home/datahup/tyh/SCRNA/coda/singleRref.Rdata")

Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )

head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']

pro='SingleR_anno_res0.3'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../fig/Step1-1/")

PCa_N.harmony = sce

save(PCa_N.harmony, file = "res0.3_PCa_harmony.Rdata")

write.csv(cellType,file = "cellType_res0.3.csv")
write.csv(PCa.markers, file =  "PCa_markers_res0.3.csv")
write.csv(top10, file = "top10gene_res0.3.csv")
rm(list = ls())
gc()


##Step 9 手动注释  SingleR
# C0,C8,C11 : T cells                              Marker:"CD3D","IL7R","CCL5"
# C1,C7 : Smooth Muscle Cells                      Marker:"ACTG2","CSRP1"
# C2,C12,C15 : Luminal Epithelial Cells            Marker:"AZGP1","MSMB","FXYD3"
# C3,C5,C13,C14,C16 : Epithelial Cells             Marker:"PPAP2A","MYO6","FOLH1"
# C4 : Macrophages                                 Marker:"C1QA","LYZ"
# C6 : Endothelial Cells                           Marker:"VWF","IFI27"
# C9 : Mast Cells                                  Marker:"CPA3","TPSAB1"
# C10 : B Cells                                    Marker:"MS4A1","CD79A"
# C17 : Proliferating Cells                        Marker:"TOP2A","CD79A"
# C18,C19 : Pericytes                              Marker:"C11orf96","RERGL"

# 各簇top基因气泡图
library(dplyr)
load(file = "res0.3_PCa_harmony.Rdata")
load(file = "res0.3_Markers.Rdata")
# cluster 0
genes_to_check=subset(top10,top10$cluster==0)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C0.png",path = "../../fig/Step1-1/Cellmaker/cluster0")
#Marker:"IL7R"
VlnPlot(PCa.harmony, features = "IL7R")
ggsave(filename = "IL7R.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster0")
FeaturePlot(PCa.harmony, features = "IL7R")
ggsave(filename = "IL7R(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster0")

# cluster 1
genes_to_check=subset(top10,top10$cluster==1)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C1.png",path = "../../fig/Step1-1/Cellmaker/cluster1")
#Marker:"ACTG2"
VlnPlot(PCa.harmony, features = "ACTG2")
ggsave(filename = "ACTG2.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster1")
FeaturePlot(PCa.harmony, features = "ACTG2")
ggsave(filename = "ACTG2(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster1")

# cluster 2
genes_to_check=subset(top10,top10$cluster==2)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C2.png",path = "../../fig/Step1-1/Cellmaker/cluster1")
#Marker:"MSMB"
VlnPlot(PCa.harmony, features = "MSMB")
ggsave(filename = "MSMB.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster2")
FeaturePlot(PCa.harmony, features = "MSMB")
ggsave(filename = "MSMB(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster2")

# cluster 3
genes_to_check=subset(top10,top10$cluster==3)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C3.png",path = "../../fig/Step1-1/Cellmaker/cluster3")
#Marker:"PPAP2A"
VlnPlot(PCa.harmony, features = "PPAP2A")
ggsave(filename = "PPAP2A.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster3")
FeaturePlot(PCa.harmony, features = "PPAP2A")
ggsave(filename = "PPAP2A(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster3")

# cluster 4
genes_to_check=subset(top10,top10$cluster==4)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C4.png",path = "../../fig/Step1-1/Cellmaker/cluster4")
#Marker:"C1QA"
VlnPlot(PCa.harmony, features = "C1QA")
ggsave(filename = "C1QA.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster4")
FeaturePlot(PCa.harmony, features = "C1QA")
ggsave(filename = "C1QA(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster4")

# cluster 5
genes_to_check=subset(top10,top10$cluster==5)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C5.png",path = "../../fig/Step1-1/Cellmaker/cluster5")

# cluster 6
genes_to_check=subset(top10,top10$cluster==6)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C6.png",path = "../../fig/Step1-1/Cellmaker/cluster6")
#Marker:"VWF"
VlnPlot(PCa.harmony, features = "VWF")
ggsave(filename = "VWF.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster6")
FeaturePlot(PCa.harmony, features = "VWF")
ggsave(filename = "VWF(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster6")

# cluster 7
genes_to_check=subset(top10,top10$cluster==7)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C7.png",path = "../../fig/Step1-1/Cellmaker/cluster7")

# cluster 8
genes_to_check=subset(top10,top10$cluster==8)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C8.png",path = "../../fig/Step1-1/Cellmaker/cluster8")

# cluster 9
genes_to_check=subset(top10,top10$cluster==9)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C9.png",path = "../../fig/Step1-1/Cellmaker/cluster9")
#Marker:"CPA3"
VlnPlot(PCa.harmony, features = "CPA3")
ggsave(filename = "CPA3.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster9")
FeaturePlot(PCa.harmony, features = "CPA3")
ggsave(filename = "CPA3(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster9")

# cluster 10
genes_to_check=subset(top10,top10$cluster==10)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C10.png",path = "../../fig/Step1-1/Cellmaker/cluster10")
#Marker:"MS4A1"
VlnPlot(PCa.harmony, features = "MS4A1")
ggsave(filename = "MS4A1.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster10")
FeaturePlot(PCa.harmony, features = "MS4A1")
ggsave(filename = "MS4A1(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster10")

# cluster 11
genes_to_check=subset(top10,top10$cluster==11)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C11.png",path = "../../fig/Step1-1/Cellmaker/cluster11")
#Marker:"TRAC"
VlnPlot(PCa.harmony, features = "TRAC")
ggsave(filename = "TRAC.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster11")
FeaturePlot(PCa.harmony, features = "TRAC")
ggsave(filename = "TRAC(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster11")

# cluster 12
genes_to_check=subset(top10,top10$cluster==12)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C12.png",path = "../../fig/Step1-1/Cellmaker/cluster12")
#Marker:"LINC01088"
VlnPlot(PCa.harmony, features = "LINC01088")
ggsave(filename = "LINC01088.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster12")
FeaturePlot(PCa.harmony, features = "LINC01088")
ggsave(filename = "LINC01088(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster12")

# cluster 13
genes_to_check=subset(top10,top10$cluster==13)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C13.png",path = "../../fig/Step1-1/Cellmaker/cluster13")
#Marker:"C1orf64"
VlnPlot(PCa.harmony, features = "C1orf64")
ggsave(filename = "C1orf64.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster13")
FeaturePlot(PCa.harmony, features = "C1orf64")
ggsave(filename = "C1orf64(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster13")

# cluster 14
genes_to_check=subset(top10,top10$cluster==14)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C14.png",path = "../../fig/Step1-1/Cellmaker/cluster14")

# cluster 15
genes_to_check=subset(top10,top10$cluster==15)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C15.png",path = "../../fig/Step1-1/Cellmaker/cluster15")
#Marker:"FABP5"
VlnPlot(PCa.harmony, features = "FABP5")
ggsave(filename = "FABP5.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster15")
FeaturePlot(PCa.harmony, features = "FABP5")
ggsave(filename = "FABP5(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster15")

# cluster 15
genes_to_check=subset(top10,top10$cluster==15)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C15.png",path = "../../fig/Step1-1/Cellmaker/cluster15")
#Marker:"FABP5"
VlnPlot(PCa.harmony, features = "FABP5")
ggsave(filename = "FABP5.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster15")
FeaturePlot(PCa.harmony, features = "FABP5")
ggsave(filename = "FABP5(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster15")

# cluster 16
genes_to_check=subset(top10,top10$cluster==16)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C16.png",path = "../../fig/Step1-1/Cellmaker/cluster16")
#Marker:"ERG"
VlnPlot(PCa.harmony, features = "ERG")
ggsave(filename = "ERG.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster16")
FeaturePlot(PCa.harmony, features = "ERG")
ggsave(filename = "ERG(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster16")

# cluster 17
genes_to_check=subset(top10,top10$cluster==17)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C17.png",path = "../../fig/Step1-1/Cellmaker/cluster17")
#Marker:"TOP2A"
VlnPlot(PCa.harmony, features = "TOP2A")
ggsave(filename = "TOP2A.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster17")
FeaturePlot(PCa.harmony, features = "TOP2A")
ggsave(filename = "TOP2A(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster17")

# cluster 18
genes_to_check=subset(top10,top10$cluster==18)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C18.png",path = "../../fig/Step1-1/Cellmaker/cluster18")
#Marker:"C11orf96"
VlnPlot(PCa.harmony, features = "C11orf96")
ggsave(filename = "C11orf96.png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster18")
FeaturePlot(PCa.harmony, features = "C11orf96")
ggsave(filename = "C11orf96(1).png",width = 18,height = 12,path = "../../fig/Step1-1/Cellmaker/cluster18")

# cluster 19
genes_to_check=subset(top10,top10$cluster==19)
DotPlot(PCa.harmony,features = unique(genes_to_check$gene)) + RotatedAxis()
ggsave(filename = "C19.png",path = "../../fig/Step1-1/Cellmaker/cluster19")

# 可视化
library(Seurat)
library(plyr)
library(ggplot2)

# 设置工作目录
setwd("/home/datahup/xzp/paper/PCa/data/step1")

# 加载数据
load("res0.3_PCa_harmony.Rdata")

# 根据注释表格定义细胞类型映射
new.cluster.ids <- c(
  "0" = "T Cells",                     # C0 : T cells
  "1" = "Smooth Muscle Cells",         # C1 : Smooth Muscle Cells
  "2" = "Luminal Epithelial Cells",    # C2 : Luminal Epithelial Cells
  "3" = "Epithelial Cells",            # C3 : Epithelial Cells
  "4" = "Macrophages",                 # C4 : Macrophages
  "5" = "Epithelial Cells",            # C5 : Epithelial Cells
  "6" = "Endothelial Cells",           # C6 : Endothelial Cells
  "7" = "Smooth Muscle Cells",         # C7 : Smooth Muscle Cells
  "8" = "T Cells",                     # C8 : T cells
  "9" = "Mast Cells",                  # C9 : Mast Cells
  "10" = "B Cells",                    # C10 : B Cells
  "11" = "T Cells",                    # C11 : T cells
  "12" = "Luminal Epithelial Cells",   # C12 : Luminal Epithelial Cells
  "13" = "Epithelial Cells",           # C13 : Epithelial Cells
  "14" = "Epithelial Cells",           # C14 : Epithelial Cells
  "15" = "Luminal Epithelial Cells",   # C15 : Luminal Epithelial Cells
  "16" = "Epithelial Cells",           # C16 : Epithelial Cells
  "17" = "Proliferating Cells",        # C17 : Proliferating Cells
  "18" = "Pericytes",                  # C18 : Pericytes
  "19" = "Pericytes"                   # C19 : Pericytes
)

# 应用注释
Idents(PCa.harmony) <- PCa.harmony@meta.data$RNA_snn_res.0.3
PCa_N.harmony <- RenameIdents(PCa.harmony, new.cluster.ids)

# 给meta.data增加合并后的celltype列
PCa_N.harmony@meta.data$celltype_merged <- plyr::mapvalues(
  as.character(PCa_N.harmony@meta.data$RNA_snn_res.0.3),
  from = names(new.cluster.ids),
  to = new.cluster.ids
)

# 定义细胞类型顺序
desired_order <- c(
  "T Cells",
  "Smooth Muscle Cells",
  "Luminal Epithelial Cells",
  "Epithelial Cells",
  "Macrophages",
  "Endothelial Cells",
  "Mast Cells",
  "B Cells",
  "Proliferating Cells",
  "Pericytes"
)

PCa_N.harmony@meta.data$celltype_merged <- factor(
  PCa_N.harmony@meta.data$celltype_merged,
  levels = desired_order
)

# 根据注释表格定义标志基因
marker_list <- list(
  "T Cells" = c("CD3D", "IL7R", "CCL5"),
  "Smooth Muscle Cells" = c("ACTG2", "CSRP1"),
  "Luminal Epithelial Cells" = c("AZGP1", "MSMB", "FXYD3"),
  "Epithelial Cells" = c("PPAP2A", "MYO6", "FOLH1"),
  "Macrophages" = c("C1QA", "LYZ"),
  "Endothelial Cells" = c("VWF", "IFI27"),
  "Mast Cells" = c("CPA3", "TPSAB1"),
  "B Cells" = c("MS4A1", "CD79A"),
  "Proliferating Cells" = c("TOP2A", "CD79A"),
  "Pericytes" = c("C11orf96", "RERGL")
)

all_markers <- unique(unlist(marker_list))

# 绘制整体标志基因点图
dp <- DotPlot(
  PCa_N.harmony,
  features = unlist(marker_list),
  group.by = "celltype_merged"
) +
  scale_color_gradientn(colors = c("#F0F0F0", "#A1D99B", "#006D2C")) +
  RotatedAxis() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    strip.text.x = element_blank()
  ) +
  guides(size = guide_legend(override.aes = list(color = "grey40")))

# 保存为PDF
ggsave(
  filename = "DotPlot_Cellmarkers.pdf",
  plot = dp,
  width = 14,
  height = 8,
  path = "../../fig/Step1-2/"
)

# 保存为PNG
ggsave(
  filename = "DotPlot_Cellmarkers.png",
  plot = dp,
  width = 14,
  height = 8,
  path = "../../fig/Step1-2/"
)

# 保存注释后的数据
save(PCa_N.harmony, file = "res0.3_PCa_N.harmony1.Rdata")