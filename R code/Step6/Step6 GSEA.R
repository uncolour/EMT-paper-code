### Step6 Patients GSEA GSVA
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/GSEA')

# 创建输出目录
fig_dir <- '/home/datahup/xzp/paper/PCa/fig/GSEA'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat("结果将保存到:", fig_dir, "\n")

# 加载包
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(GSVA)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(KEGGREST)
library(ggplot2)
library(dplyr)

# 1. 加载建模基因和系数
selected_genes_df <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/core_emt_genes.csv")
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes
print(selected_genes_df)

# 加载三组数据的完整表达矩阵
# TCGA数据 - count矩阵
tcga_counts <- readRDS("../AUcell-Model/TCGA_PRAD_counts_matrix.rds")
cat("TCGA Count矩阵维度:", dim(tcga_counts), "\n")

# GSE62116数据 - 芯片数据
gse62116_expr <- read.csv("../raw_data_GEO/GSE62116/GSE62116_expression_matrix_annotated.csv", row.names = 1)
cat("GSE62116 表达矩阵维度:", dim(gse62116_expr), "\n")

# GSE220095数据 - count矩阵
load("../raw_data_GEO/GSE220095/count_matrix.RData")
cat("GSE220095 Count矩阵维度:", dim(count_matrix), "\n")

# 数据标准化
# 简单的CPM标准化函数
normalize_counts <- function(count_matrix) {
  # 转置为基因×样本（如果需要）
  if(nrow(count_matrix) < ncol(count_matrix)) {
    count_matrix <- t(count_matrix)
  }
  
  # 过滤低表达基因
  keep_genes <- rowSums(count_matrix > 0) >= 0.1 * ncol(count_matrix)
  count_matrix <- count_matrix[keep_genes, ]
  
  # CPM标准化
  lib_sizes <- colSums(count_matrix)
  normalized <- t(t(count_matrix) / lib_sizes) * 1e6
  log_normalized <- log2(normalized + 1)
  
  return(log_normalized)
}

# 标准化各数据集
tcga_normalized <- normalize_counts(tcga_counts)
gse220095_normalized <- normalize_counts(count_matrix)

# GSE62116已经是标准化数据，直接使用
if(nrow(gse62116_expr) < ncol(gse62116_expr)) {
  gse62116_normalized <- t(gse62116_expr)
} else {
  gse62116_normalized <- as.matrix(gse62116_expr)
}

cat("标准化后数据维度:\n")
cat("TCGA:", dim(tcga_normalized), "\n")
cat("GSE62116:", dim(gse62116_normalized), "\n")
cat("GSE220095:", dim(gse220095_normalized), "\n")

# 寻找共同基因
common_genes <- Reduce(intersect, list(
  rownames(tcga_normalized),
  rownames(gse62116_normalized), 
  rownames(gse220095_normalized)
))

cat("共同基因数:", length(common_genes), "\n")

# 统一计算EMT评分
calculate_emt_score <- function(expr_matrix, genes, coefficients) {
  # 找到可用的基因
  available_genes <- intersect(rownames(expr_matrix), genes)
  cat("  可用基因数:", length(available_genes), "/", length(genes), "\n")
  
  if(length(available_genes) == 0) {
    stop("没有可用的建模基因")
  }
  
  # 提取基因表达数据
  gene_expr <- expr_matrix[available_genes, , drop = FALSE]
  
  # 获取对应的系数
  available_coef <- coefficients[available_genes]
  
  # 计算EMT评分
  emt_scores <- as.numeric(t(gene_expr) %*% available_coef)
  
  # 创建结果数据框
  result <- data.frame(
    sample_id = colnames(expr_matrix),
    emt_score = emt_scores,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# 计算各数据集的EMT评分
tcga_emt <- calculate_emt_score(tcga_normalized, selected_genes, selected_coef)
gse62116_emt <- calculate_emt_score(gse62116_normalized, selected_genes, selected_coef)
gse220095_emt <- calculate_emt_score(gse220095_normalized, selected_genes, selected_coef)

cat("EMT评分计算完成:\n")
cat("TCGA:", nrow(tcga_emt), "个样本\n")
cat("GSE62116:", nrow(gse62116_emt), "个样本\n")
cat("GSE220095:", nrow(gse220095_emt), "个样本\n")

# 整合数据和统一分组
# 提取共同基因的表达数据
tcga_common <- tcga_normalized[common_genes, ]
gse62116_common <- gse62116_normalized[common_genes, ]
gse220095_common <- gse220095_normalized[common_genes, ]

# 合并表达矩阵
combined_expr <- cbind(tcga_common, gse62116_common, gse220095_common)
cat("整合表达矩阵维度:", dim(combined_expr), "\n")

# 合并EMT评分
combined_emt <- bind_rows(
  tcga_emt %>% mutate(dataset = "TCGA"),
  gse62116_emt %>% mutate(dataset = "GSE62116"),
  gse220095_emt %>% mutate(dataset = "GSE220095")
)

# 确保样本匹配
common_samples <- intersect(colnames(combined_expr), combined_emt$sample_id)
combined_expr_final <- combined_expr[, common_samples]
combined_emt_final <- combined_emt %>% filter(sample_id %in% common_samples)

cat("最终匹配样本数:", length(common_samples), "\n")

# 统一按所有样本的EMT评分中位数分组
emt_median <- median(combined_emt_final$emt_score, na.rm = TRUE)
combined_emt_final$emt_group <- ifelse(combined_emt_final$emt_score > emt_median, "High_EMT", "Low_EMT")
combined_emt_final$emt_group <- factor(combined_emt_final$emt_group, levels = c("Low_EMT", "High_EMT"))

cat("EMT分组分布:\n")
print(table(combined_emt_final$emt_group))
cat("数据集分布:\n")
print(table(combined_emt_final$dataset, combined_emt_final$emt_group))

# 保存预处理数据
save(combined_expr_final, file = "combined_expression_unified.RData")
save(combined_emt_final, file = "combined_emt_unified.RData")
save(common_genes, file = "common_genes_unified.RData")

# 2. 差异表达分析
# 确保样本顺序一致
expr_samples <- colnames(combined_expr_final)
emt_samples <- combined_emt_final$sample_id

cat("表达矩阵样本数:", length(expr_samples), "\n")
cat("EMT信息样本数:", length(emt_samples), "\n")

# 检查样本是否完全匹配
if(!identical(sort(expr_samples), sort(emt_samples))) {
  cat("样本不匹配，进行修正...\n")
  common_samples <- intersect(expr_samples, emt_samples)
  combined_expr_final <- combined_expr_final[, common_samples]
  combined_emt_final <- combined_emt_final %>% filter(sample_id %in% common_samples)
  cat("修正后样本数:", length(common_samples), "\n")
}

# 确保样本顺序一致
expr_matrix <- combined_expr_final
emt_info <- combined_emt_final %>% arrange(sample_id)
expr_matrix <- expr_matrix[, emt_info$sample_id]  # 按EMT信息的顺序排列表达矩阵

cat("最终确认:\n")
cat("表达矩阵维度:", dim(expr_matrix), "\n")
cat("EMT信息样本数:", nrow(emt_info), "\n")
cat("样本顺序一致:", identical(colnames(expr_matrix), emt_info$sample_id), "\n")

# 现在进行差异表达分析
design <- model.matrix(~ emt_info$emt_group)
cat("设计矩阵维度:", dim(design), "\n")
cat("设计矩阵行数应等于表达矩阵列数:", nrow(design) == ncol(expr_matrix), "\n")

# 执行limma分析
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
deg_results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
deg_results$gene_symbol <- rownames(deg_results)

cat("差异表达分析完成:\n")
cat("显著差异基因 (FDR < 0.05):", sum(deg_results$adj.P.Val < 0.05), "\n")
cat("显著上调基因:", sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC > 0), "\n")
cat("显著下调基因:", sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC < 0), "\n")

# 火山图可视化
for_volcano <- data.frame(
  'avg_log2FC' = deg_results$logFC,
  'p_val_adj' = deg_results$adj.P.Val,
  'State' = 'Not significant',
  'gene' = deg_results$gene_symbol,
  row.names = deg_results$gene_symbol,
  stringsAsFactors = FALSE
)

# 定义阈值
fc_threshold <- 0.5  # log2 fold change阈值
padj_threshold <- 0.05  # FDR阈值

# 标记显著差异基因
for_volcano$State[for_volcano$p_val_adj < padj_threshold & 
                    for_volcano$avg_log2FC > fc_threshold] <- 'Up'
for_volcano$State[for_volcano$p_val_adj < padj_threshold & 
                    for_volcano$avg_log2FC < -fc_threshold] <- 'Down'

# 转换p值为-log10
for_volcano$neg_log10_padj <- -log10(for_volcano$p_val_adj)

# 统计各类型基因数量
up_count <- sum(for_volcano$State == 'Up', na.rm = TRUE)
down_count <- sum(for_volcano$State == 'Down', na.rm = TRUE)
ns_count <- sum(for_volcano$State == 'Not significant', na.rm = TRUE)

# 创建标题
plot_title <- paste0('Differential Expression: High_EMT vs Low_EMT\n',
                     'Up: ', up_count, ' | Down: ', down_count, ' | NS: ', ns_count)

# 选择要标记的top基因（上下调各10个最显著的）
top_up_genes <- for_volcano %>%
  filter(State == "Up") %>%
  arrange(p_val_adj) %>%
  head(10) %>%
  pull(gene)

top_down_genes <- for_volcano %>%
  filter(State == "Down") %>%
  arrange(p_val_adj) %>%
  head(10) %>%
  pull(gene)

top_genes <- c(top_up_genes, top_down_genes)

# 绘制火山图
p1 <- ggplot(for_volcano, aes(x = avg_log2FC, y = neg_log10_padj)) +
  # 散点图，按状态着色
  geom_point(aes(color = State), alpha = 0.6, size = 2) +
  # 添加阈值线
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(padj_threshold), 
             linetype = "dashed", color = "gray50", linewidth = 0.6) +
  # 颜色方案
  scale_color_manual(
    name = "Expression",
    values = c("Up" = "#B2182B", "Down" = "#2166AC", "Not significant" = "gray70"),
    breaks = c("Up", "Down", "Not significant"),
    labels = c(paste0("Up (", up_count, ")"), 
               paste0("Down (", down_count, ")"), 
               paste0("NS (", ns_count, ")"))
  ) +
  # 主题设置（简洁版）
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.title = element_blank(),  # 移除标题
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15),
    # 移除标题后的额外调整
    plot.title.position = "plot"
  ) +
  # 标签（只保留坐标轴标签）
  labs(
    x = 'log₂(Fold Change)',
    y = '-log₁₀(Adjusted p-value)',
    title = NULL  # 确保标题为空
  ) +
  # 坐标轴范围
  xlim(c(-max(abs(for_volcano$avg_log2FC), na.rm = TRUE) * 1.1, 
         max(abs(for_volcano$avg_log2FC), na.rm = TRUE) * 1.1)) +
  ylim(c(0, max(for_volcano$neg_log10_padj, na.rm = TRUE) * 1.05)) +
  # 可选：在右上角添加统计信息
  annotate("text", 
           x = max(for_volcano$avg_log2FC, na.rm = TRUE) * 0.9,
           y = max(for_volcano$neg_log10_padj, na.rm = TRUE) * 0.95,
           label = paste0("Up: ", up_count, "\nDown: ", down_count),
           hjust = 1, vjust = 1,
           size = 3.5,
           color = "black",
           fontface = "bold")

# 保存基础火山图
volcano_basic_file <- file.path(fig_dir, "EMT_DEG_Volcano_Basic.pdf")
ggsave(volcano_basic_file, p1, width = 8, height = 6, dpi = 300)
cat(sprintf("基础火山图已保存: %s\n", basename(volcano_basic_file)))

# 另外创建一个完全纯净的版本（无任何统计信息）
p1 <- ggplot(for_volcano, aes(x = avg_log2FC, y = neg_log10_padj)) +
  # 散点图，按状态着色
  geom_point(aes(color = State), alpha = 0.6, size = 2) +
  # 添加阈值线
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(padj_threshold), 
             linetype = "dashed", color = "gray50", linewidth = 0.6) +
  # 颜色方案
  scale_color_manual(
    name = "Expression",
    values = c("Up" = "#B2182B", "Down" = "#2166AC", "Not significant" = "gray70"),
    breaks = c("Up", "Down", "Not significant"),
    labels = c(paste0("Up (", up_count, ")"), 
               paste0("Down (", down_count, ")"), 
               paste0("NS (", ns_count, ")"))
  ) +
  # 主题设置（简洁版）
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.title = element_blank(),  # 移除标题
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15),
    # 移除标题后的额外调整
    plot.title.position = "plot"
  ) +
  # 标签（简化坐标轴标签，不用log下标）
  labs(
    x = 'Log2 Fold Change',
    y = '-Log10 Adjusted p-value',
    title = NULL  # 确保标题为空
  ) +
  # 坐标轴范围
  xlim(c(-max(abs(for_volcano$avg_log2FC), na.rm = TRUE) * 1.1, 
         max(abs(for_volcano$avg_log2FC), na.rm = TRUE) * 1.1)) +
  ylim(c(0, max(for_volcano$neg_log10_padj, na.rm = TRUE) * 1.05)) +
  # 在右上角添加统计信息
  annotate("text", 
           x = max(for_volcano$avg_log2FC, na.rm = TRUE) * 0.9,
           y = max(for_volcano$neg_log10_padj, na.rm = TRUE) * 0.95,
           label = paste0("Up: ", up_count, "\nDown: ", down_count),
           hjust = 1, vjust = 1,
           size = 3.5,
           color = "black",
           fontface = "bold")

volcano_basic_file <- file.path(fig_dir, "EMT_DEG_Volcano_Basic.pdf")
ggsave(volcano_basic_file, p1, width = 8, height = 6, dpi = 300)

# 3. 准备GSEA分析
gene_map <- bitr(deg_results$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
cat("成功映射的基因数:", nrow(gene_map), "\n")

# 准备GSEA输入
deg_mapped <- deg_results %>%
  inner_join(gene_map, by = c("gene_symbol" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID) & !is.na(logFC))

gene_list <- deg_mapped$logFC
names(gene_list) <- deg_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("GSEA输入基因数:", length(gene_list), "\n")
cat("基因列表统计:\n")
cat("  - 范围:", range(gene_list), "\n")
cat("  - 中位数:", median(gene_list), "\n")
cat("  - 正logFC基因:", sum(gene_list > 0), "\n")
cat("  - 负logFC基因:", sum(gene_list < 0), "\n")

# 获取基因集
get_kegg_genesets_full <- function() {
  kegg_genesets <- list()
  
  # 获取所有人类KEGG通路
  pathway_list <- keggList("pathway", "hsa")
  pathway_ids <- names(pathway_list)
  
  cat("所有KEGG通路数:", length(pathway_ids), "\n")
  
  # 获取所有通路的基因集（参考脚本方法）
  for (pathway_id in pathway_ids) {
    tryCatch({
      pathway_info <- keggGet(pathway_id)
      
      if (!is.null(pathway_info[[1]]$GENE)) {
        genes <- unlist(lapply(pathway_info[[1]]$GENE, function(x) strsplit(x, ';')))
        genelist <- genes[1:length(genes) %% 3 == 2]  # 提取基因符号
        
        # 不限制基因集大小（参考脚本方法）
        pathway_name <- gsub(" - Homo sapiens \\(human\\)", "", pathway_info[[1]]$NAME)
        kegg_genesets[[pathway_name]] <- genelist
      }
    }, error = function(e) {
      # 忽略错误，继续下一个通路
    })
    
    # 进度显示
    if (length(kegg_genesets) %% 50 == 0) {
      cat("已处理", length(kegg_genesets), "个KEGG通路...\n")
    }
  }
  
  return(kegg_genesets)
}

kegg_genesets <- get_kegg_genesets_full()
cat("成功获取KEGG基因集数:", length(kegg_genesets), "\n")

# 获取完整GO基因集（BP+CC+MF，不限制数量）
get_go_genesets_full <- function() {
  # 使用参考脚本的方法：从org.Hs.eg.db获取完整GO数据
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # 方法1: 使用参考脚本中的getGO函数（如果可用）
  if(file.exists("getGoTerm.R")) {
    source("getGoTerm.R")
    GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
    
    # 转换为基因集列表
    go_genesets <- list()
    for (go_id in names(GO_DATA$PATHID2EXTID)) {
      genes <- GO_DATA$PATHID2EXTID[[go_id]]
      go_name <- GO_DATA$PATHID2NAME[[go_id]]
      go_genesets[[go_name]] <- genes
    }
    
  } else {
    # 方法2: 使用备用方法获取完整GO
    cat("使用备用方法获取GO基因集...\n")
    
    # 获取所有GO BP
    go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
    # 获取所有GO CC
    go_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
    # 获取所有GO MF  
    go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
    
    # 合并所有GO，不限制数量
    all_go <- rbind(go_bp, go_cc, go_mf)
    
    # 转换为列表格式
    go_genesets <- split(all_go$gene_symbol, all_go$gs_name)
  }
  
  return(go_genesets)
}

go_genesets <- get_go_genesets_full()
cat("成功获取GO基因集数:", length(go_genesets), "\n")

# 合并基因集并转换为GSEA格式
all_genesets_combined <- c(
  kegg_genesets,
  go_genesets
)
cat("总基因集数:", length(all_genesets_combined), "\n")

# 统计基因集大小分布
geneset_sizes <- sapply(all_genesets_combined, length)
cat("基因集大小统计:\n")
cat("  最小:", min(geneset_sizes), "\n")
cat("  最大:", max(geneset_sizes), "\n")
cat("  平均:", round(mean(geneset_sizes), 1), "\n")
cat("  中位数:", median(geneset_sizes), "\n")

# 转换为GSEA需要的TERM2GENE格式
correct_columns <- data.frame(
  gs_name = character(),
  gene_symbol = character(),  # 使用gene_symbol而不是entrez_gene
  stringsAsFactors = FALSE
)

for (i in 1:length(all_genesets_combined)) {
  gs_name <- names(all_genesets_combined)[i]
  genes <- all_genesets_combined[[gs_name]]
  
  # 直接使用基因符号，不需要转换
  temp_df <- data.frame(
    gs_name = gs_name,
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )
  correct_columns <- rbind(correct_columns, temp_df)
  
  # 进度显示
  if (i %% 100 == 0) {
    cat("已处理", i, "/", length(all_genesets_combined), "个基因集...\n")
  }
}

# 移除可能的NA值
correct_columns <- correct_columns %>% filter(!is.na(gene_symbol))
cat("清理后基因集数据维度:", dim(correct_columns), "\n")
cat("唯一基因集数:", length(unique(correct_columns$gs_name)), "\n")
cat("总基因-通路关联数:", nrow(correct_columns), "\n")

# 准备GSEA输入
gene_list_symbol <- deg_results$logFC
names(gene_list_symbol) <- deg_results$gene_symbol
gene_list_symbol <- sort(gene_list_symbol, decreasing = TRUE)

cat("GSEA输入基因数（符号）:", length(gene_list_symbol), "\n")
cat("基因列表统计:\n")
cat("  - 范围:", range(gene_list_symbol), "\n")
cat("  - 中位数:", median(gene_list_symbol), "\n")
cat("  - 正logFC基因:", sum(gene_list_symbol > 0), "\n")
cat("  - 负logFC基因:", sum(gene_list_symbol < 0), "\n")

# 执行GSEA分析（使用基因符号）
gsea_result <- GSEA(
  geneList = gene_list_symbol,  # 使用基因符号的基因列表
  exponent = 1,
  minGSSize = 10,       # 提高最小基因集大小，过滤太小基因集
  maxGSSize = 500,      # 限制最大基因集大小，过滤过大基因集
  eps = 1e-10,
  pvalueCutoff = 0.01,  # 降低p值阈值，更严格
  pAdjustMethod = "BH",
  TERM2GENE = correct_columns,  # 使用基因符号的基因集
  seed = 123,
  by = "fgsea"
)

cat("总分析通路数:", nrow(gsea_result@result), "\n")
cat("显著通路数 (p.adjust < 0.01):", sum(gsea_result@result$p.adjust < 0.01), "\n")

# 进一步筛选：结合NES值和FDR
sig_gsea_df <- gsea_result@result %>%
  filter(p.adjust < 0.01 & abs(NES) > 1.5)  # 增加NES绝对值阈值

cat("严格筛选后显著通路数 (FDR<0.01 & |NES|>1.5):", nrow(sig_gsea_df), "\n")

# 按NES绝对值排序
sig_gsea_df <- sig_gsea_df %>%
  arrange(desc(abs(NES)))

# 结果保存
write.csv(sig_gsea_df, file = "GSEA_strict_significant_results.csv", row.names = FALSE)
write.csv(gsea_result@result, file = "GSEA_complete_results.csv", row.names = FALSE)

geneset_info <- list(
  kegg_genesets = kegg_genesets,
  go_genesets = go_genesets,
  all_genesets_combined = all_genesets_combined,
  correct_columns = correct_columns
)
saveRDS(geneset_info, file = "GSEA_genesets_data.rds")

# 全GSEA可视化
supp_fig1 <- ggplot(sig_gsea_df %>%
                      arrange(NES) %>%
                      mutate(
                        short_name = sapply(ID, function(x) {
                          if(nchar(x) > 40) paste0(substr(x, 1, 37), "...") else x
                        }),
                        short_name = factor(short_name, levels = short_name),
                        direction = ifelse(NES > 0, "High-EMT enriched", "Low-EMT enriched"),
                        log_fdr = -log10(p.adjust)
                      ),
                    aes(x = NES, 
                        y = short_name,
                        size = sapply(core_enrichment, function(x) {
                          if(is.na(x) || x == "") return(0)
                          length(strsplit(x, "/")[[1]])
                        }),
                        color = log_fdr)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(
    name = "Core genes",
    range = c(2, 6),
    breaks = c(20, 40, 60)
  ) +
  scale_color_viridis_c(
    name = "-log10(FDR)",
    option = "plasma",
    direction = -1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 5)),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing = unit(0.2, "cm"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig_dir, "supplementary_figure_gsea_bubble.pdf"),
       supp_fig1, width = 10, height = 8, dpi = 300)

# 4.GSEA结果可视化
selected_paths <- c(
  "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT",
  "GOMF_IMMUNOGLOBULIN_BINDING",
  "Neuroactive ligand-receptor interaction",
  "GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY",
  "GOBP_SENSORY_PERCEPTION_OF_CHEMICAL_STIMULUS",
  "GOMF_HORMONE_ACTIVITY",
  "GOMF_CHEMOKINE_ACTIVITY",
  "GOMF_NEUROPEPTIDE_RECEPTOR_ACTIVITY",
  "GOBP_SENSORY_PERCEPTION_OF_SMELL",
  "GOMF_CYTOKINE_ACTIVITY",
  "GOCC_COSTAMERE",
  "GOCC_RNAI_EFFECTOR_COMPLEX",
  "GOBP_GABAERGIC_NEURON_DIFFERENTIATION"
)

negative_pathways <- c(
  "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT",
  "GOCC_COSTAMERE", 
  "GOCC_RNAI_EFFECTOR_COMPLEX"
)

positive_pathways <- c(
  "GOMF_IMMUNOGLOBULIN_BINDING",
  "Neuroactive ligand-receptor interaction",
  "GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY",
  "GOBP_SENSORY_PERCEPTION_OF_CHEMICAL_STIMULUS",
  "GOMF_HORMONE_ACTIVITY",
  "GOMF_CHEMOKINE_ACTIVITY",
  "GOMF_NEUROPEPTIDE_RECEPTOR_ACTIVITY",
  "GOBP_SENSORY_PERCEPTION_OF_SMELL",
  "GOMF_CYTOKINE_ACTIVITY",
  "GOBP_GABAERGIC_NEURON_DIFFERENTIATION"
)

# 创建GSEA可视化目录
gsea_viz_dir <- file.path(fig_dir, "GSEA_plots")
dir.create(gsea_viz_dir, recursive = TRUE, showWarnings = FALSE)
cat("GSEA可视化结果将保存到:", gsea_viz_dir, "\n")

# 获取GSEA结果中的所有ID
all_gsea_ids <- gsea_result@result$ID
cat("GSEA结果中共有", length(all_gsea_ids), "个通路ID\n")

# 检查每个通路
found_pathways <- character()
missing_pathways <- character()

for (pathway in selected_paths) {
  if (pathway %in% all_gsea_ids) {
    found_pathways <- c(found_pathways, pathway)
    cat(sprintf("✓ %s\n", pathway))
  } else {
    missing_pathways <- c(missing_pathways, pathway)
    cat(sprintf("✗ %s\n", pathway))
  }
}

cat(sprintf("\n匹配结果: %d个成功, %d个失败\n", 
            length(found_pathways), length(missing_pathways)))

# 绘图函数
draw_pathway_simple <- function(pathway_id, index) {
  cat(sprintf("[%d/%d] %s\n", index, length(selected_paths), pathway_id))
  
  # 获取通路信息
  pathway_info <- gsea_result@result[gsea_result@result$ID == pathway_id, ]
  
  # 提取信息
  nes_value <- pathway_info$NES
  p_value <- pathway_info$pvalue
  fdr_value <- pathway_info$p.adjust
  description <- pathway_info$Description
  
  cat(sprintf("  NES: %.3f\n", nes_value))
  
  # 确定颜色和方向
  if (nes_value > 0) {
    plot_color <- "#D62728"  # 红色
    direction <- "Positive (Up in High-EMT)"
  } else {
    plot_color <- "#1F77B4"  # 蓝色
    direction <- "Negative (Down in High-EMT)"
  }
  
  # 创建安全的文件名
  safe_name <- gsub("[^[:alnum:]_]", "_", pathway_id)
  
  # 使用绝对路径
  output_file <- file.path(fig_dir, paste0(safe_name, ".png"))
  output_file_abs <- normalizePath(output_file, mustWork = FALSE)
  
  cat(sprintf("  保存到: %s\n", output_file_abs))
  
  # 简化p值和FDR显示
  # p值显示为 <0.001 或具体值
  if (p_value < 0.001) {
    p_display <- "p < 0.001"
  } else {
    p_display <- sprintf("p = %.3f", p_value)
  }
  
  # FDR显示为 <0.01 或具体值
  if (fdr_value < 0.01) {
    fdr_display <- "FDR < 0.01"
  } else if (fdr_value < 0.05) {
    fdr_display <- "FDR < 0.05"
  } else {
    fdr_display <- sprintf("FDR = %.3f", fdr_value)
  }
  
  # 构建简洁的标题
  plot_title <- paste0(
    description, 
    "\nNES = ", round(nes_value, 2), " (", direction, ")",
    "\n", p_display, ", ", fdr_display
  )
  
  # 绘图
  png(filename = output_file_abs, 
      width = 1400, 
      height = 900, 
      res = 150)
  
  # 绘制GSEA图
  result_plot <- gseaplot2(
    gsea_result,
    geneSetID = pathway_id,
    title = plot_title,
    color = plot_color,
    base_size = 12
  )
  
  # 打印图形
  print(result_plot)
  
  # 关闭设备
  dev.off()
  
  cat("  ✓ 绘图完成\n")
  return(TRUE)
}

plot_list <- list()
for (pathway in selected_paths) {
  if (pathway %in% gsea_result@result$ID) {
    # 生成并保存单个GSEA图
    p <- gseaplot2(gsea_result, geneSetID = pathway, title = pathway)
    plot_list[[pathway]] <- p
  }
}

# 绘制所有通路
success_count <- 0
for (i in seq_along(selected_paths)) {
  pathway_id <- selected_paths[i]
  success <- draw_pathway_simple(pathway_id, i)
  if (success) success_count <- success_count + 1
  cat("\n")
}
cat(sprintf("\n绘制完成: %d/%d 个通路成功\n", 
            success_count, length(selected_paths)))

# 输出拼接PDF
pdf_file_clean <- file.path(fig_dir, "GSEA_compact_clean.pdf")
pdf(pdf_file_clean, width = 12, height = 8)

# 设置紧凑参数
par(mfrow = c(5, 3),
    mar = rep(0, 4),      # 无边距
    oma = rep(0, 4),      # 无外边界
    xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n",
    bty = "n")

# 重新初始化计数器
plot_counter <- 1
negative_plot_counter <- 1

for (row in 1:5) {
  for (col in 1:3) {
    # 清空绘图区域
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)
    
    # 确定绘制哪个图片
    if (col %in% 1:2) {
      # 前两列：正富集
      if (plot_counter <= 10) {
        # 绘制图片（无边框，无标签）
        if (!is.null(plot_list[[plot_counter]])) {
          rasterImage(plot_list[[plot_counter]], 0, 0, 1, 1)
        }
        plot_counter <- plot_counter + 1
      }
    } else if (col == 3 && row <= 3) {
      # 第三列：负富集（只在前3行）
      if (negative_plot_counter <= 3) {
        plot_index <- 10 + negative_plot_counter
        
        if (!is.null(plot_list[[plot_index]])) {
          rasterImage(plot_list[[plot_index]], 0, 0, 1, 1)
        }
        negative_plot_counter <- negative_plot_counter + 1
      }
    }
    # 其他情况（第三列后2行）保持空白
  }
}

dev.off()

# 创建PNG版本
png_file_clean <- file.path(fig_dir, "GSEA_compact_clean.png")
cat(sprintf("创建纯图片PNG: %s\n", png_file_clean))

png(png_file_clean, width = 2400, height = 1600, res = 300)

# 同样参数
par(mfrow = c(5, 3),
    mar = rep(0, 4),
    oma = rep(0, 4),
    xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n",
    bty = "n")

plot_counter <- 1
negative_plot_counter <- 1

for (row in 1:5) {
  for (col in 1:3) {
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)
    
    if (col %in% 1:2) {
      if (plot_counter <= 10) {
        if (!is.null(plot_list[[plot_counter]])) {
          rasterImage(plot_list[[plot_counter]], 0, 0, 1, 1)
        }
        plot_counter <- plot_counter + 1
      }
    } else if (col == 3 && row <= 3) {
      if (negative_plot_counter <= 3) {
        plot_index <- 10 + negative_plot_counter
        if (!is.null(plot_list[[plot_index]])) {
          rasterImage(plot_list[[plot_index]], 0, 0, 1, 1)
        }
        negative_plot_counter <- negative_plot_counter + 1
      }
    }
  }
}
dev.off()

# 5.GSVA数据准备
sig_gsea_df <- read.csv("GSEA_strict_significant_results.csv")
cat(sprintf("GSEA显著通路数: %d\n", nrow(sig_gsea_df)))

# 取前38个通路
top_gsva_pathways <- sig_gsea_df %>%
  arrange(desc(abs(NES))) %>%
  head(38) %>%
  pull(ID)

# 从correct_columns中提取这些通路的基因集
gsva_genesets <- list()
for (pathway_id in top_gsva_pathways) {
  if (pathway_id %in% correct_columns$gs_name) {
    pathway_genes <- correct_columns %>%
      filter(gs_name == pathway_id) %>%
      pull(gene_symbol) %>%
      unique()
    
    if (length(pathway_genes) >= 10) {
      # 获取通路描述
      pathway_info <- gsea_result@result[gsea_result@result$ID == pathway_id, ]
      pathway_name <- ifelse(nrow(pathway_info) > 0, 
                             pathway_info$Description, 
                             pathway_id)
      
      gsva_genesets[[pathway_name]] <- pathway_genes
    }
  }
}
cat(sprintf("准备 %d 个基因集用于GSVA分析\n", length(gsva_genesets)))

gsva_expr_matrix <- as.matrix(expr_matrix)
cat("GSVA表达矩阵维度:", dim(gsva_expr_matrix), "\n")

# 移除重复样本
gsva_expr_matrix <- as.matrix(expr_matrix)
cat(sprintf("表达矩阵维度: %d 基因 × %d 样本\n", 
            nrow(gsva_expr_matrix), ncol(gsva_expr_matrix)))

duplicate_samples <- emt_info$sample_id[duplicated(emt_info$sample_id)]
cat(sprintf("发现 %d 个重复样本，已自动处理\n", length(duplicate_samples)))

emt_info_clean <- emt_info[!duplicated(emt_info$sample_id), ]

# 创建样本注释
sample_annotation <- data.frame(
  row.names = emt_info_clean$sample_id,
  emt_group = emt_info_clean$emt_group,
  stringsAsFactors = FALSE
)

cat(sprintf("样本注释: %d 个样本\n", nrow(sample_annotation)))

# 6.执行分析
gsva_result <- gsva(
  expr = gsva_expr_matrix,
  gset.idx.list = gsva_genesets,
  method = "gsva",
  kcdf = "Gaussian",
  abs.ranking = FALSE,
  min.sz = 5,
  max.sz = 500,
  parallel.sz = 1,
  mx.diff = TRUE,
  verbose = TRUE
)
cat(sprintf("GSVA分析完成！结果维度: %d 通路 × %d 样本\n", 
            nrow(gsva_result), ncol(gsva_result)))

# 保存结果
saveRDS(gsva_result, file = "GSVA_scores_matrix.rds")
write.csv(as.data.frame(t(gsva_result)), file = "GSVA_scores.csv", row.names = TRUE)

# 匹配样本
gsva_samples <- colnames(gsva_result)
annotation_samples <- rownames(sample_annotation)
common_samples <- intersect(gsva_samples, annotation_samples)

cat(sprintf("匹配样本数: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
  stop("错误：GSVA结果和样本注释没有共同样本！")
}

# 过滤数据
gsva_result_filtered <- gsva_result[, common_samples, drop = FALSE]
sample_annotation_filtered <- sample_annotation[common_samples, , drop = FALSE]

# 确保顺序一致
gsva_result_filtered <- gsva_result_filtered[, rownames(sample_annotation_filtered), drop = FALSE]

cat(sprintf("过滤后数据: %d 通路 × %d 样本\n", 
            nrow(gsva_result_filtered), ncol(gsva_result_filtered)))

# 7.可视化

# 获取前38个显著通路
top_38_pathways <- sig_gsea_df %>%
  arrange(desc(abs(NES))) %>%
  head(38) %>%
  pull(ID)

# 获取通路信息
pathway_info_38 <- data.frame()
for (pathway_id in top_38_pathways) {
  if (pathway_id %in% gsea_result@result$ID) {
    info <- gsea_result@result[gsea_result@result$ID == pathway_id, ]
    pathway_info_38 <- rbind(pathway_info_38, data.frame(
      ID = pathway_id,
      Description = info$Description,
      NES = info$NES,
      pvalue = info$pvalue,
      p.adjust = info$p.adjust,
      stringsAsFactors = FALSE
    ))
  }
}

# 按NES值排序（从负到正）
pathway_info_38 <- pathway_info_38 %>%
  arrange(NES)

# 简化通路名称函数
simplify_pathway_name <- function(full_name) {
  # 移除常见前缀
  simplified <- gsub("^GOBP_", "", full_name)
  simplified <- gsub("^GOMF_", "", simplified)
  simplified <- gsub("^GOCC_", "", simplified)
  
  # 如果是长描述，提取关键词
  if (nchar(simplified) > 50) {
    # 尝试提取主要部分
    words <- unlist(strsplit(simplified, " "))
    if (length(words) > 5) {
      simplified <- paste(words[1:4], collapse = " ")
      simplified <- paste0(simplified, "...")
    } else {
      simplified <- substr(simplified, 1, 47)
      if (nchar(simplified) == 47) {
        simplified <- paste0(simplified, "...")
      }
    }
  }
  
  return(simplified)
}

# 检查通路匹配
available_pathways <- intersect(pathway_info_38$Description, rownames(gsva_result_filtered))
cat(sprintf("在GSVA结果中找到 %d/%d 个通路\n", 
            length(available_pathways), nrow(pathway_info_38)))

# 提取可用的通路数据
gsva_38_pathways <- pathway_info_38 %>%
  filter(Description %in% rownames(gsva_result_filtered)) %>%
  arrange(NES)

# 添加简化名称
gsva_38_pathways$Simple_Name <- sapply(gsva_38_pathways$Description, simplify_pathway_name)

# 获取GSVA评分
gsva_38_scores <- gsva_result_filtered[gsva_38_pathways$Description, ]

# 使用简化名称作为行名
rownames(gsva_38_scores) <- gsva_38_pathways$Simple_Name

cat(sprintf("用于可视化的数据:\n"))
cat(sprintf("  通路数: %d\n", nrow(gsva_38_scores)))
cat(sprintf("  样本数: %d\n", ncol(gsva_38_scores)))

# 排序通路和样本
pathway_order <- gsva_38_pathways$Simple_Name

# 样本排序：先按EMT分组，再在组内随机排序
sample_order_df <- data.frame(
  Sample = colnames(gsva_38_scores),
  EMT_Group = sample_annotation_filtered[colnames(gsva_38_scores), "emt_group"]
) %>%
  group_by(EMT_Group) %>%
  arrange(EMT_Group) %>%
  ungroup()

sample_order <- sample_order_df$Sample

# 排序矩阵
gsva_38_scores_ordered <- gsva_38_scores[pathway_order, sample_order]
sample_annotation_ordered <- sample_annotation_filtered[sample_order, , drop = FALSE]

# 创建热图
gsva_heatmap_file <- file.path(fig_dir, "GSVA_38_pathways_heatmap_simple.pdf")

pdf(gsva_heatmap_file, width = 12, height = 10)

pheatmap(
  gsva_38_scores_ordered,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  # 保留annotation_col以显示上方彩色方块
  annotation_col = annotation_for_display,
  annotation_colors = list(
    emt_group = c("Low_EMT" = "#1F77B4", "High_EMT" = "#D62728")
  ),
  fontsize_row = 7,
  gaps_row = sum(gsva_38_pathways$NES < 0),
  gaps_col = sum(sample_annotation_ordered$emt_group == "Low_EMT"),
  border_color = NA,
  cellheight = 8,
  treeheight_row = 0,
  treeheight_col = 0,
  main = NULL,
  # 关键设置：
  annotation_names_col = FALSE,  # 不显示"emt_group"文字标注
  annotation_legend = FALSE,     # ❌ 不显示样本分组图例
  legend = TRUE,                 # ✅ 只显示颜色图例
  fontsize = 8.5
)

dev.off()

# 创建通路活性差异条形图
# 计算每个通路在高/低EMT组的t检验差异
# 计算每个通路在高/低EMT组的t检验差异
pathway_differences <- data.frame()

for (i in 1:nrow(gsva_38_pathways)) {
  pathway_desc <- gsva_38_pathways$Description[i]
  pathway_simple <- gsva_38_pathways$Simple_Name[i]
  
  # 使用简化名称作为索引
  pathway_scores <- gsva_38_scores[pathway_simple, ]
  
  # 提取两组数据
  high_emt_scores <- pathway_scores[sample_annotation_filtered$emt_group == "High_EMT"]
  low_emt_scores <- pathway_scores[sample_annotation_filtered$emt_group == "Low_EMT"]
  
  # 计算t检验
  if (length(high_emt_scores) > 1 && length(low_emt_scores) > 1) {
    t_test_result <- t.test(high_emt_scores, low_emt_scores)
    
    # 计算均值差异
    mean_diff <- mean(high_emt_scores) - mean(low_emt_scores)
    
    pathway_differences <- rbind(pathway_differences, data.frame(
      Pathway = pathway_desc,
      Simple_Name = pathway_simple,
      Description = ifelse(nchar(pathway_simple) > 80, paste0(substr(pathway_simple, 1, 77), "..."), pathway_simple),
      NES = gsva_38_pathways$NES[i],
      t_value = t_test_result$statistic,
      p_value = t_test_result$p.value,
      mean_diff = mean_diff,
      high_emt_mean = mean(high_emt_scores),
      low_emt_mean = mean(low_emt_scores),
      stringsAsFactors = FALSE
    ))
  }
}

# 添加FDR校正
pathway_differences$fdr <- p.adjust(pathway_differences$p_value, method = "BH")

# 添加显著性标记（与原脚本完全一致）
pathway_differences$significance <- case_when(
  pathway_differences$fdr < 0.001 ~ "***",
  pathway_differences$fdr < 0.01 ~ "**",
  pathway_differences$fdr < 0.05 ~ "*",
  TRUE ~ "ns"
)

# 确定颜色：正t值（高风险组活性更高）用蓝色，负t值（低风险组活性更高）用绿色
pathway_differences$direction <- ifelse(pathway_differences$t_value > 0, 
                                        "Higher in High-EMT", 
                                        "Higher in Low-EMT")

# 按t值排序（从负到正）
pathway_differences <- pathway_differences %>%
  arrange(t_value)

# 确保因子顺序
pathway_differences$Description <- factor(
  pathway_differences$Description,
  levels = pathway_differences$Description
)

# 绘制差异条形图
# 获取统计信息
positive_t <- pathway_differences$t_value[pathway_differences$t_value > 0]
negative_t <- pathway_differences$t_value[pathway_differences$t_value < 0]
positive_t_max <- max(positive_t, na.rm = TRUE)
negative_t_min <- min(negative_t, na.rm = TRUE)

# 计算t绝对值并排序
pathway_differences <- pathway_differences %>%
  mutate(abs_t = abs(t_value)) %>%
  arrange(abs_t)

# 标记t绝对值最小的20个通路（与原脚本逻辑一致）
pathway_differences$is_weakest <- FALSE
pathway_differences$is_weakest[1:min(20, nrow(pathway_differences))] <- TRUE

# 检查标记是否正确
cat("Debug: Checking weakest pathway assignment...\n")
cat(sprintf("Total pathways: %d\n", nrow(pathway_differences)))
cat(sprintf("Weakest pathways marked: %d\n", sum(pathway_differences$is_weakest)))
cat(sprintf("Range of |t| for weakest: %.3f to %.3f\n", 
            min(pathway_differences$abs_t[pathway_differences$is_weakest]),
            max(pathway_differences$abs_t[pathway_differences$is_weakest])))

# 重新按t值从负到正排序（用于绘图）
pathway_differences <- pathway_differences %>%
  arrange(t_value)

# 设置非对称坐标轴范围（与原脚本逻辑一致）
negative_range <- ceiling(abs(negative_t_min) / 2) * 2
positive_range_compressed <- ceiling(positive_t_max * 0.5)  # 压缩正轴
positive_range <- ceiling(positive_range_compressed / 2) * 2

# 确保正轴范围足够大
if (positive_range < max(pathway_differences$t_value)) {
  positive_range <- ceiling(max(pathway_differences$t_value) * 1.1)
}

# 创建颜色分组逻辑（与原脚本完全一致）
pathway_differences$color_group <- ifelse(
  pathway_differences$is_weakest,  # 如果是t绝对值最小的20个
  "Weakest",  # 标记为Weakest（灰色）
  ifelse(
    pathway_differences$fdr < 0.05,  # 否则检查是否显著
    ifelse(
      pathway_differences$t_value > 0,
      "Higher in High-EMT",
      "Higher in Low-EMT"
    ),
    "Weakest"  # 不显著的也标记为Weakest
  )
)

# 检查颜色分组结果
cat("\nColor group distribution:\n")
color_counts <- table(pathway_differences$color_group)
for (group in names(color_counts)) {
  cat(sprintf("  %s: %d pathways\n", group, color_counts[group]))
}

# 颜色映射（与原脚本完全一致的颜色）
color_values <- c(
  "Higher in High-EMT" = "#4682B4",   # Steel blue
  "Higher in Low-EMT" = "#32CD32",    # Lime green
  "Weakest" = "#D3D3D3"               # Light gray
)

# 创建条形图（保持原脚本的绘图风格）
diff_bar_plot <- ggplot(
  pathway_differences, 
  aes(x = reorder(Description, t_value), 
      y = t_value, 
      fill = color_group)
) +
  geom_bar(
    stat = "identity", 
    width = 0.75, 
    color = ifelse(pathway_differences$color_group == "Weakest", 
                   "gray70", "gray50"),
    linewidth = ifelse(pathway_differences$color_group == "Weakest", 0.2, 0.5)
  ) +
  geom_hline(
    yintercept = 0, 
    linetype = "solid", 
    color = "gray40", 
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = color_values,
    name = NULL,
    breaks = c("Higher in High-EMT", "Higher in Low-EMT", "Weakest")
  ) +
  scale_y_continuous(
    breaks = function(x) {
      # 智能生成刻度（与原脚本逻辑一致）
      neg_breaks <- seq(0, -negative_range, by = -5)
      pos_breaks <- seq(0, positive_range, by = 2)
      unique(sort(c(neg_breaks, pos_breaks)))
    },
    limits = c(-negative_range, positive_range),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "t-value"
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(
      size = 10,
      color = ifelse(pathway_differences$color_group == "Weakest", 
                     "gray50", "black")
    ),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(
      face = "bold", 
      size = 12, 
      margin = margin(t = 10)
    ),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(12, "pt"),
    panel.grid.major = element_line(
      color = "gray95", 
      linewidth = 0.5
    ),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = "white", 
      color = NA
    ),
    panel.background = element_rect(
      fill = "white", 
      color = NA
    ),
    plot.margin = margin(30, 40, 30, 30)
  )

# 保存图形（使用原脚本的文件名）
diff_barplot_file <- file.path(fig_dir, "GSVA_pathway_activity_differences_optimized.pdf")
ggsave(
  diff_barplot_file,
  diff_bar_plot,
  width = 20,
  height = 16,
  dpi = 300
)

cat(sprintf("\n通路活性差异条形图已保存: %s\n", basename(diff_barplot_file)))