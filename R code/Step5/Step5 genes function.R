### Step5 genes function explore
rm(list = ls())
gc()

setwd('/home/datahup/xzp/paper/PCa/data/genes function')

# 创建图片保存目录
fig_dir <- '/home/datahup/xzp/paper/PCa/fig/genes function'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat("图片将保存到:", fig_dir, "\n")

# 1. 加载数据和EMT评分
selected_genes_df <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/core_emt_genes.csv")
selected_genes <- selected_genes_df$gene
selected_coef <- selected_genes_df$coefficient
names(selected_coef) <- selected_genes

cat("建模选中的EMT基因:\n")
print(selected_genes_df)

# 加载三组数据的EMT评分
emt_scores_all <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/emt_scores.csv")
geo_scores <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/geo_validation_scores.csv")

# 加载三组数据的表达矩阵
tcga_expr <- readRDS("../AUcell-Model/TCGA_PRAD_tpm.rds")
gse62116_expr <- read.csv("../raw_data_GEO/GSE62116/GSE62116_expression_matrix_annotated.csv", row.names = 1)
geo_final_data <- readRDS("../TCGA GEO-Model/seed_125_lasso_results/geo_validation_data.rds")

# 2. 准备三组数据的分析数据
prepare_analysis_data <- function(expr_data, emt_scores, data_name) {
  gene_expr <- as.data.frame(t(expr_data[selected_genes, ]))
  gene_expr$sample_id <- rownames(gene_expr)
  
  analysis_data <- gene_expr %>%
    inner_join(emt_scores, by = "sample_id")
  
  cat(data_name, "分析数据样本数:", nrow(analysis_data), "\n")
  return(analysis_data)
}

# GSE62116数据
gse62116_scores <- emt_scores_all %>% 
  filter(data_source == "GSE62116") %>%
  dplyr::select(sample_id, emt_score)

gse62116_analysis <- prepare_analysis_data(gse62116_expr, gse62116_scores, "GSE62116")

# GSE220095数据
geo_validation_scores <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/geo_validation_scores.csv")
geo_analysis <- geo_final_data %>%
  inner_join(geo_validation_scores, by = "sample_id")

cat("GSE220095分析数据样本数:", nrow(geo_analysis), "\n")

# TCGA数据
available_genes_tcga <- selected_genes[selected_genes %in% rownames(tcga_expr)]
cat("在TCGA中找到的EMT基因数:", length(available_genes_tcga), "/", length(selected_genes), "\n")

tcga_selected_genes <- tcga_expr[available_genes_tcga, ]
tcga_X <- as.data.frame(t(tcga_selected_genes))
selected_coef_available <- selected_coef[available_genes_tcga]
tcga_emt_scores <- as.numeric(as.matrix(tcga_X) %*% selected_coef_available)

tcga_analysis <- tcga_X
tcga_analysis$sample_id <- rownames(tcga_X)
tcga_analysis$emt_score <- tcga_emt_scores
tcga_analysis$EMTGroup <- ifelse(tcga_emt_scores > median(tcga_emt_scores), "High_EMT", "Low_EMT")

cat("TCGA分析数据样本数:", nrow(tcga_analysis), "\n")

# 3. 生存分析部分 
survival_data <- geo_analysis %>%
  mutate(
    # 事件指示：yes=1（事件发生），no=0（删失）
    bcr_event = ifelse(bcr_status == "yes", 1, 0),
    bcr_time = time_to_bcr
  ) %>%
  filter(!is.na(bcr_time))

cat("GSE220095生存分析有效样本数:", nrow(survival_data), "\n")
cat("BCR事件分布:\n")
print(table(survival_data$bcr_event))
cat("事件率:", round(mean(survival_data$bcr_event) * 100, 1), "%\n")

if(nrow(survival_data) > 0) {
  cat("生存时间统计:\n")
  cat("中位生存时间:", round(median(survival_data$bcr_time, na.rm = TRUE), 1), "月\n")
  cat("生存时间范围:", round(range(survival_data$bcr_time, na.rm = TRUE), 1), "月\n")
  
  # 使用与原始代码相同的KM生存分析函数
  perform_gene_km_analysis <- function(data, gene) {
    gene_data <- data[!is.na(data[[gene]]), ]
    
    # 使用中位数分组
    median_value <- median(gene_data[[gene]], na.rm = TRUE)
    gene_level <- ifelse(gene_data[[gene]] < median_value, "Low", "High")
    gene_level <- factor(gene_level, levels = c("Low", "High"))
    gene_data$gene_level <- gene_level
    
    # 检查分组是否都有数据
    group_counts <- table(gene_data$gene_level)
    cat("  ", gene, "分组样本数: Low=", group_counts["Low"], ", High=", group_counts["High"], "\n")
    
    # 检查每组事件数是否足够（至少2个事件）
    event_summary <- gene_data %>%
      group_by(gene_level) %>%
      summarise(events = sum(bcr_event))
    
    if(any(event_summary$events < 2)) {
      cat("  ", gene, ": 某些分组事件数不足2个，跳过生存分析\n")
      return(NULL)
    }
    
    # 生存分析 - 使用与原始代码相同的参数
    sfit <- survfit(Surv(bcr_time, bcr_event) ~ gene_level, data = gene_data)
    
    # 绘制KM曲线 - 使用与原始代码相同的风格
    km_plot <- ggsurvplot(
      sfit,
      data = gene_data,
      palette = "jco",  # 使用原始代码的配色
      conf.int = TRUE,
      conf.int.style = "step",
      pval = TRUE,
      pval.method = TRUE,
      risk.table = FALSE,  # 原始代码中设置为FALSE
      legend = c(0.85, 0.85),  # 原始代码中的位置
      legend.title = gene,
      legend.labs = c("Low", "High"),
      title = paste("BCR Survival Curve for", gene),  # 原始代码的标题
      xlab = "Time to BCR (Months)",  # 修改为BCR
      surv.median.line = "hv",
      ggtheme = theme_bw(base_size = 12)  # 原始代码的主题
    )
    
    return(km_plot)
  }
  
  # 对每个基因进行生存分析
  km_plots <- list()
  for (gene in selected_genes) {
    cat("正在分析基因:", gene, "\n")
    km_result <- perform_gene_km_analysis(survival_data, gene)
    
    if(!is.null(km_result)) {
      km_plots[[gene]] <- km_result
      
      # 保存单个基因的KM图 - 使用原始代码的尺寸
      ggsave(file.path(fig_dir, paste0("KM_", gene, ".pdf")),
             km_result$plot, width = 6.5, height = 6)
    }
  }
  
  # 组合所有KM图 - 使用原始代码的组合方式
  if(length(km_plots) > 0) {
    plot_list <- lapply(km_plots, function(x) x$plot)
    ncol <- ifelse(length(plot_list) <= 3, length(plot_list), 3)
    nrow <- ceiling(length(plot_list) / ncol)
    
    combined_plot <- plot_grid(plotlist = plot_list, ncol = ncol, nrow = nrow)
    ggsave(file.path(fig_dir, "Combined_KM_Curves.pdf"), 
           combined_plot, width = 6.5 * ncol, height = 6 * nrow)
    
    cat("KM生存分析完成，生成", length(km_plots), "个基因的生存曲线\n")
  } else {
    cat("没有生成有效的生存曲线\n")
  }
} else {
  cat("没有足够的生存数据进行分析\n")
}

# 4. 基因表达与EMT评分相关性分析
# 合并所有样本数据
combined_data <- bind_rows(
  tcga_analysis %>% mutate(dataset = "TCGA"),
  gse62116_analysis %>% mutate(dataset = "GSE62116"),
  geo_analysis %>% mutate(dataset = "GSE220095")
) %>%
  mutate(
    EMT_Group = ifelse(emt_score > median(emt_score, na.rm = TRUE), "High_EMT", "Low_EMT"),
    EMT_Group = factor(EMT_Group, levels = c("Low_EMT", "High_EMT"))
  )

cat("合并数据总样本数:", nrow(combined_data), "\n")
cat("EMT分组分布:\n")
print(table(combined_data$EMT_Group))

# 可视化
create_correlation_plot <- function(data, gene) {
  # 计算相关性
  cor_test <- cor.test(data$emt_score, data[[gene]], method = "pearson")
  cor_coef <- round(cor_test$estimate, 3)
  cor_pvalue <- ifelse(cor_test$p.value < 0.001, "<0.001", 
                       sprintf("%.3f", cor_test$p.value))
  
  # 创建散点图 - 用颜色区分高低EMT组
  p <- ggplot(data, aes(x = emt_score, y = .data[[gene]], color = EMT_Group)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black", fill = "lightgray") +
    scale_color_manual(values = c("Low_EMT" = "#2E8B57", "High_EMT" = "#CD5C5C"),
                       name = "EMT Group") +
    labs(
      title = gene,
      x = "EMT Score",
      y = "Expression Level",
      subtitle = paste("r =", cor_coef, ", p =", cor_pvalue)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 18, colour = "black", face = "bold", 
                                vjust = 1.9, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

# 生成所有基因的相关性图
cor_plots <- list()
for (i in 1:length(selected_genes)) {
  gene <- selected_genes[i]
  show_y_label <- (i == 1)  # 只有第一个图显示y轴标签
  
  p <- create_correlation_plot(combined_data, gene)
  if(!show_y_label) {
    p <- p + theme(axis.title.y = element_blank())
  }
  cor_plots[[i]] <- p
}

# 组合相关性图
combined_cor_plots <- plot_grid(plotlist = cor_plots, nrow = 2, ncol = 3)
ggsave(file.path(fig_dir, "Gene_EMT_Correlation.pdf"), 
       combined_cor_plots, width = 16, height = 10, dpi = 300)

# 5. 高低EMT组基因表达差异分
gene_data_list <- list()
for(gene in selected_genes) {
  # 正确提取基因表达数据
  gene_expr <- combined_data[[gene]]
  gene_df <- data.frame(
    Expression = gene_expr,
    Gene = rep(gene, length(gene_expr)),
    EMT_Group = combined_data$EMT_Group,
    sample_id = combined_data$sample_id
  )
  gene_data_list[[gene]] <- gene_df
}

df_combined <- do.call(rbind, gene_data_list)
rownames(df_combined) <- NULL

head(df_combined)
cat("组合数据样本数:", nrow(df_combined), "\n")
cat("每个基因的样本数:", nrow(df_combined) / length(selected_genes), "\n")

# 设置基因顺序
df_combined$Gene <- factor(df_combined$Gene, levels = selected_genes)

# 创建分半小提琴图 - 完全参考之前脚本
mypalette = c("#2E8B57", "#CD5C5C")  # 直接使用颜色向量

p_split_violin <- ggplot(df_combined, aes(x = Gene, y = Expression, fill = EMT_Group)) +
  # 分半小提琴 - 使用完全相同的参数
  geom_split_violin(alpha = .5, trim = FALSE, color = NA, width = 1) +
  # 均值点
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.2)) +
  # 误差线
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3, position = position_dodge(0.2)) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, color = 'black', hjust = 1),
    legend.position = 'top'
  ) +
  # 使用scale_fill_manual并明确指定颜色
  scale_fill_manual(values = mypalette, name = 'EMT Group') +
  labs(y = "Expression Level", x = NULL) +
  # 添加显著性标记 - 使用完全相同的参数
  stat_compare_means(
    aes(group = EMT_Group),
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "NS")
    ),
    label.y = max(df_combined$Expression, na.rm = TRUE) * 1.18,
    size = 5
  )

# 保存分半小提琴图
ggsave(file.path(fig_dir, "EMT_Group_Split_Violin.pdf"), 
       p_split_violin, width = 8, height = 6, dpi = 300)


# 创建传统的箱线图版本
create_gene_boxplot <- function(data, gene, show_y_label = TRUE) {
  
  # 计算p值
  high_data <- data[[gene]][data$EMT_Group == "High_EMT"]
  low_data <- data[[gene]][data$EMT_Group == "Low_EMT"]
  p_value <- t.test(high_data, low_data)$p.value
  p_label <- ifelse(p_value < 0.001, "***", 
                    ifelse(p_value < 0.01, "**",
                           ifelse(p_value < 0.05, "*", "NS")))
  
  p <- ggplot(data, aes(x = EMT_Group, y = .data[[gene]], fill = EMT_Group)) + 
    geom_boxplot(alpha = 0.7, size = 1.2, width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.2, aes(color = EMT_Group)) +
    scale_fill_manual(values = c("Low_EMT" = "#2E8B57", "High_EMT" = "#CD5C5C")) +
    scale_color_manual(values = c("Low_EMT" = "#2E8B57", "High_EMT" = "#CD5C5C")) +
    geom_signif(
      comparisons = list(c("Low_EMT", "High_EMT")),
      map_signif_level = FALSE,
      annotations = p_label,
      tip_length = 0.01,
      size = 1,
      textsize = 6,
      y_position = max(data[[gene]], na.rm = TRUE) * 1.1
    ) +
    labs(
      title = gene,
      y = ifelse(show_y_label, "Expression Level", ""),
      x = ""
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, colour = "black", face = "bold", 
                                vjust = 1.9, hjust = 0.5),
      axis.title.y = element_text(size = 12, color = "black", face = "bold.italic", 
                                  vjust = 1.9, hjust = 0.5, angle = 90),
      axis.text.x = element_text(size = 13, color = "black", face = "bold",
                                 vjust = 0.5, hjust = 0.5, angle = 0),
      axis.text.y = element_text(size = 13, color = "black", face = "bold.italic",
                                 vjust = 0.5, hjust = 0.5, angle = 0)
    )
  
  return(p)
}

# 生成所有基因的箱线图
gene_plots <- list()
for (i in 1:length(selected_genes)) {
  gene <- selected_genes[i]
  show_y_label <- (i == 1)  # 只有第一个图显示y轴标签
  gene_plots[[i]] <- create_gene_boxplot(combined_data, gene, show_y_label)
}

# 组合所有基因的箱线图
combined_boxplots <- plot_grid(plotlist = gene_plots, nrow = 1)
ggsave(file.path(fig_dir, "EMT_Group_Gene_Expression.pdf"), 
       combined_boxplots, width = 18, height = 6, dpi = 300)
