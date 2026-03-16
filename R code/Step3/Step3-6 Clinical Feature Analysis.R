# Step 3-6: Clinical Feature Analysis
rm(list = ls())
gc()

library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd('/home/datahup/xzp/paper/PCa/data/Clinical Feature Analysis')

# 1. 加载GSE220095数据
load("/home/datahup/xzp/paper/PCa/data/raw_data_GEO/GSE220095/clinical_df.RData")
geo_validation_scores <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/geo_validation_scores.csv")

emt_scores_only <- geo_validation_scores %>%
  select(sample_id, emt_score)
gse220095_data <- clinical_df %>%
  left_join(emt_scores_only, by = c("GSM" = "sample_id")) %>%
  filter(!is.na(emt_score))

cat("GSE220095有效样本数:", nrow(gse220095_data), "\n")

# 数据清理和变量处理 - 基于参考脚本模式
gse220095_clean <- gse220095_data %>%
  mutate(
    # T分期处理 - 使用clinical_df中的t_stage，保持原始分类
    T_Stage = case_when(
      grepl("T1", t_stage, ignore.case = TRUE) ~ "T1",
      grepl("T2", t_stage, ignore.case = TRUE) ~ "T2", 
      grepl("T3", t_stage, ignore.case = TRUE) ~ "T3",
      grepl("T4", t_stage, ignore.case = TRUE) ~ "T4",
      TRUE ~ "Unknown"
    ),
    
    # N分期处理 - 保持原始分类
    N_Stage = case_when(
      grepl("N0", n_stage, ignore.case = TRUE) ~ "N0",
      grepl("N1", n_stage, ignore.case = TRUE) ~ "N1",
      grepl("N2", n_stage, ignore.case = TRUE) ~ "N2",
      TRUE ~ "Unknown"
    ),
    
    # Gleason评分处理 - 从"X=Y+Z"格式中提取X作为总分
    Gleason_Total = as.numeric(str_extract(gleason_score, "^[0-9]+")),
    
    # Gleason评分分组 - 类似年龄的二分法
    Gleason_Group = ifelse(Gleason_Total <= 7, "≤7", "≥8"),
    
    # PSA水平分组 - 类似年龄的多分组
    PSA_Group = case_when(
      as.numeric(psa_level) < 10 ~ "<10",
      as.numeric(psa_level) >= 10 & as.numeric(psa_level) < 20 ~ "10-20",
      as.numeric(psa_level) >= 20 ~ "≥20", 
      TRUE ~ "Unknown"
    ),
    
    # EMT评分分组 - 类似风险评分的中位数分组
    EMT_Group = ifelse(emt_score > median(emt_score, na.rm = TRUE), "High", "Low"),
    
    # 数据集标识
    dataset = "GSE220095",
    
    # 为后续分析准备数值型变量
    PSA_Numeric = as.numeric(psa_level),
    
    # 保留原始样本ID
    sample_id = GSM
  )
cat("总样本数:", nrow(gse220095_clean), "\n")


# 2. 加载GSE62116数据
gse62116_dir <- '/home/datahup/xzp/paper/PCa/data/raw_data_GEO/GSE62116'
gse62116_expr <- read.csv(file.path(gse62116_dir, "GSE62116_expression_matrix_annotated.csv"), row.names = 1)
gse62116_clinical <- read.csv(file.path(gse62116_dir, "GSE62116_clinical_simple.csv"))

# 读取GSE62116的EMT评分
all_emt_scores <- read.csv("../TCGA GEO-Model/seed_125_lasso_results/emt_scores.csv")
gse62116_emt_scores <- all_emt_scores %>%
  filter(data_source == "GSE62116") %>%
  select(sample_id, emt_score)

# 首先对临床数据进行去重
gse62116_clinical_dedup <- gse62116_clinical %>%
  distinct(geo_accession, .keep_all = TRUE)  # 按geo_accession去重，保留第一条

cat("GSE62116去重前临床数据样本数:", nrow(gse62116_clinical), "\n")
cat("GSE62116去重后临床数据样本数:", nrow(gse62116_clinical_dedup), "\n")

# 数据清理和变量处理 - 与GSE220095保持一致的格式
gse62116_clean <- gse62116_clinical_dedup %>%
  left_join(gse62116_emt_scores, by = c("geo_accession" = "sample_id")) %>%
  filter(!is.na(emt_score)) %>%
  mutate(
    # T分期处理 - 基于pstage字段
    T_Stage = case_when(
      grepl("T1", pstage, ignore.case = TRUE) ~ "T1",
      grepl("T2", pstage, ignore.case = TRUE) ~ "T2", 
      grepl("T3", pstage, ignore.case = TRUE) ~ "T3",
      grepl("T4", pstage, ignore.case = TRUE) ~ "T4",
      TRUE ~ "Unknown"
    ),
    
    # N分期处理 - 从pstage中提取
    N_Stage = case_when(
      grepl("N0", pstage, ignore.case = TRUE) ~ "N0",
      grepl("N1", pstage, ignore.case = TRUE) ~ "N1", 
      grepl("N\\+", pstage, ignore.case = TRUE) ~ "N+",
      grepl("N2", pstage, ignore.case = TRUE) ~ "N2",
      TRUE ~ "Unknown"
    ),
    
    # Gleason评分处理 - 使用pathgs字段
    Gleason_Total = as.numeric(pathgs),
    Gleason_Group = ifelse(Gleason_Total <= 7, "≤7", "≥8"),
    
    # PSA水平处理 - 使用preop_psa字段
    PSA_Numeric = as.numeric(preop_psa),
    PSA_Group = case_when(
      PSA_Numeric < 10 ~ "<10",
      PSA_Numeric >= 10 & PSA_Numeric < 20 ~ "10-20",
      PSA_Numeric >= 20 ~ "≥20",
      TRUE ~ "Unknown"
    ),
    
    # 年龄处理 - 参考脚本中的年龄分组
    Age_Numeric = as.numeric(age),
    Age_Group = ifelse(Age_Numeric <= 60, '≤60', '>60'),
    
    # EMT评分分组
    EMT_Group = ifelse(emt_score > median(emt_score, na.rm = TRUE), "High", "Low"),
    
    # 数据集标识
    dataset = "GSE62116",
    
    # 保留原始样本ID
    sample_id = geo_accession,
    
    # 添加其他临床变量占位符
    M_Stage = "Unknown"  # GSE62116通常没有M分期信息
  )

cat("GSE62116有效样本数:", nrow(gse62116_clean), "\n")
cat("GSE62116数据预览:\n")
print(head(gse62116_clean))

# 3. 加载和处理TCGA数据
# 读取TCGA临床数据和EMT评分
tcga_clinical <- read.delim("/home/datahup/xzp/paper/PCa/data/raw_data_TCGA/clinical.tsv")

# 读取TCGA的EMT评分
tcga_emt_scores <- all_emt_scores %>%
  filter(data_source == "TCGA") %>%
  select(sample_id, emt_score)

# 首先处理TCGA临床数据去重
tcga_clinical_dedup <- tcga_clinical %>%
  mutate(
    patient_id = substr(cases.submitter_id, 1, 12),
    patient_id = gsub("\\.", "-", patient_id)
  ) %>%
  distinct(patient_id, .keep_all = TRUE)  # 按patient_id去重，保留第一个

cat("去重前TCGA临床数据样本数:", nrow(tcga_clinical), "\n")
cat("去重后TCGA临床数据样本数:", nrow(tcga_clinical_dedup), "\n")

# 数据清理和变量处理 - 与其他数据集保持一致的格式
tcga_clean <- tcga_clinical_dedup %>%
  left_join(tcga_emt_scores, by = c("patient_id" = "sample_id")) %>%
  filter(!is.na(emt_score)) %>%
  mutate(
    # T分期处理 - 使用ajcc_pathologic_t
    T_Stage = case_when(
      grepl("^T1", diagnoses.ajcc_pathologic_t) | 
        diagnoses.ajcc_pathologic_t %in% c("T1", "T1a", "T1b", "T1c") ~ "T1",
      grepl("^T2", diagnoses.ajcc_pathologic_t) | 
        diagnoses.ajcc_pathologic_t %in% c("T2", "T2a", "T2b", "T2c") ~ "T2", 
      grepl("^T3", diagnoses.ajcc_pathologic_t) | 
        diagnoses.ajcc_pathologic_t %in% c("T3", "T3a", "T3b") ~ "T3",
      grepl("^T4", diagnoses.ajcc_pathologic_t) | 
        diagnoses.ajcc_pathologic_t %in% c("T4") ~ "T4",
      TRUE ~ "Unknown"
    ),
    
    # N分期处理 - 使用ajcc_pathologic_n
    N_Stage = case_when(
      grepl("^N0", diagnoses.ajcc_pathologic_n) ~ "N0",
      grepl("^N1", diagnoses.ajcc_pathologic_n) ~ "N1",
      grepl("^N2", diagnoses.ajcc_pathologic_n) ~ "N2",
      diagnoses.ajcc_pathologic_n == "'--" ~ "Unknown",
      is.na(diagnoses.ajcc_pathologic_n) ~ "Unknown",
      TRUE ~ "Unknown"
    ),
    
    # M分期处理 - 使用ajcc_pathologic_m
    M_Stage = case_when(
      grepl("^M0", diagnoses.ajcc_pathologic_m) ~ "M0",
      grepl("^M1", diagnoses.ajcc_pathologic_m) ~ "M1",
      diagnoses.ajcc_pathologic_m == "'--" ~ "Unknown",
      is.na(diagnoses.ajcc_pathologic_m) ~ "Unknown",
      TRUE ~ "Unknown"
    ),
    
    # Gleason评分处理 - 使用gleason_score
    Gleason_Total = as.numeric(ifelse(diagnoses.gleason_score == "'--", NA, diagnoses.gleason_score)),
    Gleason_Group = ifelse(Gleason_Total <= 7, "≤7", "≥8"),
    
    # 年龄处理 - 使用demographic.age_at_index
    Age_Numeric = as.numeric(demographic.age_at_index),
    Age_Group = ifelse(Age_Numeric <= 60, '≤60', '>60'),
    
    # PSA水平 - TCGA通常没有PSA数据，设为NA
    PSA_Numeric = NA,
    PSA_Group = "Unknown",
    
    # EMT评分分组
    EMT_Group = ifelse(emt_score > median(emt_score, na.rm = TRUE), "High", "Low"),
    
    # 数据集标识
    dataset = "TCGA",
    
    # 保留原始样本ID
    sample_id = patient_id
  )

cat("TCGA有效样本数:", nrow(tcga_clean), "\n")
cat("TCGA数据预览:\n")
print(head(tcga_clean %>% select(sample_id, T_Stage, N_Stage, M_Stage, Gleason_Total, Age_Numeric, emt_score)))

# 4. 合并三组数据
combined_clinical_data <- bind_rows(
  gse220095_clean %>% 
    mutate(
      M_Stage = "Unknown",           # 添加M分期
      Age_Numeric = NA,              # 添加年龄数值
      Age_Group = "Unknown"          # 添加年龄分组
    ) %>% 
    select(sample_id, dataset, T_Stage, N_Stage, M_Stage, 
           Gleason_Total, Gleason_Group, Age_Numeric, Age_Group,
           PSA_Numeric, PSA_Group, emt_score, EMT_Group),
  
  gse62116_clean %>% 
    select(sample_id, dataset, T_Stage, N_Stage, M_Stage,
           Gleason_Total, Gleason_Group, Age_Numeric, Age_Group,
           PSA_Numeric, PSA_Group, emt_score, EMT_Group),
  
  tcga_clean %>% 
    select(sample_id, dataset, T_Stage, N_Stage, M_Stage,
           Gleason_Total, Gleason_Group, Age_Numeric, Age_Group,
           PSA_Numeric, PSA_Group, emt_score, EMT_Group)
)

cat("合并后总样本数:", nrow(combined_clinical_data), "\n")
cat("各数据集样本分布:\n")
print(table(combined_clinical_data$dataset))

# 5. 临床特征分析可视化 - 参考脚本风格
library(see)
library(ggplot2)
library(ggpubr)
library(cowplot)


# 设置图形保存目录
fig_dir <- '/home/datahup/xzp/paper/PCa/fig/Clinical Feature Analysis'
if(!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# 创建数据处理和绘图函数
create_proper_half_violin_plot <- function(data, x_var, title = "", filename = "") {
  
  # 过滤数据：跳过NA和Unknown
  plot_data <- data %>% 
    filter(!is.na(.data[[x_var]]), 
           .data[[x_var]] != "Unknown",
           !is.na(emt_score))
  
  cat("分析变量:", x_var, "- 有效样本:", nrow(plot_data), "\n")
  
  if(nrow(plot_data) < 5) {
    cat("跳过", x_var, "- 数据量不足\n")
    return(NULL)
  }
  
  # 获取分组数量
  n_groups <- length(unique(plot_data[[x_var]]))
  
  # 设置颜色 - 完全参考原脚本
  if(n_groups == 2) {
    fill_colors <- c("#D55E00A0", "#0072B2A0")
  } else if(n_groups == 3) {
    fill_colors <- c("#D55E00A0", "#0072B2A0", '#008080')
  } else if(n_groups == 4) {
    fill_colors <- c("#D55E00A0", "#0072B2A0", '#008080', '#2E8B57')
  } else {
    fill_colors <- scales::hue_pal()(n_groups)
  }
  
  # 计算统计检验p值
  if(n_groups == 2) {
    groups <- unique(plot_data[[x_var]])
    group1_data <- plot_data$emt_score[plot_data[[x_var]] == groups[1]]
    group2_data <- plot_data$emt_score[plot_data[[x_var]] == groups[2]]
    p_value <- t.test(group1_data, group2_data)$p.value
    p_label <- ifelse(p_value < 0.001, "***", 
                      ifelse(p_value < 0.01, "**",
                             ifelse(p_value < 0.05, "*", "NS")))
  } else {
    p_value <- kruskal.test(emt_score ~ get(x_var), data = plot_data)$p.value
    p_label <- paste("p =", format.pval(p_value, digits = 3))
  }
  
  # 创建图形 - 使用标准violin图但通过coord_flip实现横向效果
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = emt_score, 
                             fill = .data[[x_var]])) +
    # 使用violin图和boxplot的组合
    geom_violin(alpha = 0.7, trim = TRUE, width = 0.8, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, 
                 position = position_dodge(0.9)) +
    geom_jitter(aes(color = .data[[x_var]]), width = 0.1, shape = "|", size = 1, alpha = 0.6) +
    coord_flip() + 
    theme_light() +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    labs(title = paste(title, " (", p_label, ")"), x = "", y = "EMT Score") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
  
  # 保存图形
  if(filename != "") {
    ggsave(file.path(fig_dir, filename), p, width = 8, height = 4, dpi = 300)
    cat("已保存:", filename, "\n")
  }
  
  return(p)
}

# 重新运行分析 - 使用标准violin图
combined_final <- combined_clinical_data %>%
  mutate(
    # T分期重新分类：T1/T2 vs T3/T4
    T_Stage_Binary = case_when(
      T_Stage %in% c("T1", "T2") ~ "T1-T2",
      T_Stage %in% c("T3", "T4") ~ "T3-T4",
      TRUE ~ "Unknown"
    ),
    
    # N分期重新分类：N0 vs N1 (合并N1, N2, N+为N1)
    N_Stage_Binary = case_when(
      N_Stage == "N0" ~ "N0",
      N_Stage %in% c("N1", "N2", "N+") ~ "N1",
      TRUE ~ "Unknown"
    ),
    
    # 替换特殊字符避免编码警告
    PSA_Group = gsub("≥", ">=", PSA_Group),
    Gleason_Group = gsub("≤", "<=", Gleason_Group) %>% gsub("≥", ">=", .),
    Age_Group = gsub("≤", "<=", Age_Group)
  )

# 年龄分组分析
p_age <- create_proper_half_violin_plot(combined_final, "Age_Group",
                                        "EMT Score by Age Group", "Age_Analysis.pdf")

# Gleason评分分析
p_gleason <- create_proper_half_violin_plot(combined_final, "Gleason_Group",
                                            "EMT Score by Gleason Score", "Gleason_Analysis.pdf")

# T分期二分类分析
p_tstage <- create_proper_half_violin_plot(combined_final, "T_Stage_Binary",
                                           "EMT Score by T Stage", "T_Stage_Binary.pdf")

# N分期二分类分析
p_nstage <- create_proper_half_violin_plot(combined_final, "N_Stage_Binary",
                                           "EMT Score by N Stage", "N_Stage_Binary.pdf")

# PSA分析
p_psa <- create_proper_half_violin_plot(combined_final, "PSA_Group",
                                        "EMT Score by PSA Level", "PSA_Analysis.pdf")

# 数据集比较
p_dataset <- create_proper_half_violin_plot(combined_final, "dataset",
                                            "EMT Score by Dataset", "Dataset_Comparison.pdf")

# 生成组合图
plots_list <- list(p_age, p_gleason, p_tstage, p_nstage, p_psa)
plots_list <- plots_list[!sapply(plots_list, is.null)]

if(length(plots_list) > 0) {
  combined_plot <- plot_grid(plotlist = plots_list, nrow = 2)
  ggsave(file.path(fig_dir, "Combined_Clinical_Features_Final.pdf"), combined_plot, 
         width = 16, height = 12, dpi = 300)
  cat("已保存组合图: Combined_Clinical_Features_Final.pdf\n")
}
