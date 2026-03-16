# Step 3-0: GEO data
rm(list = ls())
gc()

library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)

# GSE220095数据处理
setwd('/home/datahup/xzp/paper/PCa/data/raw_data_GEO/GSE220095')
#system("tar -xvf GSE220095_RAW.tar")
files <- list.files(pattern = "htseq_counts.txt.gz$")

count_list <- lapply(files, function(f) {
  dat <- read.table(gzfile(f), header = FALSE, stringsAsFactors = FALSE)
  colnames(dat) <- c("gene_id", sub("_RIB.*", "", f))  # 提取样本名
  return(dat)
})

count_matrix <- reduce(count_list, full_join, by = "gene_id")
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[ , -1]
count_matrix <- as.data.frame(lapply(count_matrix, as.integer), row.names = rownames(count_matrix))
head(count_matrix)

# ID映射函数
process_gene_annotation <- function(gene_ids) {
  # 判断基因类型
  is_ensembl <- grepl("^ENSG", gene_ids)
  is_xloc <- grepl("^XLOC", gene_ids)
  
  annotation_df <- data.frame(
    gene_id = gene_ids,
    gene_type = ifelse(is_ensembl, "Ensembl", 
                       ifelse(is_xloc, "XLOC_novel", "Other")),
    symbol = NA,
    entrezid = NA,
    stringsAsFactors = FALSE
  )
  
  # 处理 Ensembl 基因，去掉版本号和后缀
  if (sum(is_ensembl) > 0) {
    clean_ids <- str_replace(gene_ids[is_ensembl], "\\..*|_.*$", "")
    
    symbols <- mapIds(org.Hs.eg.db,
                      keys = clean_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
    
    entrez_ids <- mapIds(org.Hs.eg.db,
                         keys = clean_ids,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    
    annotation_df$symbol[is_ensembl] <- symbols
    annotation_df$entrezid[is_ensembl] <- entrez_ids
  }
  
  # XLOC 新基因保持原样或添加前缀
  annotation_df$symbol[is_xloc] <- paste0("NOVEL_", annotation_df$gene_id[is_xloc])
  
  # 对其他未匹配的基因用原始 gene_id 填充
  na_idx <- which(is.na(annotation_df$symbol))
  if(length(na_idx) > 0){
    annotation_df$symbol[na_idx] <- annotation_df$gene_id[na_idx]
  }
  
  return(annotation_df)
}
annotation_df <- process_gene_annotation(rownames(count_matrix))

mapped_symbols <- annotation_df$symbol[match(rownames(count_matrix), annotation_df$gene_id)]
rownames(count_matrix) <- make.unique(mapped_symbols)  # 自动生成唯一行名
count_matrix_t <- t(count_matrix)
rownames(count_matrix_t) <- colnames(count_matrix)
colnames(count_matrix_t) <- rownames(count_matrix)
count_matrix <- t(count_matrix)

cat("行数（样本数）:", nrow(count_matrix), "\n")         
cat("列数（基因数）:", ncol(count_matrix), "\n")
head(rownames(count_matrix), 5)
head(colnames(count_matrix), 5)


# 提取分期信息
lines <- readLines("../GSE220095_series_matrix.txt")
char_lines_all <- lines[grepl("^!Sample_characteristics_ch1", lines)]

t_stage_line <- char_lines_all[grepl("pathological t_stage", char_lines_all)]
t_stage_values <- str_split(t_stage_line, "\t")[[1]][-1]  
t_stage_values <- str_replace_all(t_stage_values, '"', '') 

n_stage_line <- char_lines_all[grepl("pathological n_stage", char_lines_all, ignore.case = TRUE)]
n_stage_values <- str_split(n_stage_line, "\t")[[1]][-1]
n_stage_values <- str_replace_all(n_stage_values, '"', '')

gleason_line <- char_lines_all[grepl("gleason score", char_lines_all, ignore.case = TRUE)]
gleason_values <- str_split(gleason_line, "\t")[[1]][-1]
gleason_values <- str_replace_all(gleason_values, '"', '')

psa_line <- char_lines_all[grepl("total psa", char_lines_all, ignore.case = TRUE)]
psa_values <- str_split(psa_line, "\t")[[1]][-1]
psa_values <- str_replace_all(psa_values, '"', '')

bcr_line <- char_lines_all[grepl("biochemical relapse.*bcr.*: no|biochemical relapse.*bcr.*: yes", char_lines_all, ignore.case = TRUE)]
bcr_values <- str_split(bcr_line, "\t")[[1]][-1]
bcr_values <- str_replace_all(bcr_values, '"', '')

# 提取time to bcr信息
time_to_bcr_line <- char_lines_all[grepl("time to_bcr", char_lines_all, ignore.case = TRUE)]
time_to_bcr_values <- str_split(time_to_bcr_line, "\t")[[1]][-1]
time_to_bcr_values <- str_replace_all(time_to_bcr_values, '"', '')

# 提取时间数值
time_to_bcr_numeric <- str_extract(time_to_bcr_values, "[0-9]+\\.?[0-9]*")

gsm_line <- lines[grepl("^!Sample_geo_accession", lines)]
gsm_ids <- str_split(gsm_line, "\t")[[1]][-1]
gsm_ids <- str_replace_all(gsm_ids, '"', '')

clinical_df <- data.frame(
  GSM = gsm_ids,
  t_stage = str_extract(t_stage_values, "pT[0-9]+[a-z]?"),
  n_stage = str_extract(n_stage_values, "pN[0-9]+"),
  gleason_score = str_extract(gleason_values, "[0-9]+=[0-9]+\\+[0-9]+"),
  psa_level = str_extract(psa_values, "[0-9]+\\.?[0-9]*"),
  bcr_status = str_extract(bcr_values, "yes|no"),
  time_to_bcr = as.numeric(time_to_bcr_numeric),
  stringsAsFactors = FALSE
)

# 数据处理
clinical_df$t_stage <- str_replace(clinical_df$t_stage, "^p", "")
clinical_df$n_stage <- str_replace(clinical_df$n_stage, "^p", "")
clinical_df$psa_level <- as.numeric(clinical_df$psa_level)
clinical_df$gleason_total <- str_extract(clinical_df$gleason_score, "^[0-9]+")

# 为生存分析创建事件指示变量
# 如果bcr_status为"yes"，则事件发生（event = 1）；如果为"no"，则是删失数据（event = 0）
clinical_df$bcr_event <- ifelse(clinical_df$bcr_status == "yes", 1, 0)

# 检查数据
cat("BCR状态分布:\n")
table(clinical_df$bcr_status)
cat("\n事件指示变量分布:\n")
table(clinical_df$bcr_event)
cat("\n时间变量的统计摘要:\n")
summary(clinical_df$time_to_bcr)

# 查看前几行数据
head(clinical_df)

# 可选：创建适合生存分析的简化数据框
survival_df <- data.frame(
  GSM = clinical_df$GSM,
  time = clinical_df$time_to_bcr,
  event = clinical_df$bcr_event,
  stringsAsFactors = FALSE
)

# 查看生存数据
head(survival_df)

save(annotation_df, file = "annotation_df.RData")
save(count_matrix, file = "count_matrix.RData")
save(clinical_df, file = "clinical_df.RData")
save(survival_df, file = "survival_df.RData")
rm(list = ls())
gc()

# GSE62116数据处理
library(huex10sttranscriptcluster.db)
library(oligo)

setwd('/home/datahup/xzp/paper/PCa/data/raw_data_GEO/GSE62116')
#system("tar -xvf GSE62116_RAW.tar")
#system("gunzip -k *.gz")

cel_files <- list.celfiles(".", full.names = TRUE)
raw_data <- read.celfiles(cel_files)

eset <- rma(raw_data)
expr_matrix <- exprs(eset)
dim(expr_matrix)

# 获取注释表
probe_ids <- rownames(expr_matrix)
gene_symbols <- mapIds(
  huex10sttranscriptcluster.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# 添加注释
expr_matrix_annot <- expr_matrix
rownames(expr_matrix_annot) <- gene_symbols

# 去掉没有基因符号的行
expr_matrix_annot <- expr_matrix_annot[!is.na(rownames(expr_matrix_annot)), ]

# 如果多个探针对应同一基因，取平均值
expr_matrix_annot <- aggregate(. ~ rownames(expr_matrix_annot), data = as.data.frame(expr_matrix_annot), FUN = mean)
rownames(expr_matrix_annot) <- expr_matrix_annot[, 1]
expr_matrix_annot <- expr_matrix_annot[, -1]
colnames(expr_matrix_annot) <- sub("_.*", "", colnames(expr_matrix_annot))
dim(expr_matrix_annot)
head(expr_matrix_annot[, 1:5])

write.csv(expr_matrix_annot, file = "GSE62116_expression_matrix_annotated.csv")

# 提取临床信息
expression_matrix <- read.csv("GSE62116_expression_matrix_annotated.csv")
gse <- getGEO(filename = "GSE62116_family.soft")

extract_simple_clinical <- function(gse_object) {
  clinical_data <- data.frame()
  
  for(i in 1:length(gse_object@gsms)) {
    gsm <- gse_object@gsms[[i]]
    
    # 提取简化临床信息
    sample_info <- data.frame(
      geo_accession = gsm@header$geo_accession,
      pathgs = if(!is.null(gsm@header$characteristics_ch1)) {
        char <- gsm@header$characteristics_ch1
        pathgs_char <- char[grep("pathgs:", char)]
        if(length(pathgs_char) > 0) gsub("pathgs: ", "", pathgs_char) else NA
      } else NA,
      
      preop_psa = if(!is.null(gsm@header$characteristics_ch1)) {
        char <- gsm@header$characteristics_ch1
        psa_char <- char[grep("preop_psa:", char)]
        if(length(psa_char) > 0) as.numeric(gsub("preop_psa: ", "", psa_char)) else NA
      } else NA,
      
      age = if(!is.null(gsm@header$characteristics_ch1)) {
        char <- gsm@header$characteristics_ch1
        age_char <- char[grep("age:", char)]
        if(length(age_char) > 0) as.numeric(gsub("age: ", "", age_char)) else NA
      } else NA,
      
      pstage = if(!is.null(gsm@header$characteristics_ch1)) {
        char <- gsm@header$characteristics_ch1
        stage_char <- char[grep("pstage:", char)]
        if(length(stage_char) > 0) gsub("pstage: ", "", stage_char) else NA
      } else NA,
      
      stringsAsFactors = FALSE
    )
    
    clinical_data <- rbind(clinical_data, sample_info)
  }
  
  return(clinical_data)
}

clinical_simple <- extract_simple_clinical(gse)
print(clinical_simple)
write.csv(clinical_simple, "GSE62116_clinical_simple.csv", row.names = FALSE)
