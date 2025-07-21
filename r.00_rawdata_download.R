rm(list = ls())
gc()
setwd("/data/nas1/mengqingyao_OD/project/scRNA-seq")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
options(stringsAsFactors = FALSE)

library(GEOquery)
library(lance)
library(tidyverse)
library(rentrez)
library(R.utils)

# 数据处理 -------------------------------------------------------------------
GEO_data <- 'GSE163973'

options(timeout = 600)

if (!file.exists("GSE163973_RAW.tar")) {
  getGEOSuppFiles(GEO_data)
}

# 解压数据
print(file.exists("../00_rawdata/GSE163973/"))

if (!dir.exists("raw_data")) {
  dir.create("raw_data")
}

untar("../00_rawdata/GSE163973/GSE163973_RAW.tar", exdir = "raw_data")

# 解压单个gz文件
gz_files <- list.files("raw_data", pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, gunzip, overwrite = TRUE)

# 样本信息整理 ------------------------------------------------------------------
gse <- getGEO("GSE163973", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[1]]))

sample_info_clean <- sample_info %>%
  dplyr::select(geo_accession, characteristics_ch1.1, source_name_ch1) %>%
  mutate(
    sample_id = geo_accession,
    condition = case_when(
      grepl("NASH", characteristics_ch1.1) ~ "NASH",
      grepl("control", characteristics_ch1.1, ignore.case = TRUE) ~ "Control",
      TRUE ~ "Other")
  )

# 匹配文件名与样本ID
sample_files <- data.frame(
  file_name = list.files("raw_data", pattern = "norm_counts", recursive = TRUE),
  sample_id = gsub("_.*", "", list.files("raw_data", pattern = "norm_counts", recursive = TRUE))
)

sample_metadata <- merge(sample_info_clean, sample_files, by = "sample_id") %>% 
  filter(condition != "Other")

sample_metadata <- sample_metadata %>% 
  filter(characteristics_ch1.1 != "disease status: NASH with cirrhosis")

table(sample_metadata$condition)
