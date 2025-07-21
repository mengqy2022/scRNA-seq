rm(list = ls())
gc()
setwd("/data/nas1/mengqingyao_OD/project/scRNA-seq")
if (!dir.exists("./14_scRNA-seq_QC")){
  dir.create("./14_scRNA-seq_QC")
}
setwd("./14_scRNA-seq_QC")
options(stringsAsFactors = F)

#options(future.globals.maxSize = 10 * 1024^3)

# 加载包 ---------------------------------------------------------------------
library(Seurat)        # 单细胞分析
library(tidyverse)
library(patchwork)     # 图形排版
library(Matrix)        # 稀疏矩阵处理
library(R.utils)       # 文件解压
library(GEOquery)      # GEO数据下载

# 数据处理 -------------------------------------------------------------------
# 解压数据
print(file.exists("../00_rawdata/GSE163973/"))

if (!dir.exists("raw_data")) {
  dir.create("raw_data")
}

untar("../00_rawdata/GSE202379/GSE202379_RAW.tar", exdir = "raw_data")

# 解压单个gz文件
gz_files <- list.files("raw_data", pattern = "\\.gz$", full.names = TRUE)
sapply(gz_files, gunzip, overwrite = TRUE)

# 样本信息整理 ------------------------------------------------------------------
gse <- getGEO("GSE202379", GSEMatrix = TRUE)
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

# 数据读取与Seurat对象创建 ---------------------------------------------------------
seurat_list <- list()

for (i in 1:nrow(sample_metadata)) {
  sample_path <- file.path("./raw_data", sample_metadata$file_name[i])
  sample_name <- sample_metadata$sample_id[i]
  sample_group <- sample_metadata$condition[i]
  
  cat("Processing:", sample_name, "(", sample_group, ")...\n")
  
  # 读取数据
  data <- tryCatch(
    {
      counts <- read.csv(sample_path, row.names = 1)
      # 转换为稀疏矩阵以提高效率
      as(as.matrix(counts), "dgCMatrix")
    },
    error = function(e) {
      cat("Error reading", sample_name, ":", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(data)) next
  
  # 创建Seurat对象 (v5兼容方式)
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 添加元数据
  seurat_obj$sample <- sample_name
  seurat_obj$group <- sample_group
  
  # 检测线粒体基因
  mt_pattern <- ifelse(any(grepl("^MT-", rownames(seurat_obj))), "^MT-", "^mt-")
  if(any(grepl(mt_pattern, rownames(seurat_obj)))){
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  } else {
    seurat_obj[["percent.mt"]] <- 0
    warning("No mitochondrial genes found in sample ", sample_name)
  }
  
  # 核糖体基因比例
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]|^Rp[sl]")
  
  seurat_list[[sample_name]] <- seurat_obj
  cat("After filtering:", ncol(seurat_obj), "cells retained in", sample_name, "\n")
}

# 使用Seurat v5的合并方式
merged_seurat <- merge(seurat_list[[1]], 
                       y = seurat_list[-1],
                       add.cell.ids = names(seurat_list),
                       project = "merged_project")

# 查看数据层次
Layers(merged_seurat)

# 添加批次信息
merged_seurat$batch <- factor(gsub("_.*", "", names(seurat_list)))

# 设置绘图主题
theme_set(theme_classic(base_size = 12))
group_colors <- c("Control" = "#E41A1C", "NASH" = "#377EB8")

qc_violin <- VlnPlot(merged_seurat, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "group",
                     pt.size = 0,
                     ncol = 3,
                     cols = group_colors) &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))

ggsave("QC_raw.png", plot = qc_violin, width = 8, height = 6, dpi = 300)
ggsave("QC_raw.pdf", plot = qc_violin, width = 8, height = 6)

# 保存原始合并数据
saveRDS(merged_seurat, file = "GSE202379_merged_raw.rds")

# 质量控制过滤  -----------------------------------------------------------------
merged_seurat <- readRDS(file = "GSE202379_merged_raw.rds")

merged_seurat <- subset(merged_seurat,
                        subset = nFeature_RNA > 200 &    # 至少检测到200个基因
                          nFeature_RNA < 4500 &          # 最多检测到4000个基因
                          #percent.mt < 15 &              # 线粒体基因比例<15%
                          nCount_RNA < 8500)           # UMI数<10000

# 质量控制小提琴图
qc_violin <- VlnPlot(merged_seurat, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     group.by = "group",
                     pt.size = 0,
                     ncol = 3,
                     cols = group_colors) &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12))

# 保存图形
ggsave("QC.png", plot = qc_violin, width = 8, height = 6, dpi = 300)
ggsave("QC.pdf", plot = qc_violin, width = 8, height = 6)

# 数据标准化和高变异基因识别 -----------------------------------------------------------
Layers(merged_seurat)

split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)

# 检查所有样本的基因是否一致
gene_lists <- lapply(seurat_list, function(x) rownames(x))
all_genes <- Reduce(intersect, gene_lists)  # 取所有样本共有的基因
print(paste("共有基因数量:", length(all_genes)))

# 如果有基因丢失，重新合并时只保留共有基因
merged_seurat_1 <- merge(
  seurat_list[[1]][all_genes, ],
  y = lapply(seurat_list[-1], function(x) x[all_genes, ]),
  add.cell.ids = names(seurat_list))

# 标准化
merged_seurat_1 <- NormalizeData(merged_seurat_1)

# 识别高变异基因
merged_seurat <- FindVariableFeatures(merged_seurat_1, 
                                      selection.method = "vst", 
                                      normalization.method = "LogNormalize",
                                      nfeatures = 2000,
                                      verbose = FALSE)

# 标记Top10高变异基因
#top10_var_genes <- head(VariableFeatures(merged_seurat), 10)
top10_var_genes <-  HVFInfo(merged_seurat) %>%
  arrange(desc(variance.standardized)) %>%
  head(10) %>%
  rownames()

# 高变异基因可视化
var_feature_plot <- VariableFeaturePlot(merged_seurat,
                                        pt.size = 1.5,
                                        assay = "RNA") +
  theme(legend.position = "top") +
  labs(title = "Highly Variable Genes (Top 2000)")

var_feature_plot[["data"]] <- na.omit(var_feature_plot[["data"]])

var_feature_plot <- LabelPoints(plot = var_feature_plot,
                                points = top10_var_genes,
                                repel = TRUE)

# 保存可视化结果
ggsave("Highly_Variable_Genes.png", plot = var_feature_plot, 
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Highly_Variable_Genes.pdf", plot = var_feature_plot, 
       width = 8, height = 6)

# pdf(file.path("VariableFeatures.pdf"))
# plot <- VariableFeaturePlot(merged_seurat, 
#                             cols = c("gray", "red"), 
#                             pt.size = 1.5,        
#                             assay = "RNA")  
# plot + 
#   ggrepel::geom_text_repel(data = HVFInfo(merged_seurat)[top10_var_genes, ], 
#                   aes(mean, variance.standardized, 
#                       label = rownames(HVFInfo(merged_seurat)[top10_var_genes, ])),
#                   size = 4, 
#                   box.padding = 0.3, 
#                   max.overlaps = Inf) + 
#   ggtitle("Top 2000 Highly Variable Genes (HVGs)") + 
#   theme_classic(base_size = 12)
# dev.off()

# 保存处理后的数据
saveRDS(merged_seurat, file = "GSE202379_processed.rds")

# 数据标准化和PCA降维分析 ------------------------------------------------------
# 数据标准化（对高变异基因进行缩放）
#merged_seurat <- readRDS(file = "GSE202379_processed.rds")

merged_seurat <- ScaleData(merged_seurat, 
                           features = VariableFeatures(merged_seurat))

# PCA降维分析
merged_seurat <- RunPCA(merged_seurat, 
                        features = VariableFeatures(merged_seurat),
                        verbose = FALSE)

# 可视化PCA结果 - ElbowPlot
elbow_plot <- ElbowPlot(merged_seurat, ndims = 30) +
  #geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  labs(title = "PCA Elbow Plot (Identifying Significant PCs)",
       subtitle = "Red line indicates suggested cutoff (PC20)") +
  theme_classic()

ggsave("PCA_ElbowPlot.png", plot = elbow_plot, 
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("PCA_ElbowPlot.pdf", plot = elbow_plot, 
       width = 8, height = 6)

# Jackstraw分析
set.seed(0629)
merged_seurat <- JackStraw(merged_seurat, num.replicate = 100, 
                           dims = 30)
merged_seurat <- ScoreJackStraw(merged_seurat, 
                                dims = 1:30)

jackstraw_plot <- JackStrawPlot(merged_seurat, 
                                dims = 1:30) +
  labs(title = "Jackstraw Plot for Significant PCs") +
  theme_classic()

ggsave("Jackstraw_Plot.png", plot = jackstraw_plot, 
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Jackstraw_Plot.pdf", plot = jackstraw_plot, 
       width = 8, height = 6)

# # 确定显著主成分（P < 0.05）
# p.values <- merged_seurat@reductions$pca@jackstraw$overall.p.values
# sig_pcs <- which(p.values[,2] < 0.05)
# cat("Number of significant PCs (p < 0.05):", length(sig_pcs), "\n")
# 
# p_cluster <- DimPlot(merged_seurat,
#         #reduction = "pca",
#         group.by = "sample",  # 使用聚类结果
#         label = TRUE,                  # 显示聚类标签
#         label.size = 4,
#         pt.size = 1) +
#   ggtitle("PCA by Seurat Clusters") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position = "right") +
#   scale_color_discrete(name = "Cluster")
# 
# ggsave("PCA_cluster.png", plot = p_cluster,
#        width = 8, height = 6, dpi = 300, bg = "white")
# ggsave("PCA_cluster.pdf", plot = p_cluster,
#        width = 8, height = 6)

# 细胞聚类分析 ----------------------------------------------------------------
# 这段代码是在使用 Seurat 和 Harmony 对单细胞RNA测序（scRNA-seq）数据进行整合，
# 目的是消除批次效应（batch effect），使得不同样本（或实验批次）的数据能够合并分析。
merged_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by.vars = "orig.ident",
  verbose = FALSE)

merged_seurat<- RunUMAP(merged_seurat, 
                        reduction = "harmony", 
                        dims = 1:30,
                        reduction.name = "umap")

merged_seurat[["RNA"]]@layers$scale.data <- NULL

saveRDS(merged_seurat,"./merged_seurat_umap.rds")

merged_seurat <- FindNeighbors(
  merged_seurat, 
  reduction = "harmony",   
  dims = 1:30,            
  verbose = FALSE)

merged_seurat <- FindClusters(merged_seurat, resolution = seq(0,1,0.05))
saveRDS(merged_seurat, file = "./merged_seurat_unannotation.rds")

#merged_seurat <- readRDS(file = "./merged_seurat_unannotation.rds")
head(merged_seurat@meta.data$seurat_clusters)

merged_seurat@meta.data$seurat_clusters <- merged_seurat@meta.data$RNA_snn_res.0.2

# p_cluster <- DimPlot(merged_seurat,
#         reduction = "pca",
#         group.by = "seurat_clusters",  # 使用聚类结果
#         label = TRUE,                  # 显示聚类标签
#         label.size = 4,
#         pt.size = 1) +
#   ggtitle("PCA by Seurat Clusters") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position = "right") +
#   scale_color_discrete(name = "Cluster")
# 
# ggsave("PCA_cluster.png", plot = p_cluster,
#        width = 8, height = 6, dpi = 300, bg = "white")
# ggsave("PCA_cluster.pdf", plot = p_cluster,
#        width = 8, height = 6)

#  设置Seurat对象中的细胞主标识（active identity）为聚类结果
Idents(merged_seurat)<-merged_seurat@meta.data$seurat_clusters

umap_plot <- DimPlot(merged_seurat, 
                     reduction = "umap", 
                     label = TRUE, 
                     label.size = 4,
                     pt.size = 0.5) +
  labs(title = "UMAP Visualization of Cell Clusters") +
  theme_classic() #+
  #NoLegend()

ggsave("UMAP_Clusters.png", plot = umap_plot, 
       width = 6, height = 6, dpi = 300, bg = "white")
ggsave("UMAP_Clusters.pdf", plot = umap_plot, 
       width = 6, height = 6)

#DimPlot(merged_seurat, group.by="group",reduction = "umap", label = T)

# 保存最终处理后的数据
saveRDS(merged_seurat, file = "merged_seurat_res.0.2.rds")
