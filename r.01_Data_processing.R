rm(list = ls())
gc()
setwd("D:/workplace/workplace_git/scRNA-seq/")
if (! dir.exists("./01_Data_processing")){
	dir.create("./01_Data_processing")
}
setwd("./01_Data_processing")
options(stringsAsFactors = FALSE)

# 软件包安装 -------------------------------------------------------------------

# 定义需要安装的包
required_packages <- c("tidyverse","Seurat","ggsci","celldex", "SingleR","CellChat",
											 "cowplot", "DoubletFinder","doParallel", "Matrix", "harmony",
											 "GEOquery","R.utils")
# 检查哪些包未安装
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# 如果有缺失的包，安装它们
if (length(missing_packages) > 0) {
	message("正在通过 Bioconductor 安装以下包: ", paste(missing_packages, collapse = ", "))
	
	# 设置国内CRAN镜像（BiocManager会依赖CRAN安装依赖包）
	options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
	
	# 检查是否已安装 BiocManager，若未安装则先安装
	if (!requireNamespace("BiocManager", quietly = TRUE)) {
		message("正在安装 BiocManager...")
		install.packages("BiocManager")
	}
	
	# 使用 BiocManager 安装所有缺失包（即使部分包在CRAN上）
	BiocManager::install(missing_packages, ask = FALSE)  # `ask = FALSE` 避免交互确认
}

# 加载所有包
lapply(required_packages, function(pkg) {
	if (!require(pkg, character.only = TRUE)) {
		warning("无法加载包: ", pkg)
	} else {
		message("成功加载包: ", pkg)
	}
})

# 样本信息整理 ------------------------------------------------------------------
gse <- getGEO("GSE163973", GSEMatrix = TRUE)
sample_info <- pData(phenoData(gse[[1]]))

# 提取关键信息
sample_info_clean <- sample_info %>%
	dplyr::select(title, characteristics_ch1.2) %>%
	mutate(
		sample_id = rownames(sample_info) ,# gsub(".*(GSM[0-9]+).*", "\\1", title),
		condition = ifelse(grepl("sample type: keloid", characteristics_ch1.2), "Keloid", "Normal"),
		#patient_id = gsub(".*patient ([0-9]+).*", "\\1", characteristics_ch1)
	)

# 数据加载与整合 -----------------------------------------------------------------
data_dir <- "../00_rawdata/GSE163973"

# 获取所有样本目录
sample_dirs <- list.files(data_dir, pattern = "^GSM", full.names = TRUE)

sample_info_clean$sample_id
basename(sample_dirs)

# 创建一个空的列表来存储所有样本的Seurat对象
# 初始化列表存储Seurat对象
seurat_list <- list()

# 遍历所有样本
for (i in 1:nrow(sample_info_clean)) {
	sample_path <- file.path("./raw_data", sample_info_clean$file_name[i])
	sample_name <- sample_info_clean$sample_id[i]
	sample_group <- sample_info_clean$condition[i]
	
	cat("Processing:", sample_name, "(", sample_group, ")...\n")
	
	# 读取数据（根据你的数据格式选择合适的方法）
	data <- tryCatch(
		{
			# 如果是10X格式数据
			# Read10X(data.dir = sample_path)
			
			# 如果是CSV格式数据
			read.csv(sample_path, row.names = 1)
		},
		error = function(e) {
			cat("Error reading", sample_name, ":", e$message, "\n")
			NULL
		}
	)
	
	if (is.null(data)) next
	
	# 创建Seurat对象（已满足条件1：min.cells = 3）
	seurat_obj <- CreateSeuratObject(
		counts = data,
		project = sample_name,
		min.cells = 3,      # 基因至少在3个细胞中表达（条件1）
		min.features = 200  # 初始过滤：细胞至少检测到200个基因
	)
	
	# 添加样本信息
	seurat_obj$sample <- sample_name
	seurat_obj$group <- sample_group
	seurat_obj$patient <- sample_info_clean$patient_id[i]
	
	# 计算质控指标
	seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")  # 人类样本
	# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")  # 小鼠样本
	
	# 计算核糖体基因比例（可选）
	seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
	
	# 应用所有过滤条件
	seurat_obj <- subset(seurat_obj,
											 subset = nFeature_RNA > 200 &    # 条件2下限
											 	nFeature_RNA < 4000 &   # 条件2上限
											 	percent.mt < 15 &       # 条件3
											 	nCount_RNA < 10000)     # 条件4
	
	# 添加到列表
	seurat_list[[sample_name]] <- seurat_obj
	
	# 打印过滤后信息
	cat("After filtering:", ncol(seurat_obj), "cells retained in", sample_name, "\n")
}

# 合并所有样本（如果需要）
merged_seurat <- merge(seurat_list[[1]], 
											 y = seurat_list[-1],
											 add.cell.ids = names(seurat_list),
											 project = "merged_project")

seurat_list <- list()

# 循环读取每个样本的数据
for (i in seq_along(sample_dirs)) {
	sample_name <- basename(sample_dirs[i])
	cat("Processing sample:", sample_name, "\n")
	#group_name <- sample_info_clean$group[i]

	# 读取10X Genomics单细胞测序数据
	# 默认输出：稀疏矩阵格式的基因表达矩阵
	data <- Read10X(data.dir = sample_dirs[i])

	# 创建Seurat对象
	# 初步质量控制，去除低质量细胞和噪声基因
	seurat_obj <- CreateSeuratObject(counts = data,
																	 project = sample_name,
																	 min.cells = 3, # 基因至少在3个细胞中表达才会被保留
																	 min.features = 200) # 细胞至少检测到200个基因才会被保留

	# 添加样本信息
	seurat_obj$sample <- sample_name
	#seurat_obj$group <- group_name
	
	# 计算线粒体基因比例
	# 线粒体基因比例高可能指示细胞凋亡或低质量
	seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
																										 pattern = "^MT-") # 匹配所有线粒体基因

	# 添加到列表
	seurat_list[[sample_name]] <- seurat_obj
}

# 添加检查步骤 - 查看基因名称中是否包含线粒体基因
# mt_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
# if(length(mt_genes) == 0) {
# 	mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE)
# }
# if(length(mt_genes) > 0) {
# 	cat("检测到的线粒体基因:", paste(head(mt_genes), collapse = ", "), "...\n")
# } else {
# 	warning("未检测到线粒体基因，请检查基因命名或数据质量")
# }

# 在合并多个样本时（merge()），v5 会保留每个样本的原始数据层（counts/data），而不是自动合并成一个统一层。
# 这是 v5 的新特性，旨在保留样本特异性信息，但会导致 DoubletFinder 等工具报错。
merged_seurat <- merge(seurat_list[[1]],
											 y = seurat_list[-1],
											 add.cell.ids = names(seurat_list))

# 去除文件名称
new_colnames <- str_extract(colnames(merged_seurat), "[ACGT]{16}-\\d+$")
new_colnames <- make.unique(new_colnames)
colnames(merged_seurat) <- new_colnames
head(colnames(merged_seurat))
#saveRDS(merged_seurat, file = "GSE202379_merged_raw.rds")

# untar("../00_rawdata/GSE163973/GSE163973_RAW.tar", exdir = "raw_data")
# 
# # 解压单个gz文件（需要遍历所有文件）
# gz_files <- list.files("raw_data", pattern = "\\.gz$", full.names = TRUE)
# sapply(gz_files, gunzip, overwrite = TRUE)
# 
# tar_files <- list.files("raw_data", pattern = "matrix.tar", recursive = TRUE, full.names = TRUE)
# sapply(tar_files, untar, exdir = "raw_data")


# 数据合并与初步QC ---------------------------------------------------------------
# 初步QC统计
qc_stats_raw <- data.frame(
	Sample = names(seurat_list),
	NASH = sapply(seurat_list, function(x) sum(x$group == "NASH")),
	Normal = sapply(seurat_list, function(x) sum(x$group == "Normal")),
	Total_Cells = sapply(seurat_list, ncol),
	Median_Genes = sapply(seurat_list, function(x) median(x$nFeature_RNA)),
	Median_UMIs = sapply(seurat_list, function(x) median(x$nCount_RNA)),
	Median_MT = sapply(seurat_list, function(x) median(x$percent.mt))
)

# 质量控制与数据预处理 --------------------------------------------------------------
# 生成QC图
qc_plot <- VlnPlot(
	merged_seurat,
	features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
	ncol = 3,
	pt.size = 0.1,
	cols = ggsci::pal_npg()(6)
) &  # 注意这里用 & 而不是 +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
		plot.margin = margin(.5, .5, 0.5, .5, "cm")
	)
qc_plot

# 保存为PNG（透明背景可选）
# ggsave("qc_plot.png", qc_plot, 
# 			 width = 12, 
# 			 height = 6, 
# 			 dpi = 300, bg = "white")

# 保存为PDF（矢量图，适合论文）
# ggsave("qc_plot.pdf", 
# 			 qc_plot, 
# 			 width = 12, 
# 			 height = 6, 
# 			 device = pdf)

# 过滤低质量细胞
# summary(merged_seurat@meta.data[["percent.mt"]])
# FeatureScatter(merged_seurat, 
# 							 feature1 = "percent.mt", 
# 							 feature2 = "nFeature_RNA") +
# 	geom_vline(xintercept = 3, color = "red")
# filtered_seurat <- subset(merged_seurat, 
# 													subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
# 														percent.mt < 3)
# length(colnames(filtered_seurat)) # 9558

library(mclust)
model <- Mclust(merged_seurat$percent.mt, G = 2)  # 假设两类细胞群
max(model$parameters$mean) - sd(model$parameters$mean)

# 去除 UMI 计数高和低 （>6000 和 <200） 的细胞，线粒体基因的百分比 5%
# 文章中是40655 cells
filtered_seurat <- subset(merged_seurat, 
														 subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
														 	percent.mt < 11.58099)
length(colnames(filtered_seurat)) # 43224

# 标准化数据
filtered_seurat <- NormalizeData(filtered_seurat,
																 normalization.method = "LogNormalize",
																 scale.factor = 10000)

# 识别高变基因
filtered_seurat <- FindVariableFeatures(filtered_seurat,
																				selection.method = "vst",
																				nfeatures = 2000)

# 数据集中表现出高细胞间差异的特征子集
top10 <- head(VariableFeatures(filtered_seurat), 10)
plot1 <- VariableFeaturePlot(filtered_seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 标准化并缩放数据（Z-score标准化）
filtered_seurat <- ScaleData(filtered_seurat)

#  运行PCA
#  降维和提取主要变异方向
filtered_seurat <- RunPCA(filtered_seurat,
													features = VariableFeatures(object = filtered_seurat),
													npcs = 50,
													verbose = FALSE)

VizDimLoadings(filtered_seurat, dims = 1:2, reduction = "pca")
DimPlot(filtered_seurat,reduction = "pca") + NoLegend()
DimHeatmap(filtered_seurat, 
					 dims = 1:15, 
					 cells = 500, 
					 balanced = TRUE)

#  批次效应校正（batch correction）
filtered_seurat  <- RunHarmony(filtered_seurat ,
															 group.by.vars = "sample",
															 max.iter.harmony = 20)

filtered_seurat <- RunUMAP(filtered_seurat, 
													 reduction = "harmony", 
													 dims = 1:15)

filtered_seurat <- FindNeighbors(filtered_seurat, 
																 reduction = "harmony", 
																 dims = 1:15)

filtered_seurat <- FindClusters(filtered_seurat, 
																resolution = 0.5) # Louvain算法

Layers(filtered_seurat)

# 双峰检验 (DoubletFinder) ----------------------------------------------------
# 双峰细胞（Doublets）是指在单细胞RNA测序（scRNA-seq）实验中，两个或多个细胞被错误地捕获在同一个液滴或孔中，导致测序数据混合了多个细胞的基因表达谱。
# 这种技术假象会影响下游分析，因此需要通过生物信息学方法（如DoubletFinder）进行检测和去除。
# 虚假细胞类型：双峰的基因表达谱是两种细胞的混合，可能被误认为是新的细胞亚群或过渡状态。
# 干扰差异表达分析：双峰的混合信号会掩盖真实的细胞类型特异性基因表达。
# 误导细胞轨迹推断：在拟时序分析中，双峰可能被错误解释为中间状态或分化过渡细胞。

# 参数扫描（pK优化）
# 合并所有counts层
filtered_seurat_join <- JoinLayers(filtered_seurat)
Layers(filtered_seurat_join)

sweep.res <- paramSweep(filtered_seurat_join,
												PCs = 1:15,
												sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 计算每个样本的理论双峰率
sample_counts <- as.data.frame(table(filtered_seurat_join$sample))
colnames(sample_counts) <- c("Sample", "CellCount")

# 10x官方公式：双峰率 = 细胞数/1000 × 0.008
sample_counts$DoubletRate <- round(sample_counts$CellCount / 1000 * 0.008, 4)
sample_counts$ExpectedDoublets <- round(sample_counts$CellCount * sample_counts$DoubletRate)
print(sample_counts)

conservative_rate <- mean(sample_counts$DoubletRate)
nExp_global <- round(ncol(filtered_seurat_join) * conservative_rate)
cat(paste0("保守估计双峰数",mean(sample_counts$DoubletRate),"："), nExp_global, "\n")

# 运行双峰检测
filtered_seurat_join <- doubletFinder(
	filtered_seurat_join,
	PCs = 1:15,
	pN = 0.25,
	pK = pK,
	nExp = nExp_global
)

# 提取双峰预测结果
doublet_col <- grep("DF.classifications", names(filtered_seurat_join@meta.data), value = TRUE)
table(filtered_seurat_join@meta.data[, doublet_col])

# 双峰在UMAP上的分布
DimPlot(filtered_seurat_join, 
				group.by = doublet_col,
				order = c("Doublet"),  # 将双峰显示在最上层
				cols = c("gray", "red")) +
	ggtitle("Predicted Doublets")

# 双峰的QC指标特征
VlnPlot(filtered_seurat_join,
				features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
				group.by = doublet_col)

observed_doublets <- sum(filtered_seurat_join@meta.data[, doublet_col] == "Doublet")
cat("理论双峰数:", sum(sample_counts$ExpectedDoublets), 
		"\n实际预测数:", observed_doublets)

# 使用scDblFinder生成合成双峰：
# library(scDblFinder)
# simulated_doublets <- scDblFinder(filtered_seurat@assays$RNA@counts)
# table(simulated_doublets$scDblFinder.class)
# filtered_seurat$scDblFinder <- simulated_doublets$scDblFinder.class

# 提取非双峰细胞
singlets <- colnames(filtered_seurat_join)[filtered_seurat_join@meta.data[, doublet_col] == "Singlet"]
filtered_seurat_clean <- subset(filtered_seurat_join, cells = singlets)

# 验证过滤效果
dim(filtered_seurat_join)  # 原始维度
dim(filtered_seurat_clean) # 过滤后维度

# 标准化数据
filtered_seurat_clean <- NormalizeData(filtered_seurat_clean, 
																 normalization.method = "LogNormalize", 
																 scale.factor = 10000)

# 识别高变基因
filtered_seurat_clean <- FindVariableFeatures(filtered_seurat_clean, 
																				selection.method = "vst", 
																				nfeatures = 2000)

# 标准化并缩放数据（Z-score标准化）
filtered_seurat_clean <- ScaleData(filtered_seurat_clean)

#  运行PCA
#  降维和提取主要变异方向
filtered_seurat_clean <- RunPCA(filtered_seurat_clean, 
													features = VariableFeatures(object = filtered_seurat_clean),
													npcs = 50,
													verbose = FALSE)

#  批次效应校正（batch correction）
filtered_seurat_clean  <- harmony::RunHarmony(filtered_seurat_clean ,
																				group.by.vars = "sample",
																				max.iter.harmony = 20)

#  Uniform Manifold Approximation and Projection (UMAP)是一种基于流形学习的非线性降维技术
filtered_seurat_clean <- RunUMAP(filtered_seurat_clean, 
																 reduction = "harmony", 
																 dims = 1:15)

filtered_seurat_clean <- FindNeighbors(filtered_seurat_clean, 
																 reduction = "harmony", 
																 dims = 1:15)

filtered_seurat_clean <- FindClusters(filtered_seurat_clean, 
																resolution = 0.5) # Louvain算法
Layers(filtered_seurat_clean)

# 重新命名分类
current_clusters <- levels(Idents(filtered_seurat_clean))
new_labels <- paste0("C-", seq_along(current_clusters))
names(new_labels) <- current_clusters
filtered_seurat_clean <- RenameIdents(filtered_seurat_clean, new_labels)

# UMAP可视化
# DimPlot(filtered_seurat_clean,
# 				reduction = "umap",
# 				label = TRUE,
# 				pt.size = 0.5) +
# 	ggtitle("Unsupervised clustering of dermal cells") +
# 	theme(plot.title = element_text(hjust = 0.5))

# 使用SingleR进行自动注释
# 下载参考数据集
# 设置清华镜像（适用于中国用户）
options(repos = c(
	CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
	Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
))

ref <- celldex::HumanPrimaryCellAtlasData()
#saveRDS(ref, file = "ref.rds")
#readRDS(file = "ref.rds")

# 查看数据源的预期路径
# 提取表达矩阵
expr <- LayerData(filtered_seurat_clean, assay = "RNA", layer = "data")

# 进行注释
# 将聚类结果映射到已知细胞类型
annotations <- SingleR(test = expr, 
											 ref = ref, 
											 labels = ref$label.main)

# 获取唯一标签
unique_labels <- unique(annotations$labels)

# 为每个标签添加前缀
prefixed_labels <- paste0("C", 1:length(unique_labels), "-", unique_labels)

# 如果你想将这些带前缀的标签重新赋给原始数据
annotations$labels <- factor(annotations$labels, 
														 levels = unique_labels,
														 labels = prefixed_labels)

filtered_seurat_clean$celltype <- annotations$labels

# 可视化注释结果
annotated_umap <- DimPlot(filtered_seurat_clean, 
													group.by = "celltype", 
													label = TRUE, 
													repel = TRUE) +
	scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(25)) +
  ggtitle("Cell Type Annotation") +
  theme_minimal() +
  guides(color = guide_legend(  # 控制图例
    ncol = 1,                  # 单列显示
    override.aes = list(size = 5)  # 图例点大小
  )) +
  theme(
    legend.position = "right",  # 图例在右侧
    legend.text = element_text(size = 12),  # 图例文字大小
    plot.title = element_text(hjust = 0.5, face = "bold")  # 标题居中加粗
  )

# 保存注释结果
ggsave("Fig_B_UMAP_annotated.png", 
			 annotated_umap, 
			 width = 10, 
			 height = 8, 
			 dpi = 300)

ggsave("Fig_B_UMAP_annotated.pdf", 
			 annotated_umap, 
			 width = 10, 
			 height = 8, 
			 dpi = 300,
			 device = pdf)

# 细胞类型占比柱状图
celltype_prop <- as.data.frame(table(filtered_seurat_clean$celltype))
colnames(celltype_prop) <- c("CellType", "Count")

celltype_prop <- as.data.frame(table(filtered_seurat_clean$celltype, 
																		 filtered_seurat_clean$sample)) # Replace 'sample' with your actual grouping variable
colnames(celltype_prop) <- c("CellType", "Group", "Count")

# 简化样本名称
# 删除 "GSM" 开头的所有数字 + 下划线，以及 "_matrix" 后缀
celltype_prop$Group <- gsub("^GSM[0-9]+_|_matrix$", "", celltype_prop$Group)

prop_plot_proportion <- ggplot(celltype_prop, 
															 aes(x = Group, 
															 		y = Count, 
															 		fill = CellType)) +
	geom_bar(stat = "identity", 
					 position = "fill") +
	scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(25)) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, 
																	 hjust = 1, 
																	 size = 10)) +
	labs(x = "Group", 
			 y = "Proportion", 
			 title = "Cell Type Proportions by Group") +
	guides(fill = guide_legend(ncol = 1, 
														 override.aes = list(size = 3)))

ggsave("Fig_F_Celltype_proportion.png", prop_plot_proportion, width = 6, height = 8, dpi = 300,bg = "white")
ggsave("Fig_F_Celltype_proportion.pdf", prop_plot_proportion, width = 6, height = 8, dpi = 300,bg = "white",
			 device = pdf)

# 聚类生物标志物 -----------------------------------------------------------------
filtered_seurat_clean$sample <- gsub("^GSM[0-9]+_|_matrix$", "", filtered_seurat_clean$sample)
unique(filtered_seurat_clean$sample)

Idents(filtered_seurat_clean) <- "celltype"

# 查看所有细胞类型
celltypes <- unique(filtered_seurat_clean$celltype)
print(celltypes)

# 创建一个空列表存储结果
deg_results <- list()

# 循环计算每个细胞类型的DEGs
for (ct in celltypes) {
	cat("\nProcessing celltype:", ct, "\n")
	
	# 提取当前细胞类型的子集
	subset_data <- subset(filtered_seurat_clean, celltype == ct)
	
	# 设置Idents为group（Disease vs Control）
	Idents(subset_data) <- "group"
	
	# 计算DEGs（Disease vs Control）
	degs <- FindMarkers(
		subset_data,
		ident.1 = "Disease",   # KL组
		ident.2 = "Control",   # NS组
		only.pos = FALSE,      # 返回上下调基因
		logfc.threshold = 0.25, # 最小log2FC阈值
		min.pct = 0.1,         # 基因在至少10%的细胞中表达
		test.use = "wilcox"    # 使用Wilcoxon秩和检验
	)
	
	# 添加细胞类型信息
	degs$celltype <- ct
	degs$gene <- rownames(degs)
	
	# 存储结果
	deg_results[[ct]] <- degs
}

# 合并所有结果
all_degs <- do.call(rbind, deg_results)
rownames(all_degs) <- NULL

# 保存结果
write.csv(all_degs, "DEGs_per_celltype_Disease_vs_Control.csv", row.names = FALSE)

# 筛选显著DEGs
sig_degs <- all_degs %>%
	filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

# 按细胞类型和log2FC排序
sig_degs <- sig_degs %>%
	arrange(celltype, desc(avg_log2FC))

# 查看每个细胞类型的DEG数量
deg_counts <- sig_degs %>%
	group_by(celltype) %>%
	summarise(n_up = sum(avg_log2FC > 0),   # 上调基因数
						n_down = sum(avg_log2FC < 0)) # 下调基因数

print(deg_counts)

# 为每个细胞类型绘制火山图
for (ct in celltypes) {
	# 提取当前细胞类型的DEGs
	ct_degs <- all_degs %>% filter(celltype == ct)
	
	# 标记显著基因
	ct_degs$significant <- ifelse(
		ct_degs$p_val_adj < 0.05 & abs(ct_degs$avg_log2FC) > 0.5,
		"Significant", "Not significant"
	)
	
	# 绘制火山图
	volcano_plot <- ggplot(ct_degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
		geom_point(aes(color = significant), alpha = 0.6) +
		scale_color_manual(values = c("gray", "red")) +
		geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
		geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
		labs(
			title = paste("DEGs in", ct, "(Disease vs Control)"),
			x = "log2(Fold Change)",
			y = "-log10(Adjusted p-value)"
		) +
		theme_minimal()
	
	# 保存图像
	ggsave(
		paste0("Volcano_Plot_", gsub("/", "_", ct), ".png"),
		volcano_plot,
		width = 8,
		height = 6,
		dpi = 300
	)
}

# 提取每个细胞类型的前5个上调和下调DEGs
top_degs <- sig_degs %>%
	group_by(celltype, direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
	top_n(5, abs(avg_log2FC))

# 提取表达矩阵
expr_matrix <- LayerData(filtered_seurat_clean, assay = "RNA", layer = "data")

# 绘制热图
heatmap_plot <- DoHeatmap(
	filtered_seurat_clean,
	features = unique(top_degs$gene),
	group.by = "celltype",
	slot = "scale.data",
	disp.min = -2,
	disp.max = 2
) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red")

# 保存热图
ggsave(
	"Heatmap_Top_DEGs_per_Celltype.png",
	heatmap_plot,
	width = 12,
	height = 10,
	dpi = 300
)



unique(all_degs$celltype)

Fibroblasts <-top_degs %>% 
	filter(celltype == "C5-Fibroblasts")

Fibroblasts$gene

#Fig_e_expression_distribution
FeaturePlot(filtered_seurat_clean, features = Fibroblasts$gene)


#Fig_d_differentially_expressed_genes
DimHeatmap(filtered_seurat_clean, 
					 dims = 1:25, 
					 cells = 1000, 
					 balanced = TRUE)

# 从markers_group中提取显著差异表达基因
markers_group <- markers_group[markers_group$p_val_adj < 0.05, ]  # 筛选显著基因

# 选取上下调基因TOP10
top_up <- markers_group %>% 
	filter(avg_log2FC > 0) %>% 
	arrange(p_val_adj, desc(avg_log2FC)) %>% 
	head(10)

top_down <- markers_group %>% 
	filter(avg_log2FC < 0) %>% 
	arrange(p_val_adj, avg_log2FC) %>% 
	head(10)

top_genes <- rbind(top_up, top_down)

# 提取表达矩阵
expr_matrix <- as.matrix(GetAssayData(filtered_seurat_clean, slot = "scale.data"))

VlnPlot(filtered_seurat_clean, features = rownames(top_genes))




# 成纤维细胞亚群分析 ---------------------------------------------------------------
# 提取成纤维细胞亚群
# 深入研究疾病相关细胞类型的异质性
fibroblast_subset <- subset(filtered_seurat, subset = celltype == "Fibroblasts")

# 重新处理成纤维细胞数据
fibroblast_subset <- NormalizeData(fibroblast_subset,  
																	 scale.factor=10000) # LogNormalize 自然对数转换
# 聚焦于在细胞间差异表达的基因
fibroblast_subset <- FindVariableFeatures(fibroblast_subset,
																					nfeatures=2000) # "vst"（方差稳定变换）选择变异度最高的2000个基因 

fibroblast_subset <- ScaleData(fibroblast_subset)

# 使用高变基因进行主成分分析 降维同时保留数据主要变异来源
fibroblast_subset <- RunPCA(fibroblast_subset)

# 使用前15个PCs 基于KNN（k-nearest neighbors）构建细胞相似性图
fibroblast_subset <- FindNeighbors(fibroblast_subset, dims = 1:15)

# 控制聚类粒度值越大，聚类数越多  Louvain算法（基于模块度优化）
fibroblast_subset <- FindClusters(fibroblast_subset, resolution = 0.3)

# 非线性降维可视化细胞关系
fibroblast_subset <- RunUMAP(fibroblast_subset, 
														 dims = 1:15,
														 n.neighbors=30,
														 min.dist=0.3)

# 成纤维细胞UMAP图
fib_umap <- DimPlot(fibroblast_subset,
										reduction = "umap", 
										label = TRUE) +
	theme_minimal() +
	ggtitle("Fibroblast Subclustering")

# 成纤维细胞marker基因小提琴图
fib_markers <- c("COL1A1", "COL3A1", "ACTA2", "PDGFRA", "PDGFRB")
fib_vln <- VlnPlot(fibroblast_subset, features = fib_markers, ncol = 3, pt.size = 0) &
	theme(axis.title.x = element_blank()) &
	scale_fill_d3()

# 保存成纤维细胞分析结果
ggsave("Fibroblast_UMAP.png", fib_umap, width = 8, height = 6, dpi = 300)
ggsave("Fibroblast_markers.png", fib_vln, width = 12, height = 8, dpi = 300)

# 分样本UMAP图
fib_sample_umap <- DimPlot(fibroblast_subset, group.by = "sample") +
	theme_minimal() +
	ggtitle("Fibroblast Subclustering by Sample")

ggsave("Fibroblast_sample_UMAP.png", fib_sample_umap, width = 8, height = 6, dpi = 300)

# 细胞通讯 --------------------------------------------------------------------
# 准备CellChat分析

# 提取细胞类型信息
cell_metadata <- filtered_seurat@meta.data
cell_metadata$celltype

cell_metadata <- cell_metadata[rownames(cell_metadata) %in% expr@Dimnames[[2]],]

# 创建 CellChat 对象
cellchat <- createCellChat(object = expr, 
													 meta = cell_metadata, 
													 group.by = "celltype")

# 设置配体-受体数据库
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 预处理细胞通讯数据
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 计算通讯概率
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 计算整合的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 可视化细胞通讯网络
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
								 weight.scale = TRUE, label.edge = FALSE,
								 title.name = "Number of interactions")

# 保存细胞通讯图
png("./Cell_communication_network.png", width = 800, height = 800)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
								 weight.scale = TRUE, label.edge = FALSE,
								 title.name = "Number of interactions")
dev.off()

# 热图展示细胞通讯数量
# 揭示不同细胞类型间的信号交流模式
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")

png("./Cell_communication_heatmap.png", width = 800, height = 800)
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")
dev.off()

# 其他可视化 -------------------------------------------------------------------
# Top10基因聚类热图
top10 <- FindAllMarkers(
	filtered_seurat,
	only.pos = TRUE,
	min.pct = 0.05,      # 降低最小表达比例
	logfc.threshold = 0.01  # 降低logFC阈值
)

# 如果此时有结果，说明原阈值过高
markers <- FindAllMarkers(
	filtered_seurat,
	only.pos = FALSE,  # 包含下调基因
	min.pct = 0,       # 关闭min.pct过滤
	logfc.threshold = 0 # 关闭logFC过滤
)


heatmap_plot <- DoHeatmap(filtered_seurat, features = top10$gene) + 
	scale_fill_viridis() +
	theme(axis.text.y = element_text(size = 6))

ggsave("./Top10_genes_heatmap.png", heatmap_plot, width = 12, height = 16, dpi = 300)

# Marker基因FeaturePlot
marker_genes <- c("EPCAM", "PECAM1", "PTPRC", "COL1A1", "ACTA2")
feature_plots <- FeaturePlot(filtered_seurat, features = marker_genes, 
														 ncol = 3, pt.size = 0.1, combine = FALSE)
feature_plots <- lapply(feature_plots, function(x) x + theme_minimal() + viridis::scale_color_viridis())
combined_feature <- patchwork::wrap_plots(feature_plots, ncol = 3)

ggsave("./Marker_genes_featureplot.png", combined_feature, width = 15, height = 10, dpi = 300)
