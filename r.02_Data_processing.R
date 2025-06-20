rm(list = ls())
gc()
setwd("D:/workplace/workplace_git/scRNA-seq")
if (! dir.exists("./02_Data_processing")){
	dir.create("./02_Data_processing")
}
setwd("./02_Data_processing")
options(stringsAsFactors = FALSE)

# 软件包安装 -------------------------------------------------------------------

# 定义需要安装的包
required_packages <- c("tidyverse","Seurat","ggsci","celldex", "SingleR","CellChat",
											 "cowplot", "DoubletFinder")
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

# 数据加载与整合 -----------------------------------------------------------------
data_dir <- "../00_rawdata/GSE163973"

# 获取所有样本目录
sample_dirs <- list.files(data_dir, pattern = "^GSM", full.names = TRUE)

# 创建一个空的列表来存储所有样本的Seurat对象
seurat_list <- list()

# 循环读取每个样本的数据
for (i in seq_along(sample_dirs)) {
	sample_name <- basename(sample_dirs[i])
	cat("Processing sample:", sample_name, "\n")
	
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
	
	# 计算线粒体基因比例
	# 线粒体基因比例高可能指示细胞凋亡或低质量
	seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, 
	                                                   pattern = "^MT-") # 匹配所有线粒体基因
	
	# 添加到列表
	seurat_list[[sample_name]] <- seurat_obj
}

# 合并所有样本
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
											 add.cell.ids = names(seurat_list))

# 质量控制与数据预处理 --------------------------------------------------------------
# 质量控制可视化
qc_plot <- VlnPlot(merged_seurat, 
									 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
									 ncol = 3, 
									 pt.size = 0.1,
									 cols = ggsci::pal_npg()(6)
) +
	theme(plot.title = element_text(size=10))

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
			plot.margin = margin(0, .5, 0, .5, "cm")
		)


# 保存为PNG（透明背景可选）
ggsave("qc_plot.png", qc_plot, 
			 width = 12, 
			 height = 6, 
			 dpi = 300, bg = "white")

# 保存为PDF（矢量图，适合论文）
ggsave("qc_plot.pdf", qc_plot, 
			 width = 12, 
			 height = 6, 
			 device = pdf)

# 过滤低质量细胞
# 线粒体基因的百分比 5%
summary(merged_seurat@meta.data[["percent.mt"]])
# 去除 UMI 计数高和低 （>6000 和 <200） 的细胞
filtered_seurat <- subset(merged_seurat, 
													subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
														percent.mt < 5)

# 标准化数据
filtered_seurat <- NormalizeData(filtered_seurat, 
																 normalization.method = "LogNormalize", 
																 scale.factor = 10000)

# 识别高变基因
filtered_seurat <- FindVariableFeatures(filtered_seurat, 
																				selection.method = "vst", 
																				nfeatures = 2000)

# 标准化并缩放数据（Z-score标准化）
filtered_seurat <- ScaleData(filtered_seurat)

# 运行PCA（必须步骤，DoubletFinder依赖PCA结果）
filtered_seurat <- RunPCA(filtered_seurat, 
													features = VariableFeatures(object = filtered_seurat),
													npcs = 50,
													verbose = FALSE)



# 降维与聚类分析 -----------------------------------------------------------------
# 细胞异质性：通过聚类分析识别不同的细胞亚群，特别是成纤维细胞的亚群划分
# PCA分析
filtered_seurat <- RunPCA(filtered_seurat, 
                          features = VariableFeatures(object = filtered_seurat))

# 可视化PCA结果
pca_plot <- DimPlot(filtered_seurat, 
                    reduction = "pca", 
                    group.by = "sample") +
	theme_minimal() +
	ggtitle("PCA by Sample")


# 双峰检测 (DoubletFinder)
sweep.res <- paramSweep(filtered_seurat, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 确定主成分数量
ElbowPlot(filtered_seurat, ndims = 50)

# 选择前20个PCs进行后续分析
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:20)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:20)

# UMAP可视化
umap_plot <- DimPlot(filtered_seurat, 
                     reduction = "umap", 
                     label = TRUE, 
                     pt.size = 0.5) +
	theme_minimal() +
	ggtitle("UMAP Clustering")

# 保存降维图
ggsave("UMAP_clusters.png", umap_plot, width = 8, height = 6, dpi = 300)

# 细胞类型注释 ------------------------------------------------------------------
# 使用SingleR进行自动注释

# 下载参考数据集
# 设置清华镜像（适用于中国用户）
options(repos = c(
  CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
))

#ref <- celldex::HumanPrimaryCellAtlasData()
#saveRDS(ref, file = "ref.rds")
readRDS(file = "ref.rds")

# 查看数据源的预期路径
# 提取表达矩阵
expr <- LayerData(filtered_seurat, assay = "RNA", layer = "data")

# 进行注释
# 将聚类结果映射到已知细胞类型
annotations <- SingleR(test = expr, ref = ref, labels = ref$label.main)

# 将注释结果添加到Seurat对象
filtered_seurat$celltype <- annotations$labels

# 可视化注释结果
annotated_umap <- DimPlot(filtered_seurat, 
                          group.by = "celltype", 
                          label = TRUE, 
                          repel = TRUE) +
	scale_color_manual(values = pal_d3("category20")(20)) +
	theme_minimal() +
	ggtitle("Cell Type Annotation")

# 保存注释结果
ggsave("UMAP_annotated.png", annotated_umap, width = 10, height = 8, dpi = 300)

# 细胞类型占比柱状图
celltype_prop <- as.data.frame(table(filtered_seurat$celltype))
colnames(celltype_prop) <- c("CellType", "Count")

prop_plot <- ggplot(celltype_prop, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
	geom_bar(stat = "identity") +
	scale_fill_manual(values = pal_d3("category20")(nrow(celltype_prop))) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Cell Type", y = "Number of Cells", title = "Cell Type Composition")

ggsave("Celltype_proportion.png", prop_plot, width = 10, height = 6, dpi = 300,bg = "white")

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
