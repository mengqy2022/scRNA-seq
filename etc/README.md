[TOC]

纤维化皮肤病是全球医疗保健的主要负担，其特征是成纤维细胞过度增殖和细胞外基质过度积累。发现成纤维细胞在多种纤维化疾病中具有异质性，但在纤维化皮肤病中成纤维细胞异质性尚未得到很好的表征。在这项研究中，我们通过使用单细胞 RNA-seq 探索了瘢痕疙瘩中的成纤维细胞异质性，这是纤维化皮肤病的一种范例。我们的结果表明，瘢痕疙瘩成纤维细胞可分为 4 个亚群：分泌-状、分泌-网状、间充质和促炎性。有趣的是，与正常疤痕相比，瘢痕疙瘩中间充质成纤维细胞亚群的百分比显着增加。功能研究表明，间充质成纤维细胞对于瘢痕疙瘩中胶原蛋白的过表达至关重要。在另一种纤维化皮肤病硬皮病中也发现了间充质成纤维细胞亚群的增加，这表明这是皮肤纤维化的广泛机制。这些发现将帮助我们更好地了解皮肤纤维化的发病机制，并为纤维化疾病治疗提供潜在的靶点。

样品制备和组织解离。本研究经南方医科大学 （2019023） 皮肤科医院医学和伦理委员会批准，每位患者在参加本研究前签署知情同意书。本研究中的所有患者均为汉族。在整形手术期间，从三名确认有瘢痕疙瘩临床证据的患者身上采集了瘢痕疙瘩组织（补充表 1）。我们在这项研究中使用的所有瘢痕疙瘩都是成熟的。我们使用了瘢痕疙瘩样品的所有内容物，包括样品的中心和边缘，并将它们混合以进行进一步分析。没有患者在手术前接受化疗、放疗或病灶内类固醇治疗。3 例接受择期瘢痕切除手术的患者获得正常瘢痕组织 （补充表 1）。瘢痕疙瘩和正常瘢痕根据其临床表现、病史、解剖位置和病理学进行诊断。将切除的皮肤浸入生理盐水中，然后立即转移到实验室。皮肤组织在 PBS 中洗涤两次。去除网状真皮下的脂肪组织后，将样品切成直径为 5 mm 的块，并与分散酶 II （Sigma） 在 37 °C 下孵育 2 小时。 剥去表皮并丢弃，将真皮切成小块，并使用胶原酶IV（YEASEN，中国）在37°C下消化2小时。将所得细胞悬液通过 70 μm 细胞过滤器 （BD Falcon） 过滤，并以 400 × g 离心 10 分钟。除去上清液，用 400 × g 的 PBS 洗涤沉淀一次，每次 5 分钟。然后将沉淀重悬于 PBS + 1% FBS 中用于流式细胞术。

单细胞 cDNA、文库制备和 3′ 末端单细胞 RNA 测序由 Novogene （中国北京） 进行。对于使用10×基因组学平台的实验，根据Chromium单细胞3′试剂盒v2用户指南中的制造商说明，使用了Chromium单细胞3′文库和凝胶珠套装v2（PN-120237），Chromium单细胞3′芯片试剂盒v2，和Chromium i7多路复用试剂盒（PN-120262）。单细胞悬液用 1× PBS + 0.04% BSA 洗涤两次。使用 TC20™ 自动细胞计数仪确认细胞数量和浓度。立即将大约 8000 个细胞置于 10× Genomics 铬控制器机器中，以生成乳状凝胶珠 （GEM）。使用 10× Genomics Chromium Single Cell 3′ 试剂盒（V2 化学试剂）制备 mRNA。在此步骤中，将细胞与包被寡核苷酸的凝胶珠一起分配到 GEM 中。这些寡核苷酸提供 poly-dT 序列，用于捕获液滴内细胞裂解后释放的 mRNA，以及细胞特异性和转录本特异性条形码（分别为 16 bp 10× 条形码和 10 bp 唯一分子标识符 （UMI）。RT-PCR 后，回收、纯化和扩增 cDNA，以产生足够的量用于文库制备。使用 Agilent Bioanalyzer 2100 评估文库质量和浓度。

3′ 末端单细胞 RNA 测序。在 Hiseq X 或 Novaseq 上运行文库，用于 Illumina PE150 测序。后处理和质量控制由 Novogene 使用 10× Cell Ranger 软件包 （v2.1.0， 10× Genomics） 进行。读数与 GRCh38 参考组装 （v2.2.0， 10× Genomics 进行比对。

# 结果与方法

由于瘢痕疙瘩是一种皮肤真皮层纤维化疾病，因此我们只使用真皮层进行 scRNA-seq 分析。经过严格的质量控制（补充图 1a、b），我们获得了 40,655 个细胞的转录组（瘢痕疙瘩：21,488 个；正常瘢痕：19,167 个）。

无监督统一表层逼近和投影（UMAP）聚类发现了 21 个细胞簇（图 1b 和补充图 1c）、

基于分层聚类（图 1c）和已建立的谱系特异性标记基因（图 1d、e），我们将这些簇分为 9 个细胞谱系。成纤维细胞谱系由 COL1A1 鉴定，内皮谱系由 ENG 鉴定 （图 1e）。

接下来，我们分析了这些细胞谱系在瘢痕疙瘩和正常疤痕中的比例。我们在比例分析中去除了表皮中的细胞，包括角质形成细胞和黑色素细胞。瘢痕疙瘩和正常瘢痕真皮的细胞谱系显示出不同的相对细胞数比（图 1f）。

由于成纤维细胞在瘢痕疙瘩纤维化过程中发生显着变化（图 1g），并且成纤维细胞对纤维化发病机制很重要，因此我们接下来对所有瘢痕疙瘩和正常瘢痕成纤维细胞进行了无监督聚类，并观察到 13 个亚簇 sC1 到 sC13 的进一步异质性（图 2a）。

图 2b、c 显示了来自瘢痕疙瘩和正常瘢痕的成纤维细胞亚簇的细胞比例。

已知成纤维细胞与皮肤中的其他细胞类型相互作用。scRNA-seq 提供了根据细胞表面受体及其相互作用配体的表达来识别通讯细胞对的机会。结果表明，成纤维细胞与我们在 scRNA-seq 中鉴定的所有细胞相互作用（图 4a）。值得注意的是，在正常疤痕中，分泌网状成纤维细胞和其他细胞之间的相互作用最丰富。然而，在纤维化皮肤中，最丰富的相互作用发生在间充质成纤维细胞和其他细胞之间（图 4a）。这些结果表明，这些间充质成纤维细胞在皮肤纤维化发展中起重要作用，这与它们在皮肤纤维化中的增加一致。我们发现与正常疤痕相比，瘢痕疙瘩中 TGFβ-TGFβ 受体相互作用显着增加（图 4b），表明 TGFβ 通路在纤维化发展中起着核心作用。我们还确定了一些以前报道的纤维化相关相互作用，例如 NOTCH 和血管生成相关的 VEGF 配体-受体相互作用（图 4b）6,38，这些相互作用在皮肤纤维化中增加。有趣的是，我们发现 POSTN-ITGAV;与正常瘢痕相比，瘢痕疙瘩中间充质细胞和其他成纤维细胞之间的 ITGB5 相互作用显著增加（图 4b、c）。



**上游：**通过质量指标检查测序读数，并使用 Cell Ranger 流程 （10× Genomics 将转录本映射到参考人类基因组 （hg38），并根据细胞特异性条形码分配给单个来源细胞。为确保仅对 PCR 扩增的转录本计数一次，仅对单个 UMI 进行基因表达水平分析41。通过这种方式，生成了用于下游分析的细胞基因 UMI 计数基质。

**下游：**通过去除 UMI 计数高和低 （>6000 和 <200） 的细胞，从每个样品中过滤不需要的变异和低质量细胞。同时，为避免双峰的影响，使用 Doubletdetection从我们的数据中删除了被鉴定为双峰的细胞。

通过总表达、乘以比例因子 （10,000） 对每个细胞的基因表达水平进行归一化，并进行对数转换。然后将批次回归出去，并将模型的缩放 Z 分数残差用作标准化表达值。我们根据它们的平均表达和分散性将前 2000 个最具可变性的基因定义为高度可变基因 （HVG）。我们通过对 HVG 进行主成分分析 （PCA） 来降低数据的维度。为了识别细胞亚群，使用 Chung 和 Storey 提出的随机化方法分配的显着 PC 对 PCA 评分进行聚类。对于这些重复（NS1、NS2、NS3 和 KF1、KF2、KF3），选择前 15 台个PCs 进行聚类。为了对细胞进行聚类，计算了在 PCA 空间中的欧几里得距离矩阵上构建的 K 最近邻 （KNN） 图，然后将其转换为共享的最近邻 （SNN） 图，以找到高度互连的细胞群落。然后使用鲁汶方法对细胞进行聚类，以最大限度地提高模块化程度。为了显示数据，将无监督均匀流形近似和投影 （UMAP） 应用于所选 PC 的单元负载，并且来自基于图的聚类的聚类分配是

细胞间通讯分析：为了确定成纤维细胞和其他真皮细胞群之间和内部的潜在相互作用，我们使用了参数阈值 = 0.25 且迭代次数 = 100027的 CellPhoneDB 2.0，其中包含配体-受体相互作用的精选存储库和用于推断谱系特异性相互作用的统计框架。使用自定义 R 脚本和 circos 软件进行分析和绘制交互图。

```R
# 加载必要的R包
library(Seurat)        # 单细胞分析核心工具
library(dplyr)         # 数据操作
library(ggplot2)       # 可视化
library(cowplot)       # 高级绘图
library(CellPhoneDB)   # 细胞通讯分析
library(harmony)       # 批次校正
library(DoubletFinder) # 双峰检测

### 1. 数据预处理和质量控制 ###
# 读取10x Genomics数据
sc_data <- Read10X(data.dir = "path/to/GSE163973/filtered_feature_bc_matrix")

# 创建Seurat对象
sc_obj <- CreateSeuratObject(counts = sc_data, 
                           project = "Keloid", 
                           min.cells = 3, 
                           min.features = 200)

# 添加样本元数据 (需根据实际样本信息调整)
sc_obj$sample <- c(rep("NS1", n1), rep("NS2", n2), rep("NS3", n3),
                   rep("KF1", m1), rep("KF2", m2), rep("KF3", m3))
sc_obj$group <- ifelse(grepl("KF", sc_obj$sample), "Keloid", "Normal")

# 计算线粒体和核糖体基因比例
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj[["percent.rb"]] <- PercentageFeatureSet(sc_obj, pattern = "^RP[SL]")

# 质量控制过滤 (参数根据补充图1a,b调整)
sc_obj <- subset(sc_obj, 
                subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 6000 & 
                         percent.mt < 20)

# 双峰检测 (DoubletFinder)
sweep.res <- paramSweep_v3(sc_obj, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp <- round(ncol(sc_obj) * 0.04)  # 预计双峰率4%
sc_obj <- doubletFinder_v3(sc_obj, PCs = 1:15, pN = 0.25, pK = pK, nExp = nExp)

# 移除双峰
sc_obj <- subset(sc_obj, cells = colnames(sc_obj)[which(sc_obj@meta.data[, grep("DF.class", colnames(sc_obj@meta.data))] == "Singlet")])

### 2. 数据标准化和批次校正 ###
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

# 使用Harmony进行批次校正
sc_obj <- ScaleData(sc_obj, vars.to.regress = c("percent.mt", "percent.rb"))
sc_obj <- RunPCA(sc_obj, npcs = 50, verbose = FALSE)
sc_obj <- RunHarmony(sc_obj, group.by.vars = "sample", max.iter.harmony = 20)

### 3. 细胞聚类和UMAP可视化 ###
sc_obj <- FindNeighbors(sc_obj, reduction = "harmony", dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.8, algorithm = 1) # Louvain算法
sc_obj <- RunUMAP(sc_obj, reduction = "harmony", dims = 1:15)

# 绘制UMAP图 (对应图1b)
DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Unsupervised clustering of dermal cells") +
  theme(plot.title = element_text(hjust = 0.5))

### 4. 细胞类型注释 ###
# 使用已知标记基因鉴定细胞类型 (对应图1d,e)
markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "DCN"),
  Endothelial = c("ENG", "PECAM1", "VWF"),
  Immune = c("PTPRC", "CD3E", "CD79A"),
  Pericytes = c("RGS5", "PDGFRB"),
  Schwann = c("S100B", "SOX10"),
  Melanocytes = c("MLANA", "PMEL")
)

# 标记基因可视化
FeaturePlot(sc_obj, features = unlist(markers), ncol = 4, pt.size = 0.1)

# 分配细胞类型标识
current.cluster.ids <- 0:20
new.cluster.ids <- c("Fibroblast_1", "Fibroblast_2", "Endothelial", ...) # 根据标记基因表达定义
sc_obj$celltype <- plyr::mapvalues(sc_obj$seurat_clusters, 
                                  from = current.cluster.ids, 
                                  to = new.cluster.ids)

### 5. 成纤维细胞亚群分析 ###
# 提取成纤维细胞子集
fibroblasts <- subset(sc_obj, idents = c("Fibroblast_1", "Fibroblast_2", ...))

# 重新聚类 (对应图2a)
fibroblasts <- FindVariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts)
fibroblasts <- RunPCA(fibroblasts, npcs = 20)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.6)
fibroblasts <- RunUMAP(fibroblasts, dims = 1:15)

# 鉴定成纤维细胞亚群标记基因 (sC1-sC13)
fibro_markers <- FindAllMarkers(fibroblasts, only.pos = TRUE, min.pct = 0.25)

### 6. 细胞比例分析 ###
# 计算各组细胞比例 (对应图1f,g)
prop_table <- prop.table(table(sc_obj$celltype, sc_obj$group), margin = 2)
melted_prop <- reshape2::melt(prop_table)

ggplot(melted_prop, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Group", y = "Proportion", fill = "Cell type") +
  theme_minimal()

### 7. 细胞通讯分析 (CellPhoneDB) ###
# 准备CellPhoneDB输入文件
write.table(as.matrix(fibroblasts@assays$RNA@data), "cpdb_counts.txt", sep = "\t", quote = F)
write.table(data.frame(Cell = rownames(fibroblasts@meta.data),
                      cell_type = fibroblasts$celltype),
            "cpdb_meta.txt", sep = "\t", quote = F, row.names = F)

# 在终端运行CellPhoneDB (参数对应文中描述)
# cellphonedb method statistical_analysis cpdb_meta.txt cpdb_counts.txt --iterations=1000 --threshold=0.25

# 可视化相互作用 (对应图4)
# 加载CellPhoneDB输出结果
interactions <- read.delim("out/significant_means.txt", check.names = F)

# TGFβ通路相互作用热图 (图4b)
tgfb_inter <- interactions[grep("TGFB", interactions$interacting_pair), ]
ggplot(tgfb_inter, aes(x = group, y = interacting_pair, fill = means)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

# POSTN-ITGAV/ITGB5相互作用circos图 (图4c)
library(circlize)
postn_inter <- interactions[interactions$interacting_pair %in% c("POSTN_ITGAV_ITGB5"), ]

# 准备circos数据
mat <- matrix(postn_inter$means, nrow = length(unique(postn_inter$group)))
rownames(mat) <- unique(postn_inter$group)
colnames(mat) <- unique(postn_inter$interacting_pair)

# 绘制circos图
chordDiagram(mat, transparency = 0.5)
```

### 关键生物学解释：

1. **质量控制**：
   - 线粒体基因比例过滤去除低质量/凋亡细胞
   - UMI计数范围(200-6000)确保捕获真实细胞
   - 双峰检测提高单细胞分辨率
2. **细胞异质性**：
   - 21个初始聚类反映皮肤真皮复杂组成
   - 成纤维细胞进一步分13个亚群(sC1-sC13)显示功能异质性
3. **纤维化相关发现**：
   - 间充质成纤维细胞比例增加与纤维化正相关
   - TGFβ通路激活是瘢痕疙瘩的核心特征
   - POSTN-ITGAV/ITGB5相互作用可能促进ECM沉积
4. **方法学要点**：
   - Harmony校正解决样本间批次效应
   - CellPhoneDB参数(1000次迭代,0.25阈值)确保互作可靠性
   - 成纤维细胞亚群分辨率(0.6)平衡过/欠聚类

### 注意事项：

1. 实际运行时需要根据GSE163973数据实际结构调整路径和参数
2. CellPhoneDB分析需要在Python环境中运行
3. 可视化细节(颜色、标签等)可根据期刊要求调整
4. 建议使用R 4.0以上版本和相应包版本以确保兼容性

如需更详细的某部分代码实现或遇到具体报错，可以提供更具体的需求或错误信息，我可以进一步优化代码。

~~~R
---
title: "瘢痕疙瘩单细胞RNA测序全流程分析"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    theme: cosmo
    highlight: tango
---

```{r full-analysis, echo=TRUE, message=FALSE, warning=FALSE, results='asis', fig.width=10, fig.height=8, cache=TRUE}
# 1. 加载所有需要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(reticulate)
library(reshape2)
library(circlize)
set.seed(1234)

# 2. 数据加载与预处理（替换为真实数据路径）
data_dir <- "path/to/GSE163973/filtered_feature_bc_matrix"
sc_data <- Read10X(data.dir = data_dir)
sc_obj <- CreateSeuratObject(
  counts = sc_data,
  project = "Keloid",
  min.cells = 3,
  min.features = 200
)

# 3. 添加样本信息（根据实际数据调整）
sc_obj$sample <- c(rep("NS1", 5000), rep("NS2", 5000), rep("NS3", 5000),
                  rep("KF1", 5000), rep("KF2", 5000), rep("KF3", 5000))
sc_obj$group <- ifelse(grepl("KF", sc_obj$sample), "Keloid", "Normal")

# 4. 质量控制
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# 5. 数据标准化
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj, vars.to.regress = c("percent.mt"))

# 6. 批次校正
sc_obj <- RunPCA(sc_obj, npcs = 50)
sc_obj <- RunHarmony(sc_obj, group.by.vars = "sample")

# 7. 细胞聚类
sc_obj <- FindNeighbors(sc_obj, reduction = "harmony", dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.8)
sc_obj <- RunUMAP(sc_obj, reduction = "harmony", dims = 1:15)

# 8. 细胞类型注释
markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "DCN"),
  Endothelial = c("ENG", "PECAM1", "VWF"),
  Immune = c("PTPRC", "CD3E", "CD79A")
)
FeaturePlot(sc_obj, features = unlist(markers), ncol = 3)
sc_obj$celltype <- plyr::mapvalues(sc_obj$seurat_clusters, 
                                  from = 0:20, 
                                  to = c(rep("Fibroblast", 5), "Endothelial", rep("Immune", 3), rep("Other", 12)))

# 9. 成纤维细胞亚群分析
fibroblasts <- subset(sc_obj, idents = c("Fibroblast"))
fibroblasts <- FindVariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts)
fibroblasts <- RunPCA(fibroblasts, npcs = 20)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.6)
fibroblasts <- RunUMAP(fibroblasts, dims = 1:15)

# 10. 细胞通讯分析准备
write.table(as.matrix(fibroblasts@assays$RNA@data), "cpdb_counts.txt", sep ="\t", quote=F)
write.table(data.frame(Cell=rownames(fibroblasts@meta.data), 
                      cell_type=fibroblasts$seurat_clusters),
            "cpdb_meta.txt", sep="\t", quote=F, row.names=F)

# 11. 可视化结果
p1 <- DimPlot(sc_obj, reduction = "umap", group.by = "celltype", label = TRUE)
p2 <- DimPlot(fibroblasts, reduction = "umap", label = TRUE)
plot_grid(p1, p2, ncol = 2)

# 12. 保存结果
saveRDS(sc_obj, "keloid_scRNAseq_processed.rds")
saveRDS(fibroblasts, "fibroblasts_subclustered.rds")
~~~

### `reticulate`

~~~R
---
title: "整合CellPhoneDB分析"
output: html_document
---

```{r setup, include=FALSE}
library(reticulate)
use_python("/path/to/your/python")  # 指定包含cellphonedb的Python环境

# 生成CellPhoneDB所需输入文件
write.table(as.matrix(seurat_obj@assays$RNA@data), "cpdb_counts.txt", sep="\t", quote=F)
write.table(data.frame(Cell=rownames(seurat_obj@meta.data), 
                      cell_type=seurat_obj$celltype),
            "cpdb_meta.txt", sep="\t", quote=F, row.names=F)

# 通过reticulate调用系统命令
system2("cellphonedb", args = c(
  "method", "statistical_analysis",
  "cpdb_meta.txt", "cpdb_counts.txt",
  "--output-path=cellphonedb_out",
  "--threads=4"
))

# 可视化结果（需确保上一步已运行）
if(file.exists("cellphonedb_out/significant_means.txt")) {
  library(ggplot2)
  sig_means <- read.delim("cellphonedb_out/significant_means.txt")
  # 示例可视化
  ggplot(subset(sig_means, means > 0.5), 
         aes(interacting_pair, pair_cell_types, fill=means)) +
    geom_tile() + theme(axis.text.x = element_text(angle=45, hjust=1))
}

###############################################################################
```{python, eval=FALSE}
# Python代码块直接运行CellPhoneDB
import os
os.system("cellphonedb method statistical_analysis cpdb_meta.txt cpdb_counts.txt --output-path=cellphonedb_out")
~~~

