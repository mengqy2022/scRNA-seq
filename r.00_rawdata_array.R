rm(list = ls())
gc()
setwd("D:/workplace/workplace_git/scRNA-seq/02_Data_processing")
if (! dir.exists("./00_rawdata")){
	dir.create("./00_rawdata")
}
setwd("./00_rawdata")
options(stringsAsFactors = FALSE)

# 软件包安装 -------------------------------------------------------------------

# 定义需要安装的包
required_packages <- c("GEOquery", "tidyverse")

# 检查哪些包未安装
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# 如果有缺失的包，安装它们
if (length(missing_packages) > 0) {
	message("正在安装以下缺失的包: ", paste(missing_packages, collapse = ", "))
	
	# 设置CRAN镜像
	options(repos = c(CRAN = "https://cloud.r-project.org/"))
	#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

	# 检查是否需要安装 BiocManager
	if ("GEOquery" %in% missing_packages && !requireNamespace("BiocManager", quietly = TRUE)) {
		message("正在安装 BiocManager...")
		install.packages("BiocManager")
	}
	
	# 安装CRAN包
	cran_packages <- missing_packages[missing_packages != "GEOquery"]
	if (length(cran_packages) > 0) {
		install.packages(cran_packages)
	}
	
	# 安装Bioconductor包
	if ("GEOquery" %in% missing_packages) {
		BiocManager::install("GEOquery")
	}
} else {
	message("所有需要的包已经安装。")
}

# 加载所有包
lapply(required_packages, function(pkg) {
	if (!require(pkg, character.only = TRUE)) {
		warning("无法加载包: ", pkg)
	} else {
		message("成功加载包: ", pkg)
	}
})

# scRNA-seq数据下载 -----------------------------------------------------------
GEO_data <- 'GSE163973'
#gene_annotation <- 'GPL24676'

# 下载GEO数据
# 在服务器下载
# wget -q -c -t 0 --content-disposition "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE163973&format=file" &
getGEOSuppFiles(GEO_data, makeDirectory = TRUE)

# 解压数据
untar("./GSE163973/GSE163973_RAW.tar", exdir = "./GSE163973")

# 获取所有样本的 .tar.gz 文件路径
tar_files <- list.files("./GSE163973", pattern = "\\.tar\\.gz$", full.names = TRUE)

# 为每个样本创建独立目录并解压
for (tar_file in tar_files) {
	# 提取样本名-GSM4994379_KL1_matrix
	sample_name <- sub("\\.tar\\.gz$", "", basename(tar_file))
	sample_dir <- file.path("./GSE163973", sample_name)
	
	# 创建样本目录（如果不存在）
	if (!dir.exists(sample_dir)) dir.create(sample_dir)
	
	# 解压 .tar.gz 到样本目录
	untar(tar_file, exdir = sample_dir)
	cat("已解压:", tar_file, "到", sample_dir, "\n")
}

if (length(tar_files) > 0) {
	file.remove(tar_files)
	message("已删除以下文件：\n", paste(tar_files, collapse = "\n"))
} else {
	message("没有找到 .tar.gz 文件可删除。")
}

# 遍历所有样本目录
#sample_dirs <- list.files("./GSE163973", pattern = "^GSM[0-9].*[^.gz]$", full.names = TRUE)
# for (sample_dir in sample_dirs) {
# 	# 列出当前样本目录下的 .gz 文件
# 	gz_files <- list.files(sample_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
# 	
#   # 解压每个 .gz 文件
# 	for (gz_file in gz_files) {
# 		# 生成解压后的文件名（移除 .gz）
# 		output_file <- sub("\\.gz$", "", gz_file)
# 
# 		# 使用 R.utils::gunzip 解压（需安装 R.utils 包）
# 		if (!requireNamespace("R.utils", quietly = TRUE)) {
# 			install.packages("R.utils")
# 		}
# 		R.utils::gunzip(gz_file, destname = output_file, remove = FALSE)
# 
# 		cat("已解压:", gz_file, "到", output_file, "\n")
# }

# 获取所有.gz文件路径
gz_files <- list.files(path = "./GSE163973", pattern = "\\.gz$", recursive = TRUE, full.names = TRUE)

# 移动文件到对应的_matrix文件夹
for (file in gz_files) {
	# 提取目标文件夹路径（_matrix文件夹）
	target_dir <- dirname(dirname(file))
	
	# 如果文件不在目标文件夹中，则移动
	if (dirname(file) != target_dir) {
		file.rename(from = file, to = file.path(target_dir, basename(file)))
	}
}

# 清除空白文件夹
remove_empty_dirs <- function(path) {
	dirs <- list.dirs(path = path, recursive = TRUE, full.names = TRUE)
	dirs <- rev(dirs)  # 从最深层的子目录开始检查
	
	for (dir in dirs) {
		if (dir == path) next  # 跳过根目录
		
		# 如果只有 "." 和 ".."（即 length <= 2），删除
		if (length(dir(dir, all.files = TRUE)) <= 2) {
			unlink(dir, recursive = TRUE)
			message("Removed empty directory: ", dir)
		}
	}
}

# 在GSE163973目录下执行
remove_empty_dirs("./GSE163973")
