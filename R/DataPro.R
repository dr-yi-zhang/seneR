#' Title
#'
#' @param data 
#' @param denoising 
#' @param threads 
#' @param savetmp 
#'
#' @return
#' @export
#'
#' @examples
scale_data <- function(data, denoising = TRUE, threads = 1, savetmp = FALSE) {
  # 1. 读取基因列表
  sene_path <- system.file("source", "seneset.txt", package = "seneR")
  sene <- read.table(sene_path, sep = "\t", header = FALSE)
  colnames(sene) <- c("Gene")
  
  # 2. 确保数据中的基因与给定基因列表中的基因有交集
  seneGenes <- intersect(rownames(data), sene$Gene)
  if (length(seneGenes) <= 10) {
    stop("Features of count data should be hgnc symbol: BIRC5, BRCA1, etc")
  }
  
  # 选择感兴趣的基因
  data_ae <- data[seneGenes, , drop = FALSE]
  
  # 3. 基因过滤，移除低表达的基因
  data_ae <- data_ae[Matrix::rowSums(data_ae) > 0, ]
  
  #  Cpm_added
  cpm_added <- Cpm_added(data_ae, sene$Gene)
  
  # Zcol
  cpm_zcol = Zcol(cpm_added)
  
  #数据标准化：Seurat 的 ScaleData 函数
  # data_ae <- ScaleData(data_ae, features = seneGenes)
  
  # 先不保存临时数据
  
  return(cpm_zcol)
}


Cpm_added <- function(count, features) {
  # 找出 count 矩阵中没有出现的基因
  addrows <- setdiff(features, rownames(count))
  # 创建一个包含缺失基因的矩阵，填充小值 0.01
  addmat <- matrix(0.01, nrow = length(addrows), ncol = ncol(count),
                   dimnames = list(addrows, colnames(count)))
  # 将缺失基因矩阵与原始 count 矩阵合并
  cpm_added <- rbind(count, addmat)
  # 将每一列标准化为每百万条计数 (CPM)
  cpm_added <- apply(cpm_added, 2, function(x) x / sum(x) * 1e6)
  # 转换回数据框形式并确保返回结果的基因顺序与 features 一致
  cpm_added <- cpm_added[features, ]
  return(as.data.frame(cpm_added))
}

# 定义 Zcol 函数
Zcol <- function(cpm_added) {
  # Log 转换
  cpm_log <- log2(cpm_added + 1)
  
  # 列标准化
  cpm_zcol <- as.data.frame(scale(cpm_log, center = TRUE, scale = TRUE))
  
  return(cpm_zcol)
}