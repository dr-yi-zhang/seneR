# seneR/R/env_Load.R

#' @title Load seneR environment
#' @description Initializes or checks the Python environment required by seneR.
#' @return Invisible NULL.
#' @export
# seneR/R/env_Load.R

#' @title Load seneR environment
#' @description Activates the reticulate environment (for local use).
#' @return Invisible NULL.
#' @export
env_Load <- function(local=TRUE) {
  
  # For deployment, this function does nothing as py_require() handles it.
  # For local use, it ensures the virtualenv is used, which is good practice.
  if(local){reticulate::virtualenv_install("seneR_env", packages = c("scikit-learn==1.5.0", "joblib",'NumPy==1.26.4'))
    }
  if (reticulate::virtualenv_exists("seneR_env")) {
    # If the env exists locally (as defined in your Rprofile), use it.
    reticulate::use_virtualenv("seneR_env", required = TRUE)
    print("Local 'seneR_env' activated.")
  } else {
    # If not running locally, this message confirms server is relying on py_require().
    print("Environment setup relies on py_require() in app.R for remote deployment.")
  }
  
  return(invisible(NULL))
}
#' Title
#'
#' @param adata 
#' @param sidnums 
#' @param denoising 
#' @param threads 
#'
#' @return
#' @export
#'
#' @examples
SenCID <- function(data,sidnums=c(1,2,3,4,5,6), binarize=FALSE, denoising,threads){
  #数据预处理
  data_scaled <- scale_data(data)
  
  # 计算SID分数
  pred_list <- lapply(sidnums, function(sidnum) Pred(data_scaled, sidnum, binarize))
  
  # 创建 SID 分数字典
  sid_list <- paste0("SID", sidnums)
  pred_dict <- setNames(pred_list, sid_list)
  
  # 计算推荐 SID
  recSID <- Recommend(data_scaled)
  # score_SID <- names(which.max(table(recSID$RecSID)))
  score_res <- NULL
  for (i in 1:length(recSID$RecSID)) {
    recScoreLine <- pred_dict[[recSID$RecSID[i]]][i,]
    recScoreLine$RecSID <- recSID$RecSID[i]
    score_res <- rbind(score_res,recScoreLine)
  }
  
  cat("Finished. Giving SID scores and SID Recommendation...\n")
  
  # 返回结果
  return(list(pred_dict = pred_dict, recSID = recSID,score_res=score_res))
}

#' Title
#'
#' @param expr 
#'
#' @return
#' @export
#'
#' @examples
seneGSVA <- function(expr){
  geneset <- system.file("source", "sene_related.gmt", package = 'seneR')
  geneset <- GSEABase::getGmt(geneset)
  #创建一个参数对象
  params <- GSVA::gsvaParam(as.matrix(expr), geneset, kcdf = "Poisson")
  #执行 GSVA 分析
  gsva_result <- GSVA::gsva(params)
  gsva_result <- as.data.frame(gsva_result)
}

#' Title
#'
#' @param seurat_obj 
#' @param color 
#' @param cores 
#'
#' @returns
#' @export
#'
#' @examples
sene_traject <- function(seurat_obj, color="SID_Score", cores=3){
  expr_matrix <- GetAssayData(seurat_obj,layer = "counts")
  pd <- new("AnnotatedDataFrame", data = seurat_obj@meta.data)
  fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(expr_matrix),
                                                    row.names = rownames(expr_matrix)))
  CDS_obj <- newCellDataSet(expr_matrix,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=negbinomial.size())
  CDS_obj <- estimateSizeFactors(CDS_obj)
  CDS_obj <- estimateDispersions(CDS_obj)
  CDS_obj <- detectGenes(CDS_obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(CDS_obj),
                                      num_cells_expressed >= 10))
  diff_test_res <- differentialGeneTest(CDS_obj[expressed_genes,],
                                        fullModelFormulaStr = "~sm.ns(SID_Score)",
                                        cores = cores)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  CDS_obj <- setOrderingFilter(CDS_obj, ordering_genes)
  CDS_obj <- reduceDimension(CDS_obj, max_components = 2,
                             method = 'DDRTree')
  pData(CDS_obj)$Pseudotime <- pData(CDS_obj)$SID_Score
  p <- plot_cell_trajectory(CDS_obj, color_by = color)
  print(p)
  return(CDS_obj)
}


#' Title
#'
#' @param CDS_obj 
#' @param num 
#' @param num_clusters 
#' @param cores 
#'
#' @returns
#' @export
#'
#' @examples
sene_heatmap <- function(CDS_obj,num = 20,num_clusters = 3, cores=3){
  expressed_genes <- row.names(subset(fData(CDS_obj),
                                      num_cells_expressed >= 10))
  diff_test_res <- differentialGeneTest(CDS_obj[expressed_genes,],
                                        fullModelFormulaStr = "~sm.ns(SID_Score)",
                                        cores = cores)
  sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
  plot_pseudotime_heatmap(CDS_obj[sig_gene_names[1:num],],
                          num_clusters = num_clusters,
                          cores = cores,
                          show_rownames = T)
}
