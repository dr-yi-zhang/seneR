#' Title
#'
#' @param libname 
#' @param pkgname 
#'
#' @return
#' @export
#'
#' @examples
env_Load <- function(libname, pkgname) {
  # 检查虚拟环境是否存在，如果不存在则创建并安装依赖
  if (!reticulate::virtualenv_exists("myenv")) {
    print("Creating virtual environment 'myenv'...")
    reticulate::virtualenv_create("myenv")
  }
  # 激活虚拟环境
  print("Setting up virtual environment...")
  reticulate::use_virtualenv("myenv", required = TRUE)

  reticulate::virtualenv_install("myenv", packages = c("scikit-learn==1.5.0", "joblib",'NumPy==1.26.4'))
  print("Virtual environment setup and dependencies installed.")

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
  score_SID <- names(which.max(table(recSID$RecSID)))
  score_res <- pred_dict[[score_SID]]
  
  cat("Finished. Giving SID scores and SID Recommendation...\n")
  
  # 返回结果
  return(list(pred_dict = pred_dict, recSID = recSID,score_res=score_res,score_SID=score_SID))
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
