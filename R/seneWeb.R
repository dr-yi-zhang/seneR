
# webAnalysis <- function(exprFilePath,metaFilePath){
#   Expr <- read.csv(exprFilePath,row.names = 1)
#   Meta <- read.csv(metaFilePath,row.names = 1)
#   SID_res <- SenCID(Expr)
#   SID_table <- cbind(SID_res$score_res,Meta)
#   SID_plot <- plot_group(SID_table,'group','SID_Score')
#   gsva_res <- seneGSVA(Expr)
#   GSVA_table <- rbind(gsva_res,t(SID_res$score_res))
#   GSVA_plot <- plot_violin(GSVA_table,Meta,'group')
#   webList <- list(SID_table,SID_plot,GSVA_table,GSVA_plot)
# }


#' Title
#'
#' @param FilePath 
#'
#' @returns
#' @export
#'
#' @examples
webAnalysis <- function(uploadDir){
  library(reticulate)
  env_Load()
  exprFilePath <- paste0(uploadDir,'/exprFile.csv')
  Expr <- read.csv(exprFilePath,row.names = 1)
  SID_res <- SenCID(Expr,binarize = T)
  SID_res <- SID_res$score_res
  gsva_res <- seneGSVA(Expr)
  gsva_res <- t(gsva_res)
  sum_res <- cbind(SID_res,gsva_res)
  res_file <- paste0(uploadDir,'/result.csv')
  write.csv(sum_res,res_file)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetBashFullPath
#' @export 
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}


RowName <- function(uploadDir){
  table <- read.csv(paste0(uploadDir,'/result.csv'))
  return(rownames(table));
}

ColName <- function(uploadDir){
  table <- read.csv(paste0(uploadDir,'/result.csv'))
  return(colnames(table));
}

ResCol <- function(uploadDir,colInx){
  table <- read.csv(paste0(uploadDir,'/result.csv'))
  res <- table[[colInx]];
  hit.inx <- is.na(res) | res == ""; # note, must use | for element-wise operation
  res[hit.inx] <- "N/A";
  return(res);
}