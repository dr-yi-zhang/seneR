#' Title
#'
#' @param data 
#' @param group 
#' @param comparisons 
#' @param map_signif_level 
#'
#' @return
#' @export
#'
#' @examples
plot_group <- function(data, group, feature = "SID_Score",  
                       method = "t.test", comparisons = NULL) {
  # 转换为 symbol 类型
  group_sym <- sym(group)
  feature_sym <- sym(feature)
  
  # 动态计算 y 的范围
  max_y <- max(data[[feature]], na.rm = TRUE)
  step_y <- (max_y - min(data[[feature]], na.rm = TRUE)) * 0.1
  
  # 初始化 ggplot 对象
  p <- ggplot(data, aes(x = !!group_sym, y = !!feature_sym)) +
    stat_boxplot(geom = "errorbar", aes(color = !!group_sym), width = 0.2) +
    geom_boxplot(aes(color = !!group_sym), outlier.shape = NA) +
    geom_jitter(aes(color = !!group_sym), width = 0.1, size = 1) +
    theme_bw() +
    ylab(feature)
  
  # 如果用户指定了 comparisons，则添加显著性标记
  if (!is.null(comparisons)) {
    label_y_positions <- seq(max_y + step_y, by = step_y, length.out = length(comparisons))
    p <- p + ggpubr::stat_compare_means(
      method = method,
      comparisons = comparisons, # 使用用户指定的比较组
      label = "p.signif",
      label.y = label_y_positions # 动态分配显著性标签的 y 位置
    )
  }
  
  return(p)
}




#' Title
#'
#' @param group1 
#' @param group2 
#'
#' @return
#' @export
#'
#' @examples
plot_cellratio <- function(group1,group2){
  library(reshape2)
  prop_df <- table(group1,group2) %>% melt()
  colnames(prop_df) <- c("Cluster","Sample","Number")
  prop_df$Cluster <- factor(prop_df$Cluster)
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #处理后有73种差异还比较明显的颜色，基本够用
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
  
  sample_color <- col_vector[1:10] 
  
  ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
    geom_bar(stat = "identity", width=0.8,position="fill")+
    scale_fill_manual(values=col_vector[1:20]) +
    theme_bw()+
    theme(panel.grid =element_blank()) +
    labs(x="",y="Ratio")+
    ####用来将y轴移动位置
    theme(axis.text.y = element_text(size=12, colour = "black"))+
    theme(axis.text.x = element_text(size=12, colour = "black"))+
    theme(
      axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
    ) 
}

#' Title
#'
#' @param expr 
#' @param meta 
#'
#' @return
#' @export
#'
#' @examples
plot_violin <- function(expr,meta,group){
  expr_long <- rownames_to_column(expr, var = "Feature") %>% 
    pivot_longer(cols = -Feature,names_to = "Sample", values_to = "Expression")
  meta_table <- rownames_to_column(meta, var = "Sample")
  plot_data <- expr_long %>%
    left_join(meta_table, by = "Sample")
  # 绘制分组小提琴图
  ggplot(plot_data, aes(x = Feature, y = Expression, fill = !!sym(group))) +
    geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
    geom_vline(xintercept = seq(1.5, length(unique(plot_data$Feature)) - 0.5, 1), 
               linetype = "dashed", color = "gray", size = 0.5) +
    stat_summary(fun.data = mean_se, geom = "pointrange", 
                 position = position_dodge(width = 0.9), color = "black") +
    theme_minimal() +
    labs(x = "Feature", y = "Expression")+ 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
      panel.grid.major.x = element_blank(),  # 隐藏默认网格线
      panel.grid.minor.x = element_blank()
    )+
    ggpubr::stat_compare_means(aes(group = !!sym(group)), method = "t.test", label = "p.signif") # 添加显著性标记
}

plot_heatmap <- function(gsva_result,meta){
  pheatmap(gsva_result, annotation_col = meta,
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}


