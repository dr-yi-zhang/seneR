#GSE249229
library(Seurat)
library(seneR)
env_Load()
seurat_obj <- readRDS('GSE249229_count.rds')
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")
sene_meta <- SenCID(counts_matrix,binarize = T)
identical(rownames(sene_meta$recSID),rownames(seurat_obj@meta.data))
table(sene_meta$recSID$RecSID)
seurat_obj@meta.data <- cbind(seurat_obj@meta.data,sene_meta$pred_dict$SID3)

#Single-cell data preprocessing ####
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
DimPlot(seurat_obj, reduction = "pca")
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)  # dims根据ElbowPlot选择
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(Idents(seurat_obj))

chon_marker <- c('CILP2', 'OGN','BNIP3','FAM162A','SPP1',
                 'IBSP','COL10A1','EPYC','FRZB','MMP2',
                 'COL1A1','COL1A2','HMGB2','CDK1','UBE2C',
                 'CCNB1')

library(ggplot2)
DotPlot(seurat_obj, 
        features = chon_marker,
        group.by = "seurat_clusters") + coord_flip()

C_celltype=c('ProC', #0
             'RepC',
             'HTC',
             'EMatC',
             'FibC',
             'FibC',#5
             'EMatC',
             'FibC',
             'ApoC',
             'StrC',
             'HTC',#10
             'ProC',
             'ApoC',
             'ProC'#13
)
Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
names(C_celltype) <- levels(seurat_obj)
seurat_obj<- RenameIdents(seurat_obj, C_celltype)
seurat_obj@meta.data$C_celltype <- Idents(seurat_obj)
Idents(seurat_obj)=seurat_obj@meta.data$C_celltype

colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#bc3c29')

matchID <- match(rownames(seurat_obj@meta.data),rownames(sene_meta$recSID))
matchInfo <- sene_meta$recSID[matchID,]
identical(rownames(sene_meta$recSID),rownames(seurat_obj@meta.data))
seurat_obj@meta.data <- cbind(seurat_obj@meta.data,matchInfo$RecSID)
colnames(seurat_obj@meta.data)


# Reproduction of Figures ####
#F3a
DimPlot(seurat_obj, group.by="C_celltype", label=T, label.size=5,
        pt.size = 1)

#F3b
seurat_obj$tissue_type <- ifelse(substr(seurat_obj$orig.ident, 1, 4)=="CTRL","CTRL","PM")
seurat_obj$Binarization <- ifelse(seurat_obj$Binarization=="0","non-senescent","senescent")
DimPlot(seurat_obj, group.by="Binarization", label=T, label.size=5,
        pt.size = 1,split.by = 'tissue_type')

#F3c
meta <- seurat_obj@meta.data
colnames(meta)
plot_group(meta[,c("SID_Score","Binarization")], meta[,c("orig.ident","tissue_type")], 'orig.ident','SID_Score')


#F3e
slpi_expression <- FetchData(seurat_obj, vars = "SLPI")
slpi_stat <- ifelse(slpi_expression$SLPI==0,'neg','pos')
seurat_obj@meta.data$slpi_exp <- slpi_expression$SLPI
seurat_obj@meta.data$slpi_stat <- slpi_stat
meta <- seurat_obj@meta.data
plot_group(meta[,c("SID_Score","Binarization")], meta[,c("slpi_stat","tissue_type")],group = 'slpi_pn',feature = 'SID_Score')

#F3g
VlnPlot(seurat_obj, features = "SLPI", group.by = "Binarization", pt.size = 0)


#F3h
library(monocle)
library(SeuratData)
seneCDS <- sene_traject(seurat_obj,color = "SID_Score",cores = 7)

#F3i
seneCDS <- sene_traject(seurat_obj,color = "C_celltype",cores = 7)

#F3j
sene_heatmap(seneCDS, num = 20, num_clusters = 3)


