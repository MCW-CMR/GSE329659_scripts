library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)

################################################################################
## Sample 63
################################################################################
uninfect_63 <- Load10X_Spatial(data.dir = "../count-63-uninfected/manAln/",
                               filename = "filtered_feature_bc_matrix.h5",
                               assay = "spatial",
                               slice = "A1",
                               filter.matrix = TRUE,
                               to.upper = FALSE
)

uninfect_63
uninfect_63@meta.data$orig.ident <- "uninfect_63"

plot1 <- VlnPlot(uninfect_63, features = "nCount_Spatial")
plot1 <- VlnPlot(uninfect_63)
uninfect_63$nCount_spatial
plot1 <- VlnPlot(uninfect_63, features = "nCount_spatial")
plot1
plot2 <- SpatialFeaturePlot(uninfect_63, features = "nCount_spatial") +
  theme(legend.position = "right")
wrap_plots(plot1, plot2)

uninfect_63 <- SCTransform(uninfect_63, assay = "Spatial", verbose = FALSE)
uninfect_63 <- FindNeighbors(uninfect_63, reduction = "pca", dims = 1:30)
uninfect_63 <- FindClusters(uninfect_63, verbose = FALSE)
uninfect_63 <- RunUMAP(uninfect_63, reduction = "pca", dims = 1:30)
p1 <- DimPlot(uninfect_63, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(uninfect_63, label = TRUE, label.size = 3)
p1 + p2
p2 <- SpatialDimPlot(uninfect_63, label = TRUE, label.size = 3, pt.size.factor = 5)
p1 + p2
p2 <- SpatialDimPlot(uninfect_63, label = TRUE, label.size = 3, pt.size.factor = 3)
p1 + p2
p2 <- SpatialDimPlot(uninfect_63, label = TRUE, label.size = 3, pt.size.factor = 4)
p1 + p2


saveRDS(uninfect_63,
        file = file.path(".",paste0("uninfect_63", "_SeuratObject.rds")))


################################################################################
## Sample 227
################################################################################
uninfect_227 <- Load10X_Spatial(data.dir = "../count-227-uninfected/manAln/",
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "A1",
                                filter.matrix = TRUE,
                                to.upper = FALSE
)

uninfect_227
uninfect_227@meta.data$orig.ident <- "uninfect_227"

uninfect_227$nCount_Spatial

plot1 <- VlnPlot(uninfect_227, features = "nCount_Spatial")
plot2 <- SpatialFeaturePlot(uninfect_227, features = "nCount_Spatial") +
  theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot2

uninfect_227 <- SCTransform(uninfect_227, assay = "Spatial", verbose = FALSE) #Due to a bug in Seurat at the time this will be switched to log norm in the merge
uninfect_227 <- RunPCA(uninfect_227, assay = "SCT", verbose = FALSE)
uninfect_227 <- FindNeighbors(uninfect_227, reduction = "pca", dims = 1:30)
uninfect_227 <- FindClusters(uninfect_227, verbose = FALSE)
uninfect_227 <- RunUMAP(uninfect_227, reduction = "pca", dims = 1:30)
p1 <- DimPlot(uninfect_227, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(uninfect_227, label = TRUE, label.size = 3, pt.size.factor = 4)
p1 + p2

saveRDS(uninfect_227,
        file = file.path(".",paste0("uninfect_227", "_SeuratObject.rds")))

################################################################################
## Sample 48555
################################################################################
infect_48555 <- Load10X_Spatial(data.dir = "../count-48555-infected/manAln/",
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "D1",
                                filter.matrix = TRUE,
                                to.upper = FALSE
)

infect_48555
infect_48555@meta.data$orig.ident <- "infect_48555"

infect_48555$nCount_Spatial
plot1 <- VlnPlot(infect_48555, features = "nCount_Spatial")
plot2 <- SpatialFeaturePlot(infect_48555, features = "nCount_Spatial") +
  theme(legend.position = "right")
wrap_plots(plot1, plot2)


infect_48555 <- SCTransform(infect_48555, assay = "Spatial", verbose = FALSE) #Due to a bug in Seurat at the time this will be switched to log norm in the merge
infect_48555 <- RunPCA(infect_48555, assay = "SCT", verbose = FALSE)
infect_48555 <- FindNeighbors(infect_48555, reduction = "pca", dims = 1:30)
infect_48555 <- FindClusters(infect_48555, verbose = FALSE)
infect_48555 <- RunUMAP(infect_48555, reduction = "pca", dims = 1:30)
p1 <- DimPlot(infect_48555, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(infect_48555, label = TRUE, label.size = 3, pt.size.factor = 4)
p1 + p2
p2 <- SpatialDimPlot(infect_48555, label = TRUE, label.size = 3, pt.size.factor = 2)
p1 + p2

saveRDS(infect_48555,
        file = file.path(".",paste0("infect_48555", "_SeuratObject.rds")))

################################################################################
## Sample 602
################################################################################

infect_602 <- Load10X_Spatial(data.dir = "../count-602-infected/manAln/",
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "D1",
                              filter.matrix = TRUE,
                              to.upper = FALSE
)

infect_602
infect_602@meta.data$orig.ident <- "infect_602"

infect_602$nCount_Spatial
plot1 <- VlnPlot(infect_602, features = "nCount_Spatial")
plot2 <- SpatialFeaturePlot(infect_602, features = "nCount_Spatial") +
  theme(legend.position = "right")
wrap_plots(plot1, plot2)


infect_602 <- SCTransform(infect_602, assay = "Spatial", verbose = FALSE) #Due to a bug in Seurat at the time this will be switched to log norm in the merge
infect_602 <- RunPCA(infect_602, assay = "SCT", verbose = FALSE)
infect_602 <- FindNeighbors(infect_602, reduction = "pca", dims = 1:30)
infect_602 <- FindClusters(infect_602, verbose = FALSE)
infect_602 <- RunUMAP(infect_602, reduction = "pca", dims = 1:30)
p1 <- DimPlot(infect_602, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(infect_602, label = TRUE, label.size = 3, pt.size.factor = 2)
p1 + p2
p2 <- SpatialDimPlot(infect_602, label = TRUE, label.size = 3, pt.size.factor = 4)
p1 + p2

saveRDS(infect_602,
        file = file.path(".",paste0("infect_602", "_SeuratObject.rds")))












