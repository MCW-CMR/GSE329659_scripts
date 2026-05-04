remotes::install_github("bnprks/BPCells/r")
install.packages("scCustomize")
install.packages("devtools")
install.packages("rlang")
devtools::install_github('immunogenomics/presto')

library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(BPCells)
library(scCustomize)
library(writexl)


all_int <- readRDS("../RDS_new_int.rds")
all_int <- JoinLayers(all_int)

write_matrix_dir(mat = all_int[["Spatial"]]$counts, dir = "merged_counts")
counts.mat <- open_matrix_dir(dir = "merged_counts/")
all_int[["Spatial"]]$counts <- counts.mat

################################################################################
##Integrate using CcA, hoping for better integration
################################################################################
DefaultAssay(all_int)
all_int[["Spatial"]] <- split(all_int[["Spatial"]], f = all_int$orig.ident) #Can't integrate an object with joined layers...duh

all_int <- IntegrateLayers(object = all_int, method = CCAIntegration,
                           orig.reduction = "pca", new.reduction = "integrated.cca")
all_int  <- FindNeighbors(all_int, reduction = "integrated.cca", dims = 1:30)
all_int <- FindClusters(all_int, resolution = 2, cluster.name = "cca_clusters", 
                        verbose = FALSE)
all_int <- RunUMAP(all_int, reduction = "integrated.cca", dims = 1:30, 
                   reduction.name = "umap.cca") 

DimPlot(all_int, reduction = "umap.cca", group.by = "orig.ident",
        combine = TRUE, label.size = 2) + 
  theme(axis.text = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  ggtitle("CCA Integration") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) 
ggsave("fig_cca_integrated_orig.ident_logNorm.pdf",
       height = 10, width = 16, 
       dpi = 300, device = "pdf")

DimPlot(all_int, reduction = "umap.cca", group.by = "cca_clusters", split.by = "infection_state",
        combine = TRUE, label.size = 2) + 
  theme(axis.text = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  ggtitle("CCA Integration") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) 
ggsave("fig_cca_integrated_orig.ident_logNorm.pdf",
       height = 10, width = 16, 
       dpi = 300, device = "pdf")

DimPlot(all_int, reduction = "umap.cca", group.by = c("cca_clusters", "broad_celltype"), 
        split.by = "infection_state",
        combine = TRUE, label.size = 2)



################################################################################
##Refine macrophages such that any that express CD4(ENSG00000010610) are relabeled as T-cells
################################################################################
all_int <- JoinLayers(all_int)
Idents(all_int) <- all_int$broad_celltype
all_int@meta.data$broad_celltype_refined <- all_int@meta.data$broad_celltype
head(all_int@meta.data)
unique(all_int@meta.data$broad_celltype)
Idents(all_int) <- all_int$broad_celltype_refined


cd4_pos_mac <- WhichCells(all_int, idents = "Macrophage", expression = ENSG00000010610 > 0)

cd4_pos_mac <- as.data.frame(cd4_pos_mac)
cd4_pos_mac$tcell <- "T-cell"
row.names(cd4_pos_mac) <- cd4_pos_mac$cd4_pos_mac
cd4_pos_mac$cd4_pos_mac <- NULL # 
head(cd4_pos_mac)
dim(cd4_pos_mac)
all_int <- AddMetaData(all_int, cd4_pos_mac, col.name = "broad_celltype_refined")
head(all_int@meta.data)

table(all_int@meta.data$broad_celltype, 
      all_int@meta.data$orig.ident)

table(all_int@meta.data$broad_celltype_refined, 
      all_int@meta.data$orig.ident)

#Write supplemental table 1 to file
write.csv(table(all_int@meta.data$broad_celltype_refined, 
      all_int@meta.data$infection_state), "rev_supp_table1.csv")

#########################################################################################
## Make Fig 5 (umap) and Supp Fig 4 (spatial dimplot) and Fig 9 (spatial CXCL10 & ISG15)
#########################################################################################

Idents(all_int) <- all_int$broad_celltype_refined
all_int$broad_celltype_refined <- factor(all_int$broad_celltype_refined, 
                                         levels = c("B-cell","CTB", "Decidual","Endothelial", "EVT", "LED",
                                                    "Macrophage", "npiCTB", "STB", "Stromal", "T-cell"))
my_groups <- levels(all_int$broad_celltype_refined)
colors <- c("#7F3C8D","#11A579","#3969AC","#A5AA99","#E73F74",
            "#80BA5A","#E68310","#008695","#661100","#f97b72","#F2B701")
names(colors) <- my_groups

##Fig 5
dim_ct_recolor <- DimPlot(all_int, reduction = "umap.cca", 
                          cols = colors, pt.size = 1,
                          label = FALSE, split.by = "infection_state", 
                          group.by = "broad_celltype",
                          repel = TRUE)
dim_ct_recolor <- dim_ct_recolor + theme(
  axis.text = element_text(size = 18),
  plot.title = element_text(size = 20)) +
  ggtitle("Broad cell types")
LabelClusters(dim_ct_recolor, id = "broad_celltype", fontface = "bold")
ggsave(filename = "rev_fig5_dimplot_broad-celltypes_infState.png", height = 7, 
       width = 12, dpi = 300, device = "png")
ggsave(filename = "rev_fig5_dimplot_broad-celltypes_infState.tiff", height = 7, 
       width = 12, dpi = 300, device = "tiff")
ggsave(filename = "rev_fig5_dimplot_broad-celltypes_infState.pdf", height = 7, 
       width = 12, dpi = 300, device = "pdf")

##Supp Fig 4
all_int <- UpdateSeuratObject(object = all_int)
ct.A1 <- SpatialDimPlot(all_int, group.by = "broad_celltype_refined", pt.size.factor = 5, images = "A1",
                        combine = T,
                        cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.A1.2 <- SpatialDimPlot(all_int, group.by = "broad_celltype_refined", pt.size.factor = 4, images = "A1.2",
                          combine = T,
                          cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.D1 <- SpatialDimPlot(all_int, group.by = "broad_celltype_refined", pt.size.factor = 2, images = "D1",
                        combine = T,
                        cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.D1.4 <- SpatialDimPlot(all_int, group.by = "broad_celltype_refined", pt.size.factor = 4, images = "D1.4",
                          combine = T,
                          cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

wrap_plots(ct.A1, ct.A1.2, ct.D1, ct.D1.4)

ggsave("spatial_plots/celltype_refined_spatialDimPlot.png", device = "png",
       height = 10, width = 15)
ggsave("spatial_plots/celltype_refined_spatialDimPlot.svg", device = "svg",
       height = 10, width = 15)


##Fig 9
cxcl10A1 <- SpatialFeaturePlot(all_int, features = "ENSG00000169245", images = "A1", #CXCL10
                               pt.size.factor = 5, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10A1.2 <- SpatialFeaturePlot(all_int, features = "ENSG00000169245", images = "A1.2", #CXCL10
                                 pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10D1 <- SpatialFeaturePlot(all_int, features = "ENSG00000169245", images = "D1", #CXCL10
                               pt.size.factor = 2, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10D1.4 <- SpatialFeaturePlot(all_int, features = "ENSG00000169245", images = "D1.4", #CXCL10
                                 pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
wrap_plots(cxcl10A1, cxcl10A1.2, cxcl10D1, cxcl10D1.4)

ggsave("spatial_plots/Fig9_CXCL10_spatialFeatures.svg", device = "svg",
       height = 10, width = 10)


isg15A1 <- SpatialFeaturePlot(all_int, features = "ENSG00000187608", images = "A1", #isg15
                              pt.size.factor = 5, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15A1.2 <- SpatialFeaturePlot(all_int, features = "ENSG00000187608", images = "A1.2", #isg15
                                pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15D1 <- SpatialFeaturePlot(all_int, features = "ENSG00000187608", images = "D1", #isg15
                              pt.size.factor = 2, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15D1.4 <- SpatialFeaturePlot(all_int, features = "ENSG00000187608", images = "D1.4", #isg15
                                pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
wrap_plots(isg15A1, isg15A1.2, isg15D1, isg15D1.4)

ggsave("spatial_plots/Fig9B_ISG15_spatialFeatures.svg", device = "svg",
       height = 10, width = 10)


#########################################################################################
## Make Fig 6 (heatmaps of ISGs) 
#########################################################################################

ISG.for.heat <- read.csv("../desired_violin_ravit/final_unique_desired_violin.csv")
dim(ISG.for.heat)
ISG.for.heat[1:5,1:2]
DoHeatmap(all_int, features = ISG.for.heat$Unique, 
          group.by = "infection_state", group.bar = TRUE, combine = TRUE,
          label = FALSE, angle = 0) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = ISG.for.heat$Unique,
    labels = ISG.for.heat$UniqueGeneID) +
  theme(axis.text.y = element_text(size = 18),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) 

ggsave("Fig6_heatmap_ISG_genes_logNorm_infection-state.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("Fig6_heatmap_ISG_genes_logNorm_infection-state.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

ggsave("Fig6_heatmap_ISG_genes_logNorm_infection-state_wLabels.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("Fig6_heatmap_ISG_genes_logNorm_infection-state_wLabels.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

DoHeatmap(all_int, features = ISG.for.heat$Unique, 
          group.by = "broad_celltype_refined", group.bar = TRUE, group.colors = colors,
          combine = TRUE, #lines.width = 100,
          label = FALSE, #angle = 90
) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = ISG.for.heat$Unique,
    labels = ISG.for.heat$UniqueGeneID) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme(axis.text.y = element_text(size = 18),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

ggsave("Fig6_heatmap_ISG_genes_logNorm_celltypes.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("Fig6_heatmap_ISG_genes_logNorm_celltypes.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

ggsave("Fig6_heatmap_ISG_genes_logNorm_celltypes_wLabels_v4.png",
       height = 22, width = 20, 
       dpi = 300, device = "png")
ggsave("Fig6_heatmap_ISG_genes_logNorm_celltypes_v2.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("Fig6_heatmap_ISG_genes_logNorm_celltypes_v2.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")
ggsave("heatmap_ISG_genes_logNorm_celltypes_wLabels.svg",
       height = 18, width = 10, 
       dpi = 300, device = "svg")

#########################################################################################
## Remake VlnPlts for Fig 7 and Fig 8
#########################################################################################
##Fig 7
#Get genes of interest 
goi_f7 <- read.csv("genes_for_vln_fig7.csv")
goi_f7

goi_f7_vect <- setNames(goi_f7$UniqueGeneID, goi_f7$Unique)

fig7 <- Stacked_VlnPlot(all_int, features = goi_f7$Unique, split.by = "infection_state",
                        pt.size = 0.5, colors_use = c("#F7766D","#00BFC4"), x_lab_rotate = 45) 
  
fig7.1 <- Reduce('+', lapply(seq_along(fig7), function(i) {
  fig7[[i]] + labs(y = unname(goi_f7_vect[i]))
}))

fig7.2 <- fig7.1 &
  theme(axis.text = element_text(size = 16))
fig7.2

ggsave("VlnPlots_revisions/rev_Fig7_vlnplts.svg",
       height = 12.5, width = 10.4, units = "in",
       dpi = 300, device = "svg")

##Fig 8
#Get genes of interest 
goi_f8 <- read.csv("genes_for_vln_fig8.csv")
goi_f8
goi_f8_vect <- setNames(goi_f8$UniqueGeneID, goi_f8$Unique)
unlist(split(goi_f8$UniqueGeneID, goi_f8$Unique))

goi_f8_vect <- goi_f8$UniqueGeneID
names(goi_f8_vect) <- goi_f8$Unique
 
fig8 <- Stacked_VlnPlot(all_int, features = goi_f8$Unique, split.by = "infection_state",
                pt.size = 0.5, colors_use = c("#F8766D","#00BFC4"), x_lab_rotate = 45) 
fig8.1 <- Reduce('+', lapply(seq_along(fig8), function(i) {
  fig8[[i]] + labs(y = unname(goi_f8_vect[i]))
}))

fig8.2 <- fig8.1 &
  theme(axis.text = element_text(size = 16))
fig8.2

ggsave("VlnPlots_revisions/rev_Fig8_vlnplts.svg",
       height = 9, width = 10.4, 
       dpi = 300, device = "svg")

##New Supp Fig 6
goi_sf5 <- read.csv("genes_for_new_vln_SuppFig6.csv")
goi_sf5

goi_sf5_vect <- setNames(goi_sf5$UniqueGeneID, goi_sf5$Unique)
suppfig6 <- Stacked_VlnPlot(all_int, features = goi_sf5$Unique, split.by = "infection_state",
                        pt.size = 0.5, colors_use = c("#F8766D","#00BFC4"), x_lab_rotate = 45) 
suppfig6.1 <- Reduce('+', lapply(seq_along(suppfig6), function(i) {
  suppfig6[[i]] + labs(y = unname(goi_sf5_vect[i]))
}))

suppfig6.2 <- suppfig6.1 &
  theme(axis.text = element_text(size = 16))
suppfig6.2

ggsave("VlnPlots_revisions/rev_new_suppFig6_vlnplts_fixed.svg",
       height = 12, width = 10.4, 
       dpi = 300, device = "svg")


#########################################################################################
## Rerun diff expression in order to remake supp table 2
#########################################################################################
Idents(all_int) <- all_int$broad_celltype_refined

de_markers_CTB <- FindMarkers(all_int, ident.1 = "infected", 
                              group.by = "infection_state", 
                              subset.ident = "CTB")

de_markers_Stromal <- FindMarkers(all_int, ident.1 = "infected", 
                                  group.by = "infection_state", 
                                  subset.ident = "Stromal")

de_markers_tCell <- FindMarkers(all_int, ident.1 = "infected", 
                                group.by = "infection_state",
                                subset.ident = "T-cell")

de_markers_Decidual <- FindMarkers(all_int, ident.1 = "infected", 
                                   group.by = "infection_state",
                                   subset.ident = "Decidual")

de_markers_Endothelial <- FindMarkers(all_int, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "Endothelial")

de_markers_STB <- FindMarkers(all_int, ident.1 = "infected", 
                              group.by = "infection_state",
                              subset.ident = "STB")

de_markers_Macrophage <- FindMarkers(all_int, ident.1 = "infected", 
                                     group.by = "infection_state",
                                     subset.ident = "Macrophage")

de_markers_EVT <- FindMarkers(all_int, ident.1 = "infected", 
                              group.by = "infection_state",
                              subset.ident = "EVT")

de_markers_LED <- FindMarkers(all_int, ident.1 = "infected", 
                              group.by = "infection_state",
                              subset.ident = "LED")

de_markers_bCell <- FindMarkers(all_int, ident.1 = "infected", ##Only have b-cells in infected, so test can't be run
                                group.by = "infection_state",
                                subset.ident = "B-cell")

de_markers_npiCTB <- FindMarkers(all_int, ident.1 = "infected", 
                                 group.by = "infection_state",
                                 subset.ident = "npiCTB")

de_markers_CTB[de_markers_CTB$p_val_adj<0.05,]
list_de_markers_updated <- list("de_markers_Endothelial", "de_markers_Decidual", "de_markers_EVT",
                                "de_markers_LED", "de_markers_Macrophage", "de_markers_STB",
                                "de_markers_tCell", "de_markers_Stromal", "de_markers_CTB")

conversion <- read.csv("../out_ensembl_geneSymbol_conversion.csv", header = TRUE)
conversion[1:5,1:2]

for (i in list_de_markers_updated) { #for file in list
  n <- get(sprintf("%s", i)) #get the gene table for each cluster and save it in n
  #print(n)
  n <- n[n$p_val_adj<0.05,]
  n$GeneSymbol <- conversion$GeneSymbol[match(rownames(n),conversion$Ensembl)]
  #print(i)
  #print(head(n))
  write.csv(n, sprintf("diff_expression/%s.csv", i)) #write n to a csv file
}




########################################################################################################
## Requested DotPlot
########################################################################################################
Idents(all_int) <- all_int$broad_celltype_refined

DotPlot(all_int, features = c("ENSG00000135346", "ENSG00000116183","ENSG00000129226", 
                              "ENSG00000010610", #"ENSG00000168685",
                              "ENSG00000110799")) + #, "ENSG00000153563"
  scale_x_discrete(breaks = c("ENSG00000135346", "ENSG00000116183","ENSG00000129226", 
                              "ENSG00000010610", #"ENSG00000168685",
                              "ENSG00000110799"), #, "ENSG00000153563"
                   labels = c("CGA","PAPPA2", "CD68", 
                              "CD4", #"IL7R", 
                              "VWF")) + #, "CD8A"
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

ggsave("rev_dotplot_distinguishing_genes_celltypes.png",
       height = 15, width = 15, 
       dpi = 300, device = "png")
ggsave("rev_dotplot_distinguishing_genes_celltypes.tiff",
       height = 15, width = 15, 
       dpi = 300, device = "tiff")
ggsave("rev_dotplot_distinguishing_genes_celltypes.pdf",
       height = 15, width = 15, 
       dpi = 300, device = "pdf")

#split by infection
all_int$infect_celltype <- paste0((Idents(all_int)), "_", all_int$infection_state)
head(all_int@meta.data)

DotPlot(all_int, features = c("ENSG00000213949","ENSG00000161638","ENSG00000116183","ENSG00000087245"),
        group.by = "infect_celltype") +
  scale_x_discrete(breaks = c("ENSG00000213949","ENSG00000161638","ENSG00000116183","ENSG00000087245"),
                   labels = c("ITGA1","ITGA5","PAPPA2","MMP2")) &
  #guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme(axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

ggsave("rev_Fig_dotplot_celltypes_infection_ITGA1_etc.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("rev_Fig_dotplot_celltypes_infection_ITGA1_etc.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

ENSG00000183878 #UTY
ENSG00000067048 #DDX3Y

VlnPlot(all_int, features = c("ENSG00000183878","ENSG00000067048"), 
        group.by = "orig.ident", combine = TRUE) 
  

FindMarkers(all_int, features = c("ENSG00000183878","ENSG00000067048"),
            group.by = "infection_state", ident.1 = "infected")
FindMarkers(all_int, features = c("ENSG00000183878","ENSG00000067048"),
            group.by = "orig.ident", ident.1 = "infect_48555")


saveRDS(all_int, "RDS_Boger_placenta_integrated.rds")
