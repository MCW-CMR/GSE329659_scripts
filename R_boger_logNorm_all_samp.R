library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)


uninfect_63 <- readRDS("../../uninfect_63_SeuratObject.rds")
uninfect_227 <- readRDS("../../uninfect_227_SeuratObject.rds")
infect_48555 <- readRDS("../../infect_48555_SeuratObject.rds")
infect_602 <- readRDS("../../infect_602_SeuratObject.rds")

rownames(uninfect_63)
seurat_list <- list(uninfect_227, uninfect_63, infect_48555, infect_602)

##Set up reference
ref <- readRDS("../../../dbGAP_study/ref_BPPV_seurat_SeuratObject.rds")
Idents(ref)
head(ref@meta.data)
length(unique(ref@meta.data$cluster_name))
rownames(ref)
Idents(ref) <- ref@meta.data$cluster_name
DimPlot(ref, reduction = "umap", group.by = "cluster_name", label = TRUE)
DimPlot(ref, reduction = "umap", group.by = "cluster_name", label = FALSE)

DefaultAssay(ref) <- "RNA"
ref[['SCT']] <- NULL
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)


##Merge Boger placental samples
merge_all <- merge(uninfect_63, list(uninfect_227, infect_48555, infect_602), 
                   add.cell.ids = c("uninfect1","uninfect2","infected1","infected2"))
head(merge_all@meta.data)

unique(merge_all@meta.data$seurat_clusters)

DefaultAssay(merge_all) <- "Spatial"
merge_all[['SCT']] <- NULL

#Since I have two treatment groups, want to work with the full merge
merge_all <- NormalizeData(merge_all, assay = "Spatial", verbose = TRUE)
merge_all <- FindVariableFeatures(merge_all)
all.genes <- rownames(merge_all)
merge_all <- ScaleData(merge_all, features = all.genes)

merge_all  <- RunPCA(merge_all, features = VariableFeatures(merge_all),
                     verbose = FALSE)
DimPlot(merge_all, reduction = "pca", group.by = "orig.ident")
merge_all  <- FindNeighbors(merge_all, reduction = "pca", dims = 1:30)
merge_all <- FindClusters(merge_all, verbose = FALSE)
merge_all  <- RunUMAP(merge_all, reduction = "pca", dims = 1:30)

#make treatment variable
infState <- merge_all@meta.data$orig.ident
infState <- gsub("_.*","",infState)
infState <- gsub("t","ted",infState)
all_meta <- merge_all@meta.data
all_meta["infection_state"] <- infState
head(all_meta)
all_meta_trim <- subset(all_meta, select = c("infection_state"))
head(all_meta_trim)
merge_all <- AddMetaData(merge_all, all_meta_trim)
head(merge_all@meta.data)

DimPlot(merge_all, reduction = "umap", label = TRUE, group.by = "orig.ident")
DimPlot(merge_all, reduction = "umap", label = TRUE, group.by = "infection_state")
SpatialDimPlot(merge_all, label = FALSE, label.size = 3, pt.size.factor = 4)

backup_merge_all <- merge_all

##Split object to integrate based on infection type
## Already done apparently
#merge_all[["Spatial"]] <- split(merge_all[["Spatial"]], f = merge_all$infection_state)

##Integrate using harmony & updated v. 5 integration method
new_int <- IntegrateLayers(object = merge_all, method = HarmonyIntegration,
                           orig.reduction = "pca", new.reduction = "harmony")
new_int  <- FindNeighbors(new_int, reduction = "harmony", dims = 1:30)
new_int <- FindClusters(new_int, resolution = 2, cluster.name = "harmony_clusters", 
                        verbose = FALSE)
new_int <- RunUMAP(new_int, reduction = "harmony", dims = 1:30, 
                   reduction.name = "umap.harmony")
DimPlot(new_int, reduction = "umap.harmony", group.by = "harmony_clusters",
        combine = FALSE, label.size = 2, split.by = "infection_state")
DimPlot(new_int, reduction = "umap.harmony", group.by = "harmony_clusters",
        combine = FALSE, label.size = 2)
DimPlot(new_int, reduction = "umap.harmony", group.by = "infection_state",
        combine = FALSE, label.size = 2)
DimPlot(new_int, reduction = "umap.harmony", group.by = "orig.ident",
        combine = FALSE, label.size = 2)

#new_merged_all <- JoinLayers(new_int) ##This doesn't work :(

##Add single cell reference labels to our dataset
anchors_new_int <- FindTransferAnchors(reference = ref, 
                                       query = new_int,
                                       query.assay = "Spatial",
)

predictions.assay_new_int <- TransferData(anchorset = anchors_new_int, 
                                          refdata = ref$cluster_name, 
                                          prediction.assay = TRUE,
                                          dims = 1:30)
unique(rownames(predictions.assay_new_int$data))


#new_int[["predictions"]] <- predictions.assay_new_int ##This still doesn't work, so I did it manually

head(predictions.assay_new_int$data) #find prediction data
predictions.assay_new_int$data[1:10,1:5] #look at prediction data

#Get prediction label to add to metadata
mat <- predictions.assay_new_int$data
mat[1:10,1:5]
#write.csv(mat, "mat.csv")
mat <- t(as.data.frame(mat))
mat[1:5,1:5]
length(colnames(mat))
mat <- as.data.frame(mat)
max_var <- colnames(mat)[max.col(mat, ties.method = "first")]
mat$max <- as.data.frame(max_var)
dim(mat)
mat[1:5,1:6]

##Add prediction label to metadata (integrated)
meta_new_int <- new_int@meta.data
meta_new_int["max_prediction"] <- mat$max
head(meta_new_int)
meta_new_int_trim <- subset(meta_new_int, select = c("max_prediction"))
head(meta_new_int_trim)
new_int <- AddMetaData(new_int, meta_new_int_trim)
head(new_int@meta.data)
head(merge_all@meta.data)

#Fix order so uninfected is first
new_int$infection_state <- factor(new_int$infection_state, 
                                  levels = c("uninfected","infected"))

#look at new umap/harmony plots
DimPlot(new_int, reduction = "umap.harmony", label = TRUE, 
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", label = TRUE, 
        group.by = "infection_state")
DimPlot(new_int, reduction = "umap.harmony", 
        label = TRUE, group.by = "max_prediction",
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", 
        label = TRUE, split.by = "infection_state", 
        group.by = "max_prediction",
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", 
        label = FALSE, split.by = "infection_state", 
        group.by = "max_prediction",
        repel = TRUE)
ggsave("integrated_labeled_infection-state_logNorm.png",
       height = 10, width = 16, 
       dpi = 300, device = "png")
ggsave("integrated_unlabeled_infection-state_logNorm.png",
       height = 10, width = 16, 
       dpi = 300, device = "png")

#Assign integrated clusters to merged_all cell names
#Assuming cell orders have remained the same -> will check if that's true later
head(meta_new_int_trim)
head(merge_all)
merge_all <- AddMetaData(merge_all, meta_new_int_trim)
head(merge_all@meta.data)

table(merge_all@meta.data$max_prediction, 
      merge_all@meta.data$orig.ident)

write.csv(table(merge_all@meta.data$max_prediction, 
                merge_all@meta.data$orig.ident), 
          "INTEGRATION_numCells_cellType_logNorm.csv")

#backup merge all 
backup_merge_all <- merge_all
#find differentially expressed genes
##Set identities to integration clusters
Idents(merge_all) <- merge_all@meta.data$max_prediction

new_merged_all <- JoinLayers(new_int) 
merge_all

#DefaultAssay(new_merged_all) <- "Spatial"
Idents(new_merged_all)
Idents(new_merged_all) <- new_merged_all@meta.data$broad_celltype
##Run differential expression test - infected vs uninfected
##Want to check these cell types: 13:CTB-2, 16:CTB-3, 15:T-cell (CD4 Naive), 32:T-cell, 
### 9:T-cell (CD4 TCM), 10:Decidual-1, 12:Decidual-2, 17:Decidual-3, 23:Endothelial,
### 18:STB, 4:Marophage-2 (Hofbauer), 6:Macrophage-1, 1:EVT-1, 24:EVT-2 14:LED

#de_markers_CTB2 <- FindMarkers(new_merged_all, ident.1 = "infected", group.by = "infection_state", subset.ident = "13:CTB-2")
#de_markers_CTB3 <- FindMarkers(new_merged_all, ident.1 = "infected", group.by = "infection_state", subset.ident = "16:CTB-3")
#de_markers_tCellCD4Naive <- FindMarkers(new_merged_all, ident.1 = "infected", group.by = "infection_state", subset.ident = "15:T-cell (CD4 Naive)")
#de_markers_Marophage2 <- FindMarkers(new_merged_all, ident.1 = "infected", group.by = "infection_state", subset.ident = "4:Marophage-2 (Hofbauer)")
#de_markers_tCellCD4TCM <- FindMarkers(new_merged_all, ident.1 = "infected", group.by = "infection_state", subset.ident = "9:T-cell (CD4 TCM)")

de_markers_CTB <- FindMarkers(new_merged_all, ident.1 = "infected", 
                              group.by = "infection_state", 
                              subset.ident = "CTB")

de_markers_Stromal <- FindMarkers(new_merged_all, ident.1 = "infected", 
                              group.by = "infection_state", 
                              subset.ident = "Stromal")

de_markers_tCell <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "T-cell")

de_markers_Decidual <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                    group.by = "infection_state",
                                    subset.ident = "Decidual")

de_markers_Endothelial <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "Endothelial")

de_markers_STB <- FindMarkers(new_merged_all, ident.1 = "infected", 
                              group.by = "infection_state",
                              subset.ident = "STB")

de_markers_Macrophage <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "Macrophage")

de_markers_EVT <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "EVT")

de_markers_LED <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                      group.by = "infection_state",
                                      subset.ident = "LED")

de_markers_bCell <- FindMarkers(new_merged_all, ident.1 = "infected", ##Error that one or both identiy groups not present in data
                              group.by = "infection_state",
                              subset.ident = "B-cell")

de_markers_npiCTB <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                group.by = "infection_state",
                                subset.ident = "npiCTB")
de_markers_npiCTB_mast <- FindMarkers(new_merged_all, ident.1 = "infected", 
                                 group.by = "infection_state",
                                 subset.ident = "npiCTB",
                                 test.use = "MAST")
de_markers_npiCTB[de_markers_npiCTB$p_val_adj<0.05,] ##no adj p-val less than 0.05
de_markers_CTB[de_markers_CTB$p_val_adj<0.05,] #no adj p-val less than 0.05
de_markers_npiCTB["ENSG00000169245",] #look for CXCL10
de_markers_npiCTB_mast["ENSG00000169245",] #look for CXCL10
de_markers_npiCTB["ENSG00000185507",] #look for IRF7 


VlnPlot(new_merged_all, features = "ENSG00000169245", split.by = "infection_state")


#list_de_markers <- list("de_markers_Endothelial", "de_markers_Decidual3", "de_markers_EVT1", "de_markers_LED", "de_markers_Macrophage1", "de_markers_STB","de_markers_tCell", "de_markers_tCellCD4TCM")
list_de_markers_updated <- list("de_markers_Endothelial", "de_markers_Decidual", "de_markers_EVT",
                        "de_markers_LED", "de_markers_Macrophage", "de_markers_STB",
                        "de_markers_tCell", "de_markers_Stromal") #, "de_markers_CTB"
head(de_markers_LED)
length(de_markers_LED[de_markers_LED$p_val_adj<0.05,])
match(rownames(de_markers_LED),conversion$Ensembl)

#make CSV files of each
for (i in list_de_markers) { #for file in list
  n <- get(sprintf("%s", i)) #get the gene table for each cluster and save it in n
  print(n)
  #write.csv(n, sprintf("diff_expression/%s.csv", i)) #write n to a csv file
}

write.csv(rownames(new_merged_all),"out_all_ensembl_ids_for_conversion.csv")
conversion <- read.csv("out_ensembl_geneSymbol_conversion.csv", header = TRUE)
conversion[1:5,1:2]


for (i in list_de_markers_updated) { #for file in list
  n <- get(sprintf("%s", i)) #get the gene table for each cluster and save it in n
  #print(n)
  n <- n[n$p_val_adj<0.05,]
  n$GeneSymbol <- conversion$GeneSymbol[match(rownames(n),conversion$Ensembl)]
  #print(i)
  #print(head(n))
  write.csv(n, sprintf("diff_expression_broad/with_stromal/%s.csv", i)) #write n to a csv file
}

for (i in list_de_markers_updated) { #for file in list
  n <- get(sprintf("%s", i)) #get the gene table for each cluster and save it in n
  #print(n)
  n <- n[n$p_val_adj<0.05,]
  n$GeneSymbol <- conversion$GeneSymbol[match(rownames(n),conversion$Ensembl)]
  #print(i)
  #print(head(n))
  print(n[n$GeneSymbol=="OAS3",])
  #write.csv(n, sprintf("diff_expression_broad/%s.csv", i)) #write n to a csv file
}




plots <- VlnPlot(new_merged_all, features = "ENSG00000213934" ,
                 group.by = "infection_state", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots)
write.csv(de_markers_0, "diff_expression/clust0_de.csv")
VlnPlot(merge_all, features = c("ISG15"),
        split.by = "orig.ident", pt.size = 0, combine = FALSE)




#requested genes
merge_all$infection_state <- factor(merge_all$infection_state, 
                                  levels = c("uninfected","infected"))
VlnPlot(merge_all, features = c("ENSG00000167207"), #NOD2
        split.by = "infection_state", pt.size = 0, combine = FALSE)
VlnPlot(new_int, features = c("ENSG00000106100"), #NOD1
        split.by = "infection_state", pt.size = 0, combine = FALSE)
VlnPlot(merge_all, features = c("ENSG00000169245"),
        split.by = "infection_state", pt.size = 0, combine = FALSE)
VlnPlot(merge_all, features = c("ENSG00000115594"),
        split.by = "infection_state", pt.size = 0, combine = FALSE)

VlnPlot(merge_all, features = c("ENSG00000089127","ENSG00000111335"),
        split.by = "infection_state", pt.size = 0, combine = FALSE)
VlnPlot(merge_all, features = c("ENSG00000111335"),
        split.by = "orig.ident", pt.size = 0, combine = FALSE)
VlnPlot(merge_all, features = c("ENSG00000146678"),
        split.by = "orig.ident", pt.size = 0, combine = FALSE)

Idents(new_int)
head(new_int@meta.data)
SpatialFeaturePlot(new_int, features = "ENSG00000106100",
                   images = "D1")
VlnPlot(new_int, features = c("ENSG00000167207"), #NOD2
        split.by = "infection_state", group.by = "max_prediction",
        pt.size = 0, combine = FALSE)
VlnPlot(new_int, features = c("ENSG00000106100"), #NOD1
        split.by = "infection_state", group.by = "max_prediction",
        pt.size = 0, combine = FALSE)
VlnPlot(new_int, features = c("ENSG00000146678"), group.by = "harmony_clusters",
        split.by = "infection_state", pt.size = 0, combine = FALSE)

innate_genes <- read.csv("innate_immunity_genes.csv")
innate_genes[1:5,1:2]
DoHeatmap(merge_all, features = innate_genes$converted_alias, 
          group.by = "infection_state", group.bar = TRUE, combine = TRUE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = innate_genes$converted_alias,
    labels = innate_genes$initial_alias)
ggsave("heatmap_innate_genes_logNorm.png",
       height = 8, width = 16, 
       dpi = 300, device = "png")

DoHeatmap(new_int, features = innate_genes$converted_alias, 
          group.by = "infection_state", group.bar = TRUE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = innate_genes$converted_alias,
    labels = innate_genes$initial_alias)

#Get genes of interest from Ravit's files
goi <- read.csv("desired_violin_ravit/final_unique_desired_violin.csv")
goi_unique <- goi$Unique[1:88]
goi
head(goi)

for (i in goi$Unique) { 
  print(i)
  p <- VlnPlot(merge_all, features = as.character(i),
               split.by = "orig.ident", group.by = "infection_state", 
               pt.size = 0, combine = TRUE)
  ggsave(sprintf("desired_violin_ravit/figs/gene_%s_violinplot.png", i), p, dpi = 300, 
         device = "png",width = 10, height = 10)
}

##Test renaming the plot
p <- VlnPlot(merge_all, features = "ENSG00000183486",
        split.by = "infection_state", group.by = "max_prediction", pt.size = 0, 
        combine = TRUE)
p <- p + ggtitle(goi$UniqueGeneID[goi$Unique == "ENSG00000183486"])
p
ggsave("desired_violin_ravit/figs/gene_ENSG00000183486_violinplot.png", p, 
              dpi = 300, device = "png",width = 15, height = 10)
##Run generation of figures
for (i in goi_unique) { #for clusters
  print(i)
  #print(goi$UniqueGeneID[goi$Unique == i])
  f <- goi$UniqueGeneID[goi$Unique == i]
  print(f)
  p <- VlnPlot(merge_all, features = i,
               split.by = "infection_state", group.by = "broad_celltype", pt.size = 0, 
               combine = TRUE)
  p <- p + ggtitle(goi$UniqueGeneID[goi$Unique == i])
  ggsave(sprintf("desired_violin_ravit/figs/by_cluster/gene_%s_violinplot.png", f), p, 
         dpi = 300, device = "png",width = 15, height = 10)
}

##Run generation of figures
for (i in goi_unique) { #for infection state
  print(i)
  #print(goi$UniqueGeneID[goi$Unique == i])
  f <- goi$UniqueGeneID[goi$Unique == i]
  print(f)
  p <- VlnPlot(merge_all, features = i,
               group.by = "infection_state", pt.size = 0, 
               combine = TRUE)
  p <- p + ggtitle(goi$UniqueGeneID[goi$Unique == i])
  ggsave(sprintf("desired_violin_ravit/figs/by_infection/gene_%s_violinplot.png", f), p, 
         dpi = 300, device = "png",width = 15, height = 10)
}

p1 <- VlnPlot(merge_all, features = c("KISS1"),
              split.by = "orig.ident", group.by = "infection_state", pt.size = 0, combine = FALSE)
p2 <- VlnPlot(merge_all, features = c("ISG15"),
              split.by = "orig.ident", pt.size = 0, combine = FALSE)
wrap_plots(p1,p2)

saveRDS(merge_all, "RDS_merge_all.rds")


##Relabel cell types to have one per name (no subtypes)
head(merge_all@meta.data)
merge_all[["broad_celltype"]] <- merge_all@meta.data$max_prediction
head(merge_all@meta.data)
merge_all@meta.data$broad_celltype <- gsub(".*:","",merge_all@meta.data$broad_celltype)
merge_all@meta.data$broad_celltype <- gsub("-[0-9]","",merge_all@meta.data$broad_celltype)
merge_all@meta.data$broad_celltype <- gsub(" .*","",merge_all@meta.data$broad_celltype)
merge_all@meta.data$broad_celltype <- gsub("Marophage","Macrophage",merge_all@meta.data$broad_celltype)
unique(merge_all@meta.data$max_prediction)
unique(merge_all@meta.data$broad_celltype)

DimPlot(merge_all, reduction = "umap", group.by = "broad_celltype",
        combine = FALSE, label.size = 2, split.by = "infection_state")


head(new_int@meta.data)


new_int[["broad_celltype"]] <- new_int@meta.data$max_prediction
head(new_int@meta.data)
new_int@meta.data$broad_celltype <- gsub(".*:","",new_int@meta.data$broad_celltype)
new_int@meta.data$broad_celltype <- gsub("-[0-9]","",new_int@meta.data$broad_celltype)
new_int@meta.data$broad_celltype <- gsub(" .*","",new_int@meta.data$broad_celltype)
new_int@meta.data$broad_celltype <- gsub("Marophage","Macrophage",new_int@meta.data$broad_celltype)
unique(new_int@meta.data$max_prediction)
unique(new_int@meta.data$broad_celltype)
#look at new umap/harmony plots
DimPlot(new_int, reduction = "umap.harmony", label = TRUE, 
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", label = TRUE, 
        group.by = "infection_state")
DimPlot(new_int, reduction = "umap.harmony", 
        label = TRUE, group.by = "broad_celltype",
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", 
        label = TRUE, split.by = "infection_state", 
        group.by = "broad_celltype",
        repel = TRUE)
DimPlot(new_int, reduction = "umap.harmony", 
        label = FALSE, split.by = "infection_state", 
        group.by = "broad_celltype",
        repel = TRUE)

dim_ct <- DimPlot(new_int, reduction = "umap.harmony", 
                  label = FALSE, split.by = "infection_state", 
                  group.by = "broad_celltype",
                  repel = TRUE)
dim_ct <- dim_ct + theme(
  axis.text = element_text(size = 16),
  plot.title = element_text(size = 20)) +
  ggtitle("Broad cell types")
LabelClusters(dim_ct, id = "broad_celltype", fontface = "bold")
ggsave(filename = "fig_dimplot_broad-celltypes_infState_bold.png", height = 7, width = 12, 
       dpi = 300, device = "png")

dim_ct_toget <- DimPlot(new_int, reduction = "umap.harmony", 
                  label = FALSE, 
                  group.by = "broad_celltype",
                  repel = TRUE)
dim_ct_toget <- dim_ct_toget + theme(
  axis.text = element_text(size = 16),
  plot.title = element_text(size = 20)) +
  ggtitle("Broad cell types")

LabelClusters(dim_ct_toget, id = "broad_celltype", fontface = "bold")
ggsave(filename = "fig_dimplot_broad-celltypes_bold.png", height = 7, width = 12, 
       dpi = 300, device = "png")

#Get genes of interest from Ravit's files
goi_pathways <- read.csv("GSEA/STRING/selected_plus_NOD2/for_NOD2_violin_plots.csv")
goi_pathways_NODRec <- goi_pathways[,1:2]
head(goi_pathways_NODRec)
length(goi_pathways_NODRec$Ensembl)
length(unique(goi_pathways_NODRec$Ensembl))

goi_pathways_IBD <- goi_pathways[1:8,4:5]
head(goi_pathways_IBD)
length(goi_pathways_IBD$Ensembl)
length(unique(goi_pathways_IBD$Ensembl))

goi_pathways_TNF <- goi_pathways[1:9,7:8]
head(goi_pathways_TNF)
length(goi_pathways_TNF$Ensembl)
length(unique(goi_pathways_TNF$Ensembl))

for (i in goi_pathways_NODRec$Ensembl) { #for clusters
  print(i)
  print(goi_pathways_NODRec$Gene[goi_pathways_NODRec$Ensembl == i])
  g <- goi_pathways_NODRec$Gene[goi_pathways_NODRec$Ensembl == i]
  p <- VlnPlot(merge_all, features = i,
               split.by = "infection_state", group.by = "broad_celltype", pt.size = 0, 
               combine = TRUE)
  p <- p + ggtitle(goi_pathways_NODRec$Gene[goi_pathways_NODRec$Ensembl == i])
  ggsave(sprintf("GSEA/STRING/selected_plus_NOD2/violin_plots/NOD_receptor/gene_%s_violinplot.png", g), p, 
         dpi = 300, device = "png",width = 15, height = 10)
}

for (i in goi_pathways_IBD$Ensembl) { #for clusters
  print(i)
  print(goi_pathways_IBD$Gene[goi_pathways_IBD$Ensembl == i])
  g <- goi_pathways_IBD$Gene[goi_pathways_IBD$Ensembl == i]
  p <- VlnPlot(merge_all, features = i,
               split.by = "infection_state", group.by = "broad_celltype", pt.size = 0, 
               combine = TRUE)
  p <- p + ggtitle(goi_pathways_IBD$Gene[goi_pathways_IBD$Ensembl == i])
  ggsave(sprintf("GSEA/STRING/selected_plus_NOD2/violin_plots/IBS/gene_%s_violinplot.png", g), p, 
         dpi = 300, device = "png",width = 15, height = 10)
}

for (i in goi_pathways_TNF$Ensembl) { #for clusters
  print(i)
  print(goi_pathways_TNF$Gene[goi_pathways_TNF$Ensembl == i])
  g <- goi_pathways_TNF$Gene[goi_pathways_TNF$Ensembl == i]
  p <- VlnPlot(merge_all, features = i,
               split.by = "infection_state", group.by = "broad_celltype", pt.size = 0, 
               combine = TRUE)
  p <- p + ggtitle(goi_pathways_TNF$Gene[goi_pathways_TNF$Ensembl == i])
  ggsave(sprintf("GSEA/STRING/selected_plus_NOD2/violin_plots/TNF/gene_%s_violinplot.png", g), p, 
         dpi = 300, device = "png",width = 15, height = 10)
}


##Testing dotplot for NOD-like receptor
##Switch to gene names instead of ENSEMBL: https://www.biostars.org/p/9596576/
Idents(merge_all) <- merge_all@meta.data$broad_celltype
gene_labels <- goi_pathways_NODRec$Gene
p1 <- DotPlot(merge_all, features = goi_pathways_NODRec$Ensembl,
        split.by = "infection_state") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1.5 <- DotPlot(merge_all, features = goi_pathways_NODRec$Ensembl, 
                split.by = "infection_state")

p1.5$data$features.plot <- factor(p1.5$data$features.plot,
                                  levels = goi_pathways_NODRec$Ensembl,
                                  labels = gene_labels)

p2 <- p1.5 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

merge_all$infect_celltype <- paste0(Idents(merge_all),"_",merge_all$infection_state)
DotPlot(merge_all, features = goi_pathways_NODRec$Ensembl,
        group.by = "infect_celltype") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 <- DotPlot(merge_all, features = goi_pathways_NODRec$Ensembl,
              group.by = "infect_celltype") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1.5 <- DotPlot(merge_all, features = goi_pathways_NODRec$Ensembl,
                group.by = "infect_celltype")

p1.5$data$features.plot <- factor(p1.5$data$features.plot,
                                  levels = goi_pathways_NODRec$Ensembl,
                                  labels = gene_labels)

p2 <- p1.5 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- p2 + theme(plot.background = element_rect(fill = "white"))

ggsave(filename = "GSEA/STRING/selected_plus_NOD2/fig_dotplot_broad-celltypes.png", 
       p2, height = 7, width = 12, 
       dpi = 300, device = "png")


write.csv(table(merge_all@meta.data$broad_celltype, 
                merge_all@meta.data$infection_state), 
          "INTEGRATION_numCells_broad-cellType_logNorm.csv")

head(new_int@meta.data)
head(rownames(new_int))
saveRDS(new_int, "RDS_new_int.rds")


new_int <- JoinLayers(new_int)
cxcl10A1 <- SpatialFeaturePlot(new_int, features = "ENSG00000169245", images = "A1", #CXCL10
                   pt.size.factor = 5, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10A1.2 <- SpatialFeaturePlot(new_int, features = "ENSG00000169245", images = "A1.2", #CXCL10
                               pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10D1 <- SpatialFeaturePlot(new_int, features = "ENSG00000169245", images = "D1", #CXCL10
                                 pt.size.factor = 2, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
cxcl10D1.4 <- SpatialFeaturePlot(new_int, features = "ENSG00000169245", images = "D1.4", #CXCL10
                               pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "CXCL10") 
wrap_plots(cxcl10A1, cxcl10A1.2, cxcl10D1, cxcl10D1.4)

ggsave("spatail_feature_plots/CXCL10_spatialFeatures.svg", device = "svg",
       height = 10, width = 10)


isg15A1 <- SpatialFeaturePlot(new_int, features = "ENSG00000187608", images = "A1", #isg15
                               pt.size.factor = 5, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15A1.2 <- SpatialFeaturePlot(new_int, features = "ENSG00000187608", images = "A1.2", #isg15
                                 pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15D1 <- SpatialFeaturePlot(new_int, features = "ENSG00000187608", images = "D1", #isg15
                               pt.size.factor = 2, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
isg15D1.4 <- SpatialFeaturePlot(new_int, features = "ENSG00000187608", images = "D1.4", #isg15
                                 pt.size.factor = 4, combine = TRUE) +
  theme(legend.position = "right") +
  labs(fill = "ISG15") 
wrap_plots(isg15A1, isg15A1.2, isg15D1, isg15D1.4)

ggsave("spatail_feature_plots/ISG15_spatialFeatures.svg", device = "svg",
       height = 10, width = 10)


"#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499","#44AA99","#999933","#882255",
"#661100","#6699CC","#888888"

"#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90",
"#f97b72","#4b4b8f","#A5AA99"

head(new_int@meta.data)
Idents(new_int) <- new_int$broad_celltype
new_int$broad_celltype <- factor(new_int$broad_celltype, 
                            levels = c("B-cell","CTB", "Decidual","Endothelial", "EVT", "LED",
                                       "Macrophage", "npiCTB", "STB", "Stromal", "T-cell"))
my_groups <- levels(new_int$broad_celltype)
colors <- c("#7F3C8D","#11A579","#3969AC","#A5AA99","#E73F74",
            "#80BA5A","#E68310","#008695","#661100","#f97b72","#F2B701")
names(colors) <- my_groups

ct.A1 <- SpatialDimPlot(new_int, group.by = "broad_celltype", pt.size.factor = 5, images = "A1",
               combine = T,
               cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.A1.2 <- SpatialDimPlot(new_int, group.by = "broad_celltype", pt.size.factor = 4, images = "A1.2",
                        combine = T,
                        cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.D1 <- SpatialDimPlot(new_int, group.by = "broad_celltype", pt.size.factor = 2, images = "D1",
                        combine = T,
                        cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ct.D1.4 <- SpatialDimPlot(new_int, group.by = "broad_celltype", pt.size.factor = 4, images = "D1.4",
                        combine = T,
                        cols = colors) +
  labs(fill = "Cell Type") +
  guides(fill = guide_legend(override.aes = list(size = 5)))

wrap_plots(ct.A1, ct.A1.2, ct.D1, ct.D1.4)

ggsave("spatail_feature_plots/celltypes_spatialFeatures_recolored_v2.png", device = "png",
       height = 10, width = 15)
ggsave("spatail_feature_plots/celltypes_spatialFeatures_recolored_v2.svg", device = "svg",
       height = 10, width = 15)

Idents(new_int) <- new_int@meta.data$broad_celltype
dim_ct_recolor <- DimPlot(new_int, reduction = "umap.harmony", 
                  cols = colors,
                  label = FALSE, split.by = "infection_state", 
                  group.by = "broad_celltype",
                  repel = TRUE)
dim_ct_recolor <- dim_ct_recolor + theme(
  axis.text = element_text(size = 18),
  plot.title = element_text(size = 20)) +
  ggtitle("Broad cell types")
LabelClusters(dim_ct_recolor, id = "broad_celltype", fontface = "bold")
ggsave(filename = "fig5_dimplot_broad-celltypes_infState_bold_recolor_v2.png", height = 7, 
       width = 12, dpi = 300, device = "png")
ggsave(filename = "fig5_dimplot_broad-celltypes_infState_bold_recolor_v2.tiff", height = 7, 
       width = 12, dpi = 300, device = "tiff")
ggsave(filename = "fig5_dimplot_broad-celltypes_infState_bold_recolor_v2.pdf", height = 7, 
       width = 12, dpi = 300, device = "pdf")



ISG.for.heat <- read.csv("desired_violin_ravit/final_unique_desired_violin.csv")
dim(ISG.for.heat)
ISG.for.heat[1:5,1:2]
innate_genes[1:5,1:2]
DoHeatmap(new_int, features = ISG.for.heat$Unique, 
          group.by = "infection_state", group.bar = TRUE, combine = TRUE,
          label = TRUE, angle = 0) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = ISG.for.heat$Unique,
    labels = ISG.for.heat$UniqueGeneID) +
  theme(axis.text.y = element_text(size = 18),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) 
ggsave("heatmap_ISG_genes_logNorm_infection-state.png",
       height = 18, width = 10, 
       dpi = 300, device = "png")
ggsave("heatmap_ISG_genes_logNorm_infection-state.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("heatmap_ISG_genes_logNorm_infection-state.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

ggsave("heatmap_ISG_genes_logNorm_infection-state_wLabels.png",
       height = 18, width = 10, 
       dpi = 300, device = "png")
ggsave("heatmap_ISG_genes_logNorm_infection-state_wLabels.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("heatmap_ISG_genes_logNorm_infection-state_wLabels.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

DoHeatmap(new_int, features = ISG.for.heat$Unique, 
          group.by = "broad_celltype", group.bar = TRUE, group.colors = colors,
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
ggsave("heatmap_ISG_genes_logNorm_celltypes.png",
       height = 18, width = 10, 
       dpi = 300, device = "png")
ggsave("heatmap_ISG_genes_logNorm_celltypes.tiff",
       height = 18, width = 10, 
       dpi = 300, device = "tiff")
ggsave("heatmap_ISG_genes_logNorm_celltypes.pdf",
       height = 18, width = 10, 
       dpi = 300, device = "pdf")

ggsave("heatmap_ISG_genes_logNorm_celltypes_wLabels_v4.png",
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

DoHeatmap(new_int, features = innate_genes$converted_alias, 
          group.by = "infection_state", group.bar = TRUE, combine = TRUE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = innate_genes$converted_alias,
    labels = innate_genes$initial_alias)
ggsave("heatmap_innate_genes_logNorm_infection-state.png",
       height = 8, width = 16, 
       dpi = 300, device = "png")

DoHeatmap(new_int, features = innate_genes$converted_alias, 
          group.by = "infection_state", group.bar = TRUE, label = FALSE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = innate_genes$converted_alias,
    labels = innate_genes$initial_alias) +
  theme(axis.text.y = element_text(size = 18),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5)))
ggsave("heatmap_innate_genes_logNorm_infection-state.png",
       height = 16, width = 16, 
       dpi = 300, device = "png")

DoHeatmap(new_int, features = innate_genes$converted_alias, 
          group.by = "broad_celltype", group.bar = TRUE, group.colors = colors,
          label = FALSE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = innate_genes$converted_alias,
    labels = innate_genes$initial_alias) +
  theme(axis.text.y = element_text(size = 18),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5)))
ggsave("heatmap_innate_genes_logNorm_celltype.png",
       height = 16, width = 16, 
       dpi = 300, device = "png")


comb.ISG.innate <- read.csv("combined_ISG_innate.csv")
dim(comb.ISG.innate)
comb.ISG.innate[1:5,1:4]
DoHeatmap(new_int, features = comb.ISG.innate$converted_alias, 
          group.by = "infection_state", group.bar = TRUE, combine = TRUE,
          label = FALSE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = comb.ISG.innate$converted_alias,
    labels = comb.ISG.innate$initial_alias) +
  theme(axis.text.y = element_text(size = 20),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) 
ggsave("heatmap_ISG-innate_genes_logNorm_infection-state.png",
       height = 26, width = 10, 
       dpi = 300, device = "png")

DoHeatmap(new_int, features = comb.ISG.innate$converted_alias, 
          group.by = "broad_celltype", group.bar = TRUE, group.colors = colors,
          combine = TRUE, label = FALSE) +
  scale_y_discrete( #https://stackoverflow.com/questions/65214766/label-y-axis-with-a-different-column-in-ggplot2
    breaks = comb.ISG.innate$converted_alias,
    labels = comb.ISG.innate$initial_alias) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme(axis.text.y = element_text(size = 20),
        legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) 
ggsave("heatmap_ISG-innate_genes_logNorm_celltypes.png",
       height = 26, width = 10, 
       dpi = 300, device = "png")

table(new_int@meta.data$broad_celltype,
      new_int@meta.data$infection_state)
write.csv(table(new_int@meta.data$broad_celltype,
                new_int@meta.data$infection_state), "out_table_celltypes_infection.csv")

DefaultAssay(new_int)
Idents(new_int)
head(new_int@meta.data)
Idents(new_int) <- new_int$broad_celltype
DotPlot(new_int, features = c("ENSG00000182393","ENSG00000183709","ENSG00000185436","ENSG00000243646")) +
  scale_x_discrete(breaks = c("ENSG00000182393","ENSG00000183709","ENSG00000185436","ENSG00000243646"),
                   labels = c("IFNL1","IFNL2","IFNLR1","IL10RB"))

VlnPlot(new_int, features = c("ENSG00000182393","ENSG00000183709","ENSG00000185436","ENSG00000243646"), 
        group.by = "broad_celltype")

Idents(new_int) <- new_int$broad_celltype
DotPlot(new_int, features = c("ENSG00000135346", "ENSG00000116183","ENSG00000129226", 
                              "ENSG00000010610", #"ENSG00000168685",
                              "ENSG00000110799")) +
  scale_x_discrete(breaks = c("ENSG00000135346", "ENSG00000116183","ENSG00000129226", 
                              "ENSG00000010610", #"ENSG00000168685",
                              "ENSG00000110799"),
                   labels = c("CGA","PAPPA2", "CD68", 
                              "CD4", #"IL7R", 
                              "VWF"))
ggsave("heatmap_ISG-innate_genes_logNorm_celltypes.png",
       height = 26, width = 10, 
       dpi = 300, device = "png")


length(WhichCells(new_int, idents = "Macrophage", expression = ENSG00000010610 > 0))
length(WhichCells(new_int, idents = "Macrophage"))

DotPlot(new_int, features = c("ENSG00000213949","ENSG00000161638","ENSG00000116183","ENSG00000087245")) +
  scale_x_discrete(breaks = c("ENSG00000213949","ENSG00000161638","ENSG00000116183","ENSG00000087245"),
                   labels = c("ITGA1","ITGA5","PAPPA2","MMP2"))
