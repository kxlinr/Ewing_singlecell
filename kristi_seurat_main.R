df_seurat=readRDS("C:/Users/KLin/Desktop/seur_obj.RDS")

#Pre-Processing Seurat Workflow
df_seurat[["percent.mt"]] <- PercentageFeatureSet(df_seurat, pattern = "^MT-") #Creating a percent.mt column for showing percentage of mitochondrial genes
VlnPlot(df_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #Visualize QC metrics as a violin plot
df_seurat <- subset(df_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #subsetting the dataset and removing bad quality cells

#Normalizing the dataset
df_seurat <- NormalizeData(df_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
df_seurat <- FindVariableFeatures(df_seurat, selection.method = "vst", nfeatures = 2000)

#Integration of Datasets
#df_seurat.anchors <- FindIntegrationAnchors(object.list = df_seurat.list, anchor.features = features)

#Scaling the data
all.genes <- rownames(df_seurat)
df_seurat <- ScaleData(df_seurat, features = all.genes)

#Perform linear dimensional reduction (PCA) and visualize using DimHeatmap
df_seurat <- RunPCA(df_seurat, features = VariableFeatures(object = df_seurat))
DimHeatmap(df_seurat, dims = 1:15, cells = 500, balanced = TRUE)

#Cluster Cells
df_seurat <- FindNeighbors(df_seurat, dims = 1:10)
df_seurat <- FindClusters(df_seurat, resolution = 0.5)

#Run non-linear dimensionality reduction
df_seurat <- RunUMAP(df_seurat, dims = 1:10)

#Finding differentially expressed features (cluster biomarkers)
df_seurat.markers <- FindAllMarkers(df_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#--------------------------------------------------------------
#Find the top 10 genes expressed in each cluster
top10 <- df_seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(df_seurat, features = top10$gene) + NoLegend()

#Assigning cell type identity to clusters and plotting UMAP
##Manually change "features" to be the gene markers of each cell type (or whatever gene we're interested in)
features <- c(#"EWSR1") #just because
  
              ##no osteocyte or osteoclast markers found
              #"MCAM", "COL1A1", "COL1A2", "SPARC") #osteoblast markers
              #"CD44", "CD151", "SOX5", "SOX9", "COL1A1", "COL1A2", "COL4A1", "COL4A2") #chondrocytes
  
              #"BMPR2", "CD44", "CD45", "MCAM", "LY6E", "NGFRAP1") #MSC markers
              #"NES", "PAX6", "TUBB3") #NEC
              ## no hematapoietic stem cell markers found

              "GFAP", "VIM", "DCX", "MKI67") #radial glial markers

              #"NES", "CD44") #Key markers (?)


         ##Try to plot metadata -- none of these work lol
              #df_seurat@meta.data)
              #df_seurat@meta.data[["orig.ident"]])
              #"ES_origin")

#RidgePlot(object = df_seurat, features = features, ncol = 2) #didnt use this bc the cluster axes overlap and hide e/o
#VlnPlot(object = df_seurat, features = features)
FeaturePlot(object=df_seurat, features=features)


########How to highlight sample source using metadata?
#-----------------------------------------------------------------
new.cluster.ids <- c("EWS", "MSC", "NEC")
names(new.cluster.ids) <- levels(df_seurat)
df_seurat <- RenameIdents(df_seurat, new.cluster.ids)
DimPlot(df_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(df_seurat, file = "X.rds")
