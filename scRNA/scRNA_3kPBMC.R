library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(pbmc$nFeature_RNA,probs = c(0.01,0.99))
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 600 & nFeature_RNA < 4000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
VariableFeatures =VariableFeatures(pbmc)
write.csv(VariableFeatures,'VariableFeatures.csv')

pbmc <- ScaleData(pbmc, features = VariableFeatures)
pbmc <- RunPCA(pbmc, features = VariableFeatures)

pbmc <- FindNeighbors(pbmc, dims = 1:8)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:8)
DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "pbmc_cluster.rds")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(pbmc.markers,'pbmc.markers.csv')

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

RNA_Count=pbmc[VariableFeatures,]@assays$RNA@layers$counts
rownames(RNA_Count)<-VariableFeatures
cell_ids <- Cells(pbmc)
colnames(RNA_Count)<-cell_ids

write.csv(RNA_Count,'RNA_Count.csv')
write.csv(pbmc@meta.data,'cluster.csv')

Naive_CD4_T=c('IL7R', 'CCR7')
CD14_Mono =c('CD14', 'LYZ')
Memory_CD4=c('IL7R', 'S100A4')
B='MS4A1'
CD8_T = 'CD8A'
FCGR3A_Mono =c('FCGR3A', 'MS4A7')
NK=c("GNLY", "NKG7")
DC=c('FCER1A', 'CST3')
Platelet='PPBP'

for (i in 0:11) {
  print(paste0(i,':',intersect(top20[top20$cluster==i,'gene']$gene,DC)))
}
#____top10_____
# 0-Naive_CD4_T->"IL7R" ;Memory_CD4->"IL7R"
# 2:CD8_T->"CD8A"; 
# 3-Naive_CD4_T->"CCR7"
# 4-CD14_Mono->"LYZ"
# 5-B->"MS4A1"
# 6-NK->"NKG7"
# 8-FCGR3A_Mono->"FCGR3A",'MS4A7'
# 9-DC->"FCER1A"

#________top20___
# Naive_CD4_T:0,2,3
# CD14_Mono :4
# Memory_CD4:0
# B:5
# NK:6,7
# DC:4,9

# FCGR3A_Mono:1->"MS4A7"
