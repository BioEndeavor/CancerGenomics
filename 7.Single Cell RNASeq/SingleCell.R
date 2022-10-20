#Import Libraries
library(Seurat)
library(tidyr)
library(dplyr)

#Read the dataset that includes barcodes,features and matrix file
pbmc.data = Read10X(data.dir = "C:/Users/joash/OneDrive/Desktop/pbmc/leukemia/")

#Create Seurat Object
pbmc = CreateSeuratObject(counts = pbmc.data,min.cells = 3,min.features = 200)
pbmc
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

#Plot ViolinPlot
pdf('Violin.pdf')
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#Plot FeatureScatter
pdf('Feature.pdf')
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc
seu <- FindVariableFeatures(object = pbmc)

#Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = seu))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#Plot DimensionPlot
pdf('DimPlot.pdf')
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#ElbowPlot(pbmc)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

pbmc = RunUMAP(pbmc, dims = 1:10)

pdf('Umap.pdf')
DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, reduction = "umap", label = T)
dev.off()

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

a = pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
genes = a %>% pull(gene)

#Plot FeaturePlot
pdf('FeaturePlot.pdf')
FeaturePlot(pbmc, features = genes[1:2])
FeaturePlot(pbmc, features = genes[1:2], cols = c("white", "red"))
dev.off()