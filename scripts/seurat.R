library(Seurat)

pbmc <- CreateSeuratObject(counts = mz_0_TIC, project = "lyso", 
                           min.cells = 0, min.features = 0,
                           meta.data = data_cell)
pbmc

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = all.genes, npcs = 20, nfeatures.print = 3)

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca", group.by = "plate")
DimPlot(pbmc, reduction = "pca", group.by = "group")

pbmc <- JackStraw(pbmc, num.replicate = 200)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

######################### Cluster the cells ###########################
# takes as input the previously defined dimensionality of the dataset (first 10 PCs)
pbmc <- FindNeighbors(pbmc, dims = dims, k.param = k.param, reduction = "pca")
pbmc <- FindClusters(pbmc, resolution = resolution)
head(Idents(pbmc), 5)

pbmc <- RunTSNE(pbmc, dims = dims, perplexity = perplexity, seed.use = seed.use)
DimPlot(pbmc, reduction = "tsne", label = T, pt.size = 1.8)
DimPlot(pbmc, reduction = "tsne", group.by = "plate", pt.size = 1.8)
DimPlot(pbmc, reduction = "tsne", group.by = "group", pt.size = 1.8)
