# monocle -----------------------------------------------------------------
library(monocle)
library(dplyr)

feature_meta <- data.frame(gene_short_name = rownames(mz_0_TIC), 
                           row.names = rownames(mz_0_TIC))
pd <- new("AnnotatedDataFrame", data = data_cell)
fd <- new("AnnotatedDataFrame", data = feature_meta)

# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(mz_0_TIC),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0,
                       expressionFamily = tobit())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Clustering cells without marker genes -----------------------------------
HSMM <- setOrderingFilter(HSMM, feature_meta$gene_short_name)

plot_pc_variance_explained(HSMM, return_all = F, max_components = 50)

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = num_dim, 
                        reduction_method = 'tSNE', verbose = T, perplexity = perplexity)

HSMM <- clusterCells(HSMM, num_clusters = 5)
table(HSMM@phenoData@data$group, HSMM@phenoData@data$Cluster)

plot_cell_clusters(HSMM, 1, 2, color = "group", cell_size = 2) +
  scale_color_hue(l = 55) + 
  scale_size(range = c(1, 2)) + 
  labs(x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic()

ggsave(filename = "tSNE.tiff", device = "tiff", path = "./")
ggsave(filename = "tSNE.eps", device = "eps", path = "./")

# Constructing Single Cell Trajectories -----------------------------------
ordering_genes <- unique(feature_meta$gene_short_name)
HSMM <- setOrderingFilter(HSMM, ordering_genes)

HSMM <- reduceDimension(HSMM, max_components = 2, reduction_method = "DDRTree")
HSMM <- orderCells(HSMM)

plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "group")

