# filter ------------------------------------------------------------------
rm(list = ls())
data_raw <- read.csv("res_raw_0.csv", row.names = 2)
data_raw <- data_raw[, c(-1:-8)]
data_raw[is.na(data_raw)] <- 1

data_signal <- data_raw[, seq(1, ncol(data_raw)-1, 2)]
data_noise <- data_raw[, seq(2, ncol(data_raw), 2)]
data_S_N <- data_signal / data_noise

data_raw <- read.csv("res_raw_0.csv", row.names = 2)
data_raw <- data_raw[, c(-1:-8)]
data_raw <- data_raw[apply(data_S_N, 1, 
                           function(x) {sum(x >= 3, na.rm = T)}) > ncol(data_raw)/2*0.2, 
                     seq(1, ncol(data_raw)-1, 2)]

library(stringr)
grp <- rep("lysosome", ncol(data_raw))
grp[grep(colnames(data_raw), pattern = "[Ee][Nn][Dd][Oo]")] <- "endosome"
grp[grep(colnames(data_raw), pattern = "[Aa][Uu][Tt]")] <- "autolysosome"
data_cell <- data.frame(group = grp, 
                        plate = str_split(colnames(data_raw), "\\.", simplify = T)[,1],
                        sample = str_split(colnames(data_raw), "\\.", simplify = T)[,2])
rownames(data_cell) <- colnames(data_raw)

# impute ------------------------------------------------------------------
my_impute <- function(x, method = "0"){
  if(method == "0"){
    x[is.na(x)] <- 0
  }
  return(x)
}

mz_0 <- my_impute(data_raw, method = "0")

# PCA ---------------------------------------------------------------------
library("stringr")
library("FactoMineR")
library("factoextra") 

test <- "group"
dat <- cbind(t(mz_0), data_cell[test])
dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca, repel =T,
             geom.ind = c("point"), # show points only (nbut not "text")
             col.ind = dat[[test]], # color by groups
             mean.point = F,
             addEllipses = F, # Concentration ellipses
             legend.title = "Groups")

# FC, t test ---------------------------------------------------------------
i <- group
df <- mz_0
FC <- apply(df, 1, function(x){
  mean(x[data_cell$group == i]) / mean(x[data_cell$group != i])
})
p.value <- apply(df, 1, function(x){
  t.test(x[data_cell$group == i], x[data_cell$group != i], 
              alternative = "two.sided", paired = FALSE)$p.value
})
marker <- data.frame(FC = FC, p.value = p.value)
rownames(marker) <- rownames(df)
marker$change <- as.factor(ifelse(marker$p.value < 0.05 & abs(log2(marker$FC)) > log2(1.2), 
                                  ifelse(log2(marker$FC) > log2(1.2), "UP", "DOWN"), "NOT"))
marker$mz <- rownames(marker)

# volcano plot ---------------------------------------------------------------
library(ggplot2)

ggplot(data = marker, aes(x = log2(FC), y = -log10(p.value), color = change)) + 
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 15))) + 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme(plot.title = element_text(size = 15, hjust = 0.5)) + 
  scale_colour_manual(values = c("blue", "black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size=1)

# plot heatmap--------------------------------------------------------------------

marker_one <- subset(marker, p.value < 0.05 & log2(FC) > log2(1.2))
marker_two <- subset(marker, p.value < 0.05 & log2(FC) < -log2(1.2))

library(pheatmap)
plot_map <- rbind(mz_0[rownames(marker_one), ], 
                  mz_0[rownames(marker_two), ])

n <- t(scale(t(plot_map)))
n[n > 1.5] = 1.5
n[n < -1.5] = -1.5

pheatmap(n[, order(data_cell$group)], 
         show_colnames = F, show_rownames = F, 
         cluster_cols = F, cluster_rows = T, 
	 annotation_col = data_cell["group"])

