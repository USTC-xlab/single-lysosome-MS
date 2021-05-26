library(ggplot2)

for(i in 0:4){
  assign(paste("cluster", i, sep = ""), 
         data.frame(mz = as.numeric(rownames(mz_0_TIC[, data_cell$seurat_clusters == i])), 
                    intensity = apply(mz_0_TIC[, data_cell$seurat_clusters == i], 1, mean)))
  
  ggplot(get(paste("cluster", i, sep = ""))) + 
    geom_linerange(aes(x = mz, ymin = 0, ymax = intensity, color = "red")) + 
    labs(x = "m/z", y = "intensity", title = paste("cluster", i, sep = "")) +
    scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0)) +
    scale_x_continuous(limits = c(50, 750), breaks = c(seq(50, 750, 100)), expand = c(0.01, 0)) +
    theme_classic()
  
  ggsave(filename = paste("cluster", i, ".eps", sep = ""), device = "eps", path = "./")
}
