#load Monocel3 object cds
Track_genes <- graph_test(cds,
                          neighbor_graph = "principal_graph",
                          cores = 6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes,"Trajectory_genes.csv",row.names = F)



library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

genes <- row.names(subset(Track_genes,q_value < 0.01 & morans_I > 0.2))
plot_matrix <- exprs(cds)[match(genes,
                                rownames(rowData(cds))),
                          order(pseudotime(cds))]

#smooth.spline,Z-score
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes
dim(plot_matrix)

plot_matrix_combin <- list()
for (i in 1:length(seq(1, 25353-53, 100))){
  num <- seq(1, 25353-53, 100)
  A <- plot_matrix[,num[i]:(100+num[i]-1)]
  a <- rowMeans(A)
  a <- as.data.frame(a)
  a <- a$a
  plot_matrix_combin[[i]] <- a
}

length(plot_matrix_combin)
tailed <- as.data.frame(rowMeans(plot_matrix[,25301:25353]))
plot_matrix_combin[[254]] <- tailed[,1]

plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
rownames(plot_matrix_combin) <- rownames(plot_matrix)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap::pheatmap(plot_matrix_combin, 
                   useRaster = T,
                   cluster_cols=FALSE, 
                   cluster_rows=T, 
                   show_rownames=T, 
                   show_colnames=F, 
                   clustering_method = "ward.D2",
                   cutree_rows=7,
                   filename="pseudotime_heatmap_genes_7C.pdf",
                   width = 10,height = 24,
                   border_color = NA,
                   fontsize_row = 4,
                   color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                   clustering_callback = callback)

