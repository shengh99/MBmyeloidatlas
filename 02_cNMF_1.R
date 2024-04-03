#load data
TAM <- sc[,sc@meta.data$cell_type %in% c("MG1","MG2","MG3","BAM")]

#filter ribosomal and mitochondrial genes
TAM_sub1 = list()
TAM_sub2 = list()
TAM_sub3 = list()
TAM_count = list()
gene <- readRDS(file = "gencode.v38.annotation.gtf.filtered.coding_gene.Rds")
for (i in unique(TAM$orig.ident)) {
  TAM_sub1[[i]] = subset(TAM,orig.ident == i)
  TAM_sub2[[i]] = TAM_sub1[[i]][intersect(row.names(TAM_sub1[[i]]),gene$gene_name),]
  TAM_sub3[[i]] = TAM_sub2[[i]][-c(grep(pattern = "^RP[SL]",x=rownames(TAM_sub2[[i]])),
                                   grep(pattern = "^MT-",x=rownames(TAM_sub2[[i]]))),]
  TAM_count[[i]] = t(TAM_sub3[[i]]@assays$RNA@counts %>% as.matrix())
  
  write.table(TAM_count[[i]],
              file = paste0(i,".count.txt"),
              sep='\t')
  saveRDS(TAM_sub3[[i]],
          file = paste0(i,".Rds"))
  
}
