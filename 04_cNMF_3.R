dir.input = "03_res1"
dir.output = "04_res2"
dir_count = "01_count_data"
usage_filter = 0.01
top_gene = 100
cor_min = 0
cor_max = 1
color = NULL
cluster_method = "complete"
library(tidyverse)

dirs = setdiff(dir(dir.input),c("k_selection.txt"))
ref.file = read.table(paste(dir.input,"/k_selection.txt",sep = ""),
                      header = F,
                      sep = ",",
                      stringsAsFactors = F)
colnames(ref.file) = c("sample","k")
rownames(ref.file) = ref.file$sample


#
for (i in dirs){
  usage.file = dir(dir.output,pattern = paste(i,".usages",sep = ""))
  usage.df = read.table(paste(dir.output,"/",usage.file,sep = ""),
                        header = T,
                        row.names = 1,
                        sep = "\t",
                        stringsAsFactors = F)
  colnames(usage.df) = paste(i,1:dim(usage.df)[2],sep = ".")
  
  #normalize
  usage.df = usage.df/rowSums(usage.df)
  write.table(usage.df,
              file = paste(dir.output,"/",i,"_program.usage.norm.txt",sep = ""),
              quote = F,
              sep = "\t",
              row.names = T,
              col.names = T)
  
  #QC
  tpmdf1 = gather(usage.df,"program","ratio")
  tpmdf1 %>% ggplot(aes(x=program,y=ratio))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(color = "red",alpha=0.4,width=0.2)+
    labs(title=i)+
    theme(axis.text.x.bottom=element_text(angle = 45,hjust=1),
          plot.title=element_text(hjust=0.5,size=20))
  ggsave(paste(dir.output,"/",i,"_program.usage.norm.QC.png",sep = ""),
         device = "png",
         width = 20,
         height = 16,
         units = c("cm"))
  
  #score
  score.file = dir(dir.output,pattern = paste(i,".gene_spectra_score",sep = ""))
  score.df = read.table(paste(dir.output,"/",score.file,sep = ""),
                        header = T,
                        row.names = 1,
                        sep = "\t",
                        stringsAsFactors = F)
  score.df=as.data.frame(t(score.df))
  colnames(score.df)=paste(i,1:dim(score.df)[2],sep=".")
  
  topn.df=as.data.frame(matrix(nrow = top_gene,ncol = ncol(score.df)))
  colnames(topn.df)=colnames(score.df)
  
  for (k in colnames(score.df)) {
    tmpv=score.df[,k]
    names(tmpv)=rownames(score.df)
    topn.df[,k]=names(rev(tail(sort(tmpv),top_gene)))
  }
  
  #save
  write.table(topn.df,
              file = paste(dir.output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = ""),
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = T)
  score.df$gene=rownames(score.df)
  write.table(score.df,
              file = paste(dir.output,"/",i,"_program.Zscore.txt",sep = ""),
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = T)
}


check.usage=data.frame()
for (i in dirs) {
  usage.file = paste(dir.output,"/",i,"_program.usage.norm.txt",sep = "")
  usage.df = read.table(usage.file,
                        header = T,
                        row.names = 1,
                        sep = "\t",
                        stringsAsFactors = F)
  check.usage = rbind(check.usage,as.data.frame(colMeans(usage.df)))
}
colnames(check.usage)=c("mean_ratio")

check.usage$sample_programs=rownames(check.usage)
check.usage=check.usage%>%arrange(mean_ratio)
check.usage$sample_programs=factor(check.usage$sample_programs,
                                   levels = check.usage$sample_programs)

linex = sum(check.usage$mean_ratio < usage_filter)
check.usage%>%ggplot(aes(x=sample_programs,y=mean_ratio))+
  geom_point()+
  geom_hline(yintercept = usage_filter,color="red")+
  geom_vline(xintercept = linex+0.5,color="red")+
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1))
ggsave(paste(dir.output,"/","check.usage.png",sep = ""),
       width = 30,
       height = 16,
       device = "png",
       units = "cm")
maybe.bg = as.character(check.usage$sample_programs[check.usage$mean_ratio<usage_filter])



library(pheatmap)
library(RColorBrewer)
library(scales)

all.score.df = data.frame()
all.score.topn.df= data.frame()

for (i in dirs) {
  score.file=paste(dir.output,"/",i,"_program.Zscore.txt",sep = "")
  score.df=read.table(score.file,
                      header = T,
                      sep = "\t",
                      stringsAsFactors = F)
  if (i==dirs[1]) {all.score.df=score.df}
  if (i!=dirs[1]) {all.score.df=all.score.df%>%inner_join(score.df,by="gene")}
  
  score.topn.file=paste(dir.output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
  score.topn.df = read.table(score.topn.file,
                             header = T,
                             sep = "\t",
                             stringsAsFactors = F)
  if (i==dirs[1]) {all.score.topn.df=score.topn.df}
  if (i!=dirs[1]) {all.score.topn.df=cbind(all.score.topn.df,score.topn.df)}
  
}

rownames(all.score.df)=all.score.df$gene
all.score.df$gene=NULL
all.score.df=all.score.df[rowSums(is.na(all.score.df))==0,]
all.score.rm.df = all.score.df[,setdiff(colnames(all.score.df),maybe.bg)]

write.table(all.score.rm.df,
            file = paste(dir.output,"/","program_genecount.txt",sep = ""),
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)

all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")


colanno=as.data.frame(colnames(all.score.rm.df.cor))
colnames(colanno) = "colnames"
colanno$Sample = str_replace(colanno$colnames,"\\..*","")
rownames(colanno)=colanno$colnames
colanno$colnames = NULL

ann_colors=list(Sample=c(CB_1="#B07AA1",CB_2="#9C755F"))

all.score.rm.df.cor_1 = all.score.rm.df.cor
all.score.rm.df.cor[all.score.rm.df.cor < -1]=-1
all.score.rm.df.cor[all.score.rm.df.cor > 1]=1

tmpp=pheatmap::pheatmap(all.score.rm.df.cor,
                        cluster_rows = T,
                        cluster_cols = T,
                        clustering_method = "complete",
                        show_colnames = F,
                        show_rownames = F,
                        treeheight_row = 0,
                        treeheight_col = 0,
                        annotation_col = colanno,
                        annotation_names_row = F,
                        annotation_names_col = F,
                        annotation_colors = ann_colors,
                        color = colorRampPalette(c("#1f78b4","white","#b21534"))(50),
                        width = 8,
                        height = 6,
                        filename = paste("CB_module_corpearson_",cluster_method,"_heatmap.pdf",sep = ""))


mat.files=dir("./04_res2/",pattern = "dt_0_1.txt$")
all.mat=data.frame()
for (fi in mat.files){
  tmp.mat=read.table(paste0("./04_res2/",fi),
                     header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  tmp.mat=as.data.frame(t(tmp.mat))
  sampleid=str_replace(fi,"\\..*$","")
  colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = ".")
  tmp.mat$gene=rownames(tmp.mat)
  if (sampleid == "CB_1"){
    all.mat=tmp.mat
  }else{all.mat=all.mat%>%full_join(tmp.mat,by="gene")}
}

#module1/signature program1[1:2]
m1=c("CB_1.1","CB_2.1")
m1.loading=all.mat[,c("gene",m1)]
used.gene=c()
for (mi in m1){
  tmp.df=m1.loading[,c("gene",mi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

m1.loading=m1.loading[m1.loading$gene %in% used.gene,]
rownames(m1.loading)=m1.loading$gene
m1.loading$gene=NULL
m1.loading[is.na(m1.loading)] <- 0
m1.loading$total_loading=rowSums(m1.loading)
m1.loading$average_loading=(m1.loading$total_loading)/length(m1)

m1.loading=m1.loading%>%arrange(desc(average_loading))

write.csv(m1.loading,file = "./04_res2/module/module1_loading.csv")
head(rownames(m1.loading),30)



#module2/signature program1[3:4]
m2=c("CB_1.5","CB_2.2")
m2.loading=all.mat[,c("gene",m2)]
used.gene=c()
for (mi in m2){
  tmp.df=m2.loading[,c("gene",mi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

m2.loading=m2.loading[m2.loading$gene %in% used.gene,]
rownames(m2.loading)=m2.loading$gene
m2.loading$gene=NULL
m2.loading[is.na(m2.loading)] <- 0
m2.loading$total_loading=rowSums(m2.loading)
m2.loading$average_loading=(m2.loading$total_loading)/length(m2)

m2.loading=m2.loading%>%arrange(desc(average_loading))

write.csv(m2.loading,file = "./04_res2/module/module2_loading.csv")
head(rownames(m2.loading),30)

#module3/signature program1[6:7]
m3=c("CB_1.3","CB_2.4")
m3.loading=all.mat[,c("gene",m3)]
used.gene=c()
for (mi in m3){
  tmp.df=m3.loading[,c("gene",mi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

m3.loading=m3.loading[m3.loading$gene %in% used.gene,]
rownames(m3.loading)=m3.loading$gene
m3.loading$gene=NULL
m3.loading[is.na(m3.loading)] <- 0
m3.loading$total_loading=rowSums(m3.loading)
m3.loading$average_loading=(m3.loading$total_loading)/length(m3)

m3.loading=m3.loading%>%arrange(desc(average_loading))

write.csv(m3.loading,file = "./04_res2/module/module3_loading.csv")
head(rownames(m3.loading),30)


#module4/signature program1[8:9]
m4=c("CB_1.6","CB_2.3")
m4.loading=all.mat[,c("gene",m4)]
used.gene=c()
for (mi in m4){
  tmp.df=m4.loading[,c("gene",mi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

m4.loading=m4.loading[m4.loading$gene %in% used.gene,]
rownames(m4.loading)=m4.loading$gene
m4.loading$gene=NULL
m4.loading[is.na(m4.loading)] <- 0
m4.loading$total_loading=rowSums(m4.loading)
m4.loading$average_loading=(m4.loading$total_loading)/length(m4)

m4.loading=m4.loading%>%arrange(desc(average_loading))

write.csv(m4.loading,file = "./04_res2/module/module4_loading.csv")
head(rownames(m4.loading),30)

#module5/signature program1[10:13]
m5=c("CB_1.4","CB_1.8","CB_1.7","CB_2.6")
m5.loading=all.mat[,c("gene",m5)]
used.gene=c()
for (mi in m5){
  tmp.df=m5.loading[,c("gene",mi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

m5.loading=m5.loading[m5.loading$gene %in% used.gene,]
rownames(m5.loading)=m5.loading$gene
m5.loading$gene=NULL
m5.loading[is.na(m5.loading)] <- 0
m5.loading$total_loading=rowSums(m5.loading)
m5.loading$average_loading=(m5.loading$total_loading)/length(m5)

m5.loading=m5.loading%>%arrange(desc(average_loading))

write.csv(m5.loading,file = "./04_res2/module/module5_loading.csv")
head(rownames(m5.loading),30)

m.top100 = data.frame(col1=c(1:100))
m.top100$module1=head(rownames(m1.loading),100)
m.top100$module2=head(rownames(m2.loading),100)
m.top100$module3=head(rownames(m3.loading),100)
m.top100$module4=head(rownames(m4.loading),100)
m.top100$module5=head(rownames(m5.loading),100)
m.top100=m.top100[,-1]

write.csv(m.top100,file = "./04_res2/module/top100/top100gene.csv")
write.table(m.top100,
            file = "./04_res2/module/top100/top100gene.txt",
            quote = F,
            sep = "\t",
            col.names = T,
            row.names = T)



################################################################################
#remove all
#UMAP cluster mapping for MB modules
#load Rds
sc <- readRDS(file = "TAM_PCA.Rds")

signature <- read.xlsx("tumor_top30gene_8module.xlsx",
                       sheetIndex = 1,
                       header = T)
signature$index = rownames(signature)
signature = signature%>%melt(id="index")
signature$index = NULL
signature = unique(signature)
colnames(signature) = c("set","gene")
signature=signature[signature$gene %in% rownames(sc),]

#calculate gene set scores
for (i in as.character(unique(signature$set))) {
  signature_small = signature%>%filter(set==i)
  genes.for.scoring <- list(signature_small$gene)
  sc <- AddModuleScore(object = sc,
                       features = genes.for.scoring,
                       name = i)
  
}
colnames(sc@meta.data)[12:19]=c("module1","module2","module3","module4","module5","module6","module7","module8")
sc@meta.data$barcode <- rownames(sc@meta.data)
anno.new <- sc@meta.data[,c("orig.ident","Subgroup","module1","module2","module3","module4","module5","module6","module7","module8","barcode")]   

#re-create Seurat object
TAM.module <- CreateSeuratObject(counts = sc[["RNA"]]@counts)
TAM.module@meta.data$barcode = rownames(TAM.module@meta.data)
TAM.module@meta.data = inner_join(TAM.module@meta.data,anno.new,by="barcode")
rownames(TAM.module@meta.data) = TAM.module@meta.data$barcode

TAM.module@meta.data$orig.ident.x=NULL
colnames(TAM.module@meta.data)[4] = "orig.ident"
TAM.module[["percent_MT"]] <- PercentageFeatureSet(TAM.module,pattern = "^MT-")

TAM.module <-
  NormalizeData(
    TAM.module,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
TAM.module <-
  FindVariableFeatures(
    TAM.module,
    selection.method = "vst",
    nfeatures = 3000
  )
TAM.module <-
  ScaleData(
    TAM.module,
    features = rownames(TAM.module),
    vars.to.regress = c("percent_MT","nCount_RNA"),
    do.center = T
  )
TAM.module <-
  RunPCA(
    TAM.module,
    features = VariableFeatures(TAM.module)
  )
ElbowPlot(TAM.module,ndims = 60)
TAM.module.harmony <- RunHarmony(
  TAM.module,
  group.by.vars = "orig.ident",
  reduction = "pca",
  max.iter.harmony = 100,
  lambda = 0.5,
  reduction.save = "harmony",
  assay.use = "RNA"
)
TAM.module.harmony <- FindNeighbors(TAM.module.harmony,
                                    reduction = "harmony",
                                    dims = 1:25) %>% FindClusters(resolution = 0.45)
TAM.module.harmony <- RunUMAP(TAM.module.harmony,
                              reduction = "harmony",
                              dims = 1:30)
#score
df.new = TAM.module.harmony@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  summarize(mean(module1),
            mean(module2),
            mean(module3),
            mean(module4),
            mean(module5),
            mean(module6),
            mean(module7),
            mean(module8))
colnames(df.new)[2:9] = c("module1","module2","module3","module4","module5","module6","module7","module8")

df.new = as.data.frame(df.new)
rownames(df.new) = as.character(df.new$seurat_clusters)
df.new$seurat_clusters = NULL

df.new = as.data.frame(scale(df.new))
df.new[df.new > 2] = 2
df.new[df.new < (-2)] = -2

df.new$cluster = rownames(df.new)
df.new = df.new %>% reshape2::melt(id.vars = c("cluster"))
colnames(df.new)[2:3] = c("module","score_scaled")

df.new %>% ggplot(aes(x=cluster,y=module,fill=score_scaled))+
  geom_tile()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete("",expand = c(0,0))+
  scale_fill_gradient2("module\nscore",low = "#1f78b4",mid = "#ffffff",high = "#b21534",breaks=-2:2,labels=-2:2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.x.bottom = element_text())

