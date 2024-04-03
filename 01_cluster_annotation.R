#library
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)
library(harmony)
library(scRNAtoolVis)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)

#load data
sc_1.data <- Read10X(data.dir = "./outs1/filtered_feature_bc_matrix/")
sc_2.data <- Read10X(data.dir = "./outs2/filtered_feature_bc_matrix/")


#CreateSeuratObject
sc_1 <- CreateSeuratObject(counts = sc_1.data,min.cells = 3,min.features = 200)
sc_2 <- CreateSeuratObject(counts = sc_2.data,min.cells = 3,min.features = 200)

#merge into a single dataset
alldata <- merge(sc_1,y=c(sc_2),add.cell.ids = c("CB_1","CB_2"))

#mitochondrial genes
alldata[["percent_MT"]] <- PercentageFeatureSet(alldata,pattern = "^MT-")
alldata[["percent_RP"]] <- PercentageFeatureSet(alldata,pattern = "^RP[SL]")
hb <-c("HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBE1",
    "HBG1",
    "HBG2",
    "HBM",
    "HBQ1",
    "HBZ")
hb <- CaseMatch(hb, rownames(alldata))
alldata[["percent_HB"]] <- PercentageFeatureSet(alldata,features = hb)
VlnPlot(alldata,features = c("nFeature_RNA","nCount_RNA","percent_MT","percent_RP","percent_HB"),ncol = 3,pt.size = 0)

#QC
alldatafiltered <- subset(alldata,subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent_MT < 25 & nCount_RNA > 1000 & nCount_RNA < 30000 & percent_HB < 1)
saveRDS(alldatafiltered,file = "./01_QC/01_QC.Rds")


#LogNormalize
alldatafiltered <- NormalizeData(alldatafiltered,normalization.method = "LogNormalize",scale.factor = 10000)
alldatafiltered <- FindVariableFeatures(alldatafiltered,selection.method = "vst",nfeatures = 3000)
alldatafiltered <- 
  ScaleData(
    alldatafiltered,
    vars.to.regress = c("percent_MT","nCount_RNA"),
    do.center = T)
alldatafiltered <- RunPCA(
  alldatafiltered,
  features = VariableFeatures(alldatafiltered))
ElbowPlot(alldatafiltered,ndims = 60)

#harmony
alldatafiltered.harmony <-
  RunHarmony(
    alldatafiltered,
    group.by.vars = "orig.ident",
    reduction = "pca",
    max.iter.harmony = 100,
    lambda = 0.5,
    reduction.save = "harmony",
    assay.use = "RNA"
  )
alldatafiltered.harmony <- 
  RunUMAP(
    alldatafiltered.harmony,
    dims = 1:30,
    reduction = "harmony")
alldatafiltered.harmony <-
  FindNeighbors(
    alldatafiltered.harmony,
    dims = 1:30,
    reduction = "harmony"
  )
alldatafiltered.harmony <- 
  FindClusters(
    alldatafiltered.harmony,
    resolution = 0.2
  )
DimPlot(
  alldatafiltered.harmony,
  reduction = "umap",
  label = T)
FeaturePlot(alldatafiltered.harmony,
            features = c("AIF1","PTPRC","C1QA","SPP1","ZIC1","PAX6","PROM1","NEUROD1","TOP2A"),
            ncol = 3)
FeaturePlot(alldatafiltered.harmony,
            features = c("CD1C","FCER1A","P2RY12","FCN1","S100A8","S100A9","CD163","GNLY","OLIG2","MBP"),
            ncol = 3)
#filter non-immune cells
immunecell = alldatafiltered.harmony[
  ,alldatafiltered.harmony@meta.data$seurat_clusters %in% c(0,1,4,6,7,8)]
##recluster
immunecell <-
  NormalizeData(
    immunecell,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
immunecell <-
  FindVariableFeatures(
    immunecell,
    selection.method = "vst",
    nfeatures = 3000
  )
immunecell <-
  ScaleData(
    immunecell,
    features = rownames(immunecell),
    vars.to.regress = c("percent_MT","nCount_RNA"),
    do.center = T
  )
immunecell <-
  RunPCA(
    immunecell,
    features = VariableFeatures(immunecell)
  )
ElbowPlot(immunecell,ndims = 60)
immunecell.harmony <-
  RunHarmony(
    immunecell,
    group.by.vars = "orig.ident",
    reduction = "pca",
    max.iter.harmony = 100,
    lambda = 0.5,
    reduction.save = "harmony",
    assay.use = "RNA"
  )
immunecell.harmony <- 
  RunUMAP(
    immunecell.harmony,
    dims = 1:20,
    reduction = "harmony")
immunecell.harmony <-
  FindNeighbors(
    immunecell.harmony,
    dims = 1:20,
    reduction = "harmony"
  )
immunecell.harmony <- 
  FindClusters(
    immunecell.harmony,
    resolution = 0.2
  )
DimPlot(
  immunecell.harmony,
  reduction = "umap",
  label = T)
sc.markers <-
  FindAllMarkers(
    immunecell.harmony,
    only.pos = T,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = F
  )
saveRDS(immunecell.harmony,file = "immunecell.Rds")

#remove all
#annotation
sc <- readRDS(file = "immunecell.Rds")
clustercelltype <- c("0"="MG1",
                     "1"="MG2",
                     "2"="MG3",
                     "3"="BAM",
                     "4"="Nt",
                     "5"="Mo",
                     "6"="NK/T")
sc[["cell_type"]] = unname(clustercelltype[sc@meta.data$seurat_clusters])
sc.annotation <- RenameIdents(sc,
                              "0"="MG1",
                              "1"="MG2",
                              "2"="MG3",
                              "3"="BAM",
                              "4"="Nt",
                              "5"="Mo",
                              "6"="NK/T")
saveRDS(sc.annotation,file = "immunecell_annotation.Rds")

#color
color_cluster = c("#66C5CC",'#F6CF71',"#F89C74","#87C55F","#9EB9F3","#FE88B1","#C9DB74")
names(color_cluster) = c("MG1","MG2","MG3","BAM","Nt","Mo","NK/T")
#UMAP plot
DimPlot(sc,
        reduction = "umap",
        group.by = "cell_type",
        label = T,
        label.size = 5,
        pt.size = 0.1,
        #label.color = color_cluster[unique(sc@meta.data$seurat_clusters)]
) + 
  scale_color_manual(values = color_cluster) + 
  theme(plot.title = element_blank())

Idents(sc) = "cell_type"
umap_data <- DimPlot(sc,reduction = "umap")
umap_data <- umap_data$data
head(umap_data)

p1 = ggplot() + 
  geom_point(data = umap_data,
             mapping = aes(UMAP_1,UMAP_2,color=ident),
             size = 0.5) +
  scale_color_manual(values = color_cluster) + 
  theme_classic()

segment.df = data.frame(x=c(-7,-7),
                        xend=c(-3.5,-7),
                        y=c(-13,-13),
                        yend=c(-13,-9))
p2=p1+geom_segment(data = segment.df,
                   mapping = aes(x=x,xend=xend,y=y,yend=yend),
                   arrow=arrow(length = unit(0.2,"cm")),
                   linewidth = 1) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x.bottom = element_text(hjust = 0.1,size = 18),
        axis.title.y.left = element_text(hjust = 0.1,size = 18))

genes <- read.table(file = "markergenes.csv",header = T)
markers <- as.character(genes$gene)
annogene = c("LRMDA","KCNQ3","C1QB","C1QC","JUN","CCL3","CD163","MRC1","S100A8","S100A9","VCAN","FCN1","NKG7","IL32")
AverageHeatmap(sc,
               markerGene = markers,
               group.by = "cell_type",
               annoCol = T,
               myanCol = color_cluster,
               column_names_rot = 0,
               row_title = "",
               border = T,
               showRowNames = F,
               markGenes = annogene)

#cluster proportion
cellratio <- prop.table(table(Idents(sc),sc$orig.ident),margin = 2)
cellratio <- as.data.frame(cellratio)
colnames(cellratio) = c("Cluster","Sample","Ratio")
cellratio %>% openxlsx::write.xlsx("sample_cluster_ratio.xlsx")
ggplot(cellratio) + 
  geom_bar(aes(x = Sample,y = Ratio,fill = Cluster),stat = "identity",width = 0.6,size = 0.5,colour = '#222222') + 
  theme_classic() + 
  RotatedAxis()+
  labs(x='',y='') + 
  scale_fill_manual(values = color_cluster) +
  scale_y_continuous(labels = percent)+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5,linetype = "solid"))



#GO enrichment analysis
for (i in unique(sc.markers$cluster)) {
  deg = sc.markers[sc.markers$cluster == i,]
  gene_up=deg[deg$p_val_adj < 0.05,]
  gene_up=gene_up[order(gene_up$avg_log2FC,decreasing = T),]
  top = gene_up[1:100,]
  df_name <- bitr(top$gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  #GO_BP
  go <- enrichGO(gene = unique(df_name$ENTREZID),
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable = T)
  go_simp = clusterProfiler::simplify(go)
  go_simp_df = data.frame(go_simp)
  write.csv(go_simp_df,file = paste0("./cluster_GO/",i,"_GOBP_simplify.csv"))
}
dt <- read.csv("Mg2Mg1_select_GOBP.csv",header = T,sep = "\t")
dt$X=NULL
dt$group <- ""
dt$group[which(dt$logpvalue > 0)] = "MG2 UP"
dt$group[which(dt$logpvalue < 0)] = "MG1 UP"
dt$Description <- factor(dt$Description,levels = rev(dt$Description))
p <- ggplot(dt,
            aes(x=logpvalue,y=Description,fill=group))+
  geom_col()+theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "grey",size = 1),
        axis.text = element_text(size = 12))
up <- dt[which(dt$logpvalue > 0),]
down <- dt[which(dt$logpvalue < 0),]
p1 <- p+geom_text(data = up,
                  aes(x=-0.2,y=Description,label=Description),
                  size=3.5,
                  hjust = 1)+
  geom_text(data = down,
            aes(x=0.2,y=Description,label=Description),
            size=3.5,
            hjust = 0)+
  scale_x_continuous(breaks = seq(-5,5,2.5))+
  theme(legend.position = "top")

#DEG volcano
deg <- FindMarkers(sc,
                   ident.1 = "MG2",
                   ident.2 = "MG1",
                   logfc.threshold = 0,
                   min.pct = 0.1)
write.csv(deg,file = "Mg2vsMg1_valcano_DEG.csv",quote = F,row.names = T)

log2FC = 0.25
padj = 0.05
k1 = deg$p_val_adj < padj & deg$avg_log2FC < (-log2FC)
k2 = deg$p_val_adj < padj & deg$avg_log2FC > log2FC
table(k1)
table(k2)
deg$change= ifelse(k1,"down",ifelse(k2,"up","ns"))

deg <- deg%>%mutate(Difference = pct.1 - pct.2)
deg <- data.frame(
  gene = rownames(deg),
  deg
)
ggplot(deg,aes(x=Difference,y=avg_log2FC,color = change))+
  geom_point(size=0.4)+
  scale_color_manual(values = c("blue","grey","red"))+
  geom_text_repel(data = subset(deg,avg_log2FC >= 1.5 & Difference >= 0.5 & p_val_adj <= 0.05),
                  mapping = aes(label=gene),
                  color="black",
                  max.overlaps = 100,
                  segment.size = 0.2,
                  size=4,
                  force = 3,
                  force_pull = 1)+
  geom_text_repel(data = subset(deg,avg_log2FC <= -0.5 & Difference <= -0.05 & p_val_adj <= 0.05),
                  mapping = aes(label=gene),
                  color="black",
                  max.overlaps = 100,
                  segment.size=0.2,
                  size=4,
                  force = 10,
                  force_pull = 1)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  ylab(bquote(avg_log[2]*FC))+
  theme_classic()









