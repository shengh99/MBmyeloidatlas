#library
library(foreach)
library(SCENIC)
library(Seurat)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(harmony)

set.seed(123)
#loading Seurat object
sc <- readRDS(file = "annotation_Idents.Rds")
sc.sub <- subset(sc,downsample = 1000)


#part1 Co-expression
#SCENIC initialization
dbs = c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
        "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs) = c("500bp","10kb")
scenicOptions <- initializeScenic(org = "hgnc",
                                  dbDir = "./mc9nr/gene_based",
                                  dbs = dbs,
                                  nCores = 4)
cellInfo <- data.frame(sc.sub@meta.data)
cellInfo <- cellInfo[,c("orig.ident","seurat_clusters","module_type")]
saveRDS(cellInfo,file = "./int/cellInfo.Rds")

scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"

saveRDS(scenicOptions,file = "./int/scenicOptions.Rds")


#filter
exprMat = sc.sub@assays$RNA@counts
exprMat = as.matrix(exprMat)
genesKept <- geneFiltering(exprMat,
                           scenicOptions = scenicOptions,
                           minCountsPerGene = 3*0.01*ncol(exprMat),
                           minSamples = ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept,]
#
scenicOptions@settings$seed = 123
runCorrelation(exprMat_filtered,scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log,
          scenicOptions,
          nParts = 50,
          resumePreviousRun = T)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions,
                                        log2(as.matrix(exprMat_filtered)+1))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

saveRDS(scenicOptions,file = "int/scenicOptions_end.Rds")

Idents(sc.sub) = sc.sub@meta.data$module_type
Info <- data.frame(Module=Idents(sc.sub))
#load AUC matrix
regulonAUC <- loadInt(scenicOptions,
                      "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
#scale
regulonActivity_byModule <- sapply(split(rownames(Info),Info$Module), 
                                   function(cells)
                                     rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byModule_Scaled <- t(scale(t(regulonActivity_byModule),center = T,scale = T))
pdf("regulonActivity_byModule_Scaled.pdf",height = 8)
ComplexHeatmap::Heatmap(regulonActivity_byModule_Scaled,
                        name = "Regulon activity")

rss2 = regulonActivity_byModule_Scaled
head(rss2)
df = do.call(rbind,
             lapply(1:ncol(rss2),function(i){
               dat=data.frame(
                 path = rownames(rss2),
                 module = colnames(rss2)[i],
                 sd.1 = rss2[,i],
                 sd.2 = apply(rss2[,-i],1,median)
               )
             }))
df$fc = df$sd.1 - df$sd.2
top6 <- df %>% group_by(module) %>% top_n(6,fc)
rowcn = data.frame(path = top6$module)
n = rss2[top6$path,]

col_ann = data.frame(Module=c(colnames(n)))
rownames(col_ann) = col_ann$Module
pheatmap::pheatmap(n,
                   cluster_rows = F,
                   cluster_cols = F,
                   show_rownames = T,
                   show_colnames = T,
                   angle_col=315,
                   annotation_col = col_ann,
                   annotation_names_col = F,
                   width = 8,
                   height = 10,
                   border_color = NA,
                   filename = "regulonActivity_byModule_Scaled_top6.pdf"
)











