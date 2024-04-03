#library
Sys.umask("002")
library(survival)
library(survminer)
library(tidyverse)
library(GSVA)
library(dplyr)

expr_count <- read.csv(file = "gene_TPM.csv",sep = "\t",header = T)
meta_data <- read.csv(file = "patient_metadata.csv",sep = "\t",header = T)
gene_select <- read.csv(file = "module_gene.csv",sep = "\t",header = T)

row_blank <- which(expr_count$HGNC_symbol_from_ensemblv77 == "")
expr_count <- expr_count[-row_blank,]
expr_count <- expr_count %>% distinct(HGNC_symbol_from_ensemblv77,.keep_all=T)
row.names(expr_count) = expr_count[,1]
expr_count$HGNC_symbol_from_ensemblv77=NULL

gene_set = list()
gene_set[["module1"]] = as.character(gene_select$module1[gene_select$module1 != ""])
gene_set[["module2"]] = as.character(gene_select$module2[gene_select$module2 != ""])
gene_set[["module3"]] = as.character(gene_select$module3[gene_select$module3 != ""])
gene_set[["module4"]] = as.character(gene_select$module4[gene_select$module4 != ""])
gene_set[["module5"]] = as.character(gene_select$module5[gene_select$module5 != ""])
gene_set[["module6"]] = as.character(gene_select$module6[gene_select$module6 != ""])
gene_set[["module7"]] = as.character(gene_select$module7[gene_select$module7 != ""])
gene_set[["module8"]] = as.character(gene_select$module8[gene_select$module8 != ""])
ES <- gsva(as.matrix(expr_count),gene_set)

ES=as.data.frame(ES)
for (i in c(1:8)) {
  es = ES[i,]
  es = as.data.frame(es)
  es = t(es)
  meta_data = merge(meta_data,es,by.x="Study_ID",by.y=0)
}


meta_data$module8 = ifelse(meta_data$module8 < median(meta_data$module8),"low","high")
ggsurvplot(survfit(Surv(OS_years,Dead)~module8,
                   data = meta_data),
           conf.int = F,
           pval = T)
