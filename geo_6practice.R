#安装kegg
setwd("E:/ac/bio/training20200515")
load(file = 'e://ac/bio/training20200515/data/geo-cesc/prepare.rdata')
load(file = 'data/geo-cesc/de.rdata')


#dir.create("R_Library",recursive = T)
#install.packages("KEGG.db_1.0.tar.gz",repos = NULL,lib = "R_Library")
library(KEGG.db,lib = "R_Library")
library(tidyverse)
deg <- dplyr::select(de_result,ENTREZ_GENE_ID,logFC,direction) %>%
  dplyr::filter(!is.na(ENTREZ_GENE_ID)) %>% 
  mutate(ENTREZ_GENE_ID = str_split(ENTREZ_GENE_ID," ///",simplify = T)[,1]) %>% 
  distinct(ENTREZ_GENE_ID,.keep_all = T)
## 思路：提取表中的ENtrer的ID,还有logFC和发展方向。排除①ID是na；②一个探针多个ID（
# ③去重id

## 为KEGG做准备，需要(direction没有ns的)基因名称还有 有pathway的基因
gene <- filter(deg,direction !="ns") %>% 
  pull(ENTREZ_GENE_ID)

genelist <- deg$logFC
  names(genelist) <- deg$ENTREZ_GENE_ID
genelist <- sort(genelist,decreasing = T)

## 开始KEGG
library(org.Hs.eg.db)
library(clusterProfiler)
de_ekp<-enrichKEGG(gene = gene,
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   use_internal_data = T)
de_ekp <- setReadable(de_ekp,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
de_ekp_df <- as.data.frame(de_ekp)
