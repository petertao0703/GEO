options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
library(tidyverse)
library(cowplot)
library(ggthemes)
library(GEOquery)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(PCAtools)
library(WGCNA)
library(genefilter)
load(file = 'E:/ac/bio/training20200515/data/geo-cesc/input.rdata')
gset = gset[[1]]
## 样本信息表
sample_info = pData(gset) %>% 
select(geo_accession,title) %>%
  mutate(group = str_split(title,"-",simplify = T)[,1]) %>% 
  mutate(group_num = case_when(
    group == 'Normal'~ 1,
    group == 'CIN1'~ 2,
    group == 'CIN2'~ 3,
    group == "CIN3"~ 4,
    group == "Cancer"~ 5)) %>% 
mutate(test1 = round(runif(128,1,100),digits = 2),
       test2 = round(runif(128,1,100),digits = 2),
       test3 = round(runif(128,1,100),digits = 2),
       test4 = round(runif(128,1,100),digits = 2))  
group_by(sample_info,group) %>% 
summarise(count=n())
## 基因矩阵
genexp <- as.data.frame(exprs(gset)) 
genexp <- rownames_to_column(genexp,var = "probe")
genexp <- dplyr::filter(genexp,!str_detect(probe,"AFFX"))
## 二表关联
genexp <- genexp[,which(colnames(genexp)
                        %in% row.names(sample_info))]
## 画图看相关性
par(mar=c(10,5,3,3))
boxplot(genexp,
        outline = F,
        las =2)
## 箱线图标准化
library(limma)
temp <- normalizeBetweenArrays(genexp)
boxplot(temp,
        outline = F,
        las=2)


