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
setwd("e://ac/bio/training20200515")
## 打开三表，开始进行差异化表达分析
load(file = "e://ac/bio/training20200515/data/geo-cesc/prepare.rdata")
## 创建设计矩阵  
design <- model.matrix(~ 0 + sample_info$group)
colnames(design) <- levels(factor(sample_info$group))
row.names(design) <- row.names(sample_info)

## 构建对比分组
contrasts<- makeContrasts(
  CN=Cancer-Normal,
  C1N=CIN1-Normal,
  C2N=CIN2-Normal,
  C3N=CIN3-Normal,
  levels = design
)
## 差异表达分析
fit <- lmFit(gene_exp, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
de_result <- topTable(fit,   ## topTable只能每次处理一组的对比
                      coef = 'CN',  ## 所以用了“coef = "CN”
                      number = Inf) ## “inf”是无限大的意思
de_result[1:5,1:5]
## 差异基因的注释
de_result<- rownames_to_column(de_result,var = "probe") %>% 
  dplyr::select(-t,-B) %>% 
  dplyr::mutate(direction = if_else(
    adj.P.Val > 0.01,"ns",if_else(
    abs(logFC) < 2,"ns",if_else(
      logFC >= 2,"up","down"## 注意括号ifelse(up,down)
  )))) %>% 
left_join(gene_info,by = c("probe"="ID")) %>% 
  left_join(rownames_to_column(gene_exp,var = "probe"),by="probe")
de_result<- arrange(de_result)


## 画火山图
top30 <- arrange(de_result,desc(abs(logFC))) %>% 
  dplyr::slice(1:30) %>% 
EnhancedVolcano::EnhancedVolcano(de_result,
                                 x = "logFC",
                                y = "P.Value",
                                lab = de_result$Gene_Symbol,
                                selectLab = top30$Gene_Symbol,
                                pCutoff = 0.01,
                                FCcutoff = 1,
                                title = "My Volance plot",
                                subtitle = "Cancer vs Normal")+
  theme_classic()

## 筛选差异最大的基因(正常组与癌症组的对比)
Name <- rownames_to_column(sample_info,var = "Gene_symbol") %>% 
  dplyr::select(Gene_symbol,group) %>% 
  dplyr::filter(group == "Normal"| group == "Cancer")

top_exp <- arrange(de_result,desc(abs(logFC))) %>% 
  slice(1:30) %>% 
  dplyr::select(Gene_Symbol,one_of(Name$Gene_symbol)) %>% #现在还是tibble
column_to_rownames(var = "Gene_Symbol")

data.frame()

## 思路：先写需要的名单，符合Normal&Cancer
##       再写需要的data.frame，记得用one of 
  
                   

## 画热图
pheatmap::pheatmap(top_exp,
                   color = colorRampPalette(c("yellow","orange","red"))(100),
                   show_colnames = F,
                   cutree_cols = 2,
                   annotation_col = dplyr::select(sample_info,group),
                   annotation_colors = list(
                     group=c(Normal = "yellow",Cancer = "red")
                   ))



install.packages("devtools")
library(usethis)
use_git_config(user.name="petertao0703", user.email="1074270679@qq.com")
