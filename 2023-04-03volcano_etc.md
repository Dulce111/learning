# 4.3 volcano & KEGG富集分析学习
1. 继续volcano: [CSDN](https://blog.csdn.net/qq_35294674/article/details/121399594)

`head(data_volcano)`   
`pvalue = 0.05`  
`data_volcano$sig[-1*log10(data_volcano$SF0M.HC.P.value) < -1*log10(pvalue) ] <- "Not sig"`  

`data_volcano$sig[-1*log10(data_volcano$SF0M.HC.P.value) > -1*log10(pvalue) & data_volcano$SF0M.HC.Ratio>0 ] <- "up"`

`data_volcano$sig[-1*log10(data_volcano$SF0M.HC.P.value) > -1*log10(pvalue) & data_volcano$SF0M.HC.Ratio<0 ] <- "down"`
`head(data_volcano)`
 - 暂时不考虑标记 

 #三段式作图结构  ：
  - 设置图片格式 不要忘记在文件名后加后缀并且与所用函数一致  
`pdf("0403volcano1")`  
 - 作图语句处于中间  
`ggplot(data_volcano, aes(SF0M.HC.Ratio, -1*log10(SF0M.HC.P.value))) + 
  geom_point(aes(color = sig)) + 
  labs(title = "SF0M V.S. HC",
       x = "log2(Fold Change)", 
       y = "-log10(pvalue)" ) + 
  scale_color_manual( values = c("#546de5", "#d2dae2","#ff4757")) + 
  geom_hline(yintercept = -log10(pvalue), linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2)`
  - 最后是结束画图语句  ps: dev.cur---useful  
`dev.off()`
 - 查看导出图片位置（默认工作目录当中）  
`getwd()`  
##pvalue线下出现红点？---（在Linux中运行可以）

## 0403 学习用R进行GO/KEGG分析
  学习来源：[知乎](https://zhuanlan.zhihu.com/p/561522453#:~:text=%E8%BD%AC%E5%BD%95%E7%BB%84%E9%89%B4%E5%AE%9A%E4%BA%86%E5%B7%AE%E5%BC%82,%E4%B8%8E%E8%A1%A8%E5%9E%8B%E8%BF%9B%E8%A1%8C%E8%81%94%E7%B3%BB%E3%80%82)  

 ### R包cluster Profiler
 - 安装 by *bioconductor*
`BiocManager::install('clusterProfiler')`
 
  - 使用R中带有的hg19等参考基因组进行分析
  暂略
`library("clusterProfiler")`
 - 无参分析 需要输入数据、背景基因集（注释好的GO）和基因列表：  
  - 读取数据  
`getwd()`  
`setwd('D:/test_GO_KEGG/GO_nopara')`  
`gene_GO <- read.delim("gene_GO.txt", stringsAsFactors = F)`  
`genes <- read.delim('gene.txt')$gene_id`   
 - 函数 ***enricher*** 的使用
`go_rich <- enricher(gene = genes, 
                    TERM2GENE = gene_GO[c('ID', 'gene_id')],
                    TERM2NAME = gene_GO[c('ID', 'Description')],
                    pAdjustMethod = 'BH', 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2)
write.table(go_rich, 'go_rich.txt', sep = '\t', row.names = F,quote = F)`  
`help(enricher)`
> A universal enrichment analyser   
>  **Useage**:   
enricher(
  gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  gson = NULL,
  TERM2GENE,
  TERM2NAME = NA
)
 
 - 读入数据
`tmp <- read.delim('go_rich.txt')`  
 - 分析
`gene_GO <- gene_GO[!duplicated(gene_GO$ID), ]`  
`tmp <- merge(tmp, gene_GO[c('ID', 'ONTOLOGY')], by = 'ID')`  
`tmp <- tmp[c(10, 1:9)]`  
`tmp <- tmp[order(tmp$pvalue), ]`  
`write.table(tmp, 'go_rich.add_Ontology.txt', sep = '\t', row.names = F, quote = F)`  
 - **可视化**：用cluster profiler的默认作图(fail)  
barplot(tmp)
help('barplot')
 - **用ggplot2**绘图  √

`tmp$term <- paste(tmp$ID, tmp$Description, sep = ': ')`  
`tmp <- tmp[order(tmp$ONTOLOGY, tmp$p.adjust), ]`   
`tmp$term <- factor(tmp$term, levels = tmp$term)`

`ggplot(tmp, aes(term, -log10(p.adjust))) + 
  geom_col(aes(fill = ONTOLOGY), width = 0.5, show.legend = FALSE) + 
  scale_fill_manual(values = c('#D06660', '#5AAD36', '#6C85F5')) + 
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') + 
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = 'black', fill = 'transparent')) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() + 
  labs(x = '', y = '-log10 P-Value\n')`   

### 无参KEGG分析:
 > 同样需要先对基因进行功能注释，背景基因集需要包含基因名称、通路和功能描述。再额外提供待富集的基因列表（差异基因）
 - 读取数据
`getwd()`  
`setwd('D:/test_GO_KEGG/KEGG_nopara')`  
`gene_KEGG <- read.delim('gene_KEGG.txt', stringsAsFactors = F)`  
`genes <- read.delim('gene.txt', stringsAsFactors = F)$gene_id`  

 - KEGG富集分析  
`kegg_rich <- enricher(gene = genes, 
                      TERM2GENE = gene_KEGG[c('Pathway', 'gene_id')],
                      TERM2NAME = gene_KEGG[c('Pathway', 'Description')],
                      pAdjustMethod = 'BH', 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)`  

`write.table(kegg_rich, 'kegg_rich.txt', sep = '\t', row.names = F, quote = F)`  
`getwd()`

#作图与之前相似