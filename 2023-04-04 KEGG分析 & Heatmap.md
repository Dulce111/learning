# 0404-KEGG&GO富集分析尝试
- 无参进行的示例数据分析  
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
`help('enricher')`

#作图与之前相似

 - 没懂KEGG Mapper分析？直接给出通路网络图？


- [GO analysis](https://github.com/eggnogdb/eggnog-mapper)

## 继续KEGG，自己数据尝试 & learning
`help("read.delim")`  
### 1. 读取处理过的数据（包含gene name & accession number &pvalue&ratio等）

 - 跳过第一行编号?  
    A: 从输出修改(`row.names = F`)，目前数据不包含行序号，可用  

 - 读取数据  
`getwd()`   
`setwd('D:/DATA_proteomics/filter_pvalue_data_differencially')`   
`filteredSF0C.HC <- read.table('SF0C.HC.csv', header = T, sep = ',')`   

### 2. 将accession number转换为EntrezID
***bitr***函数   

`help('bitr')`  

 - 查看org.Hs.eg.db中的gene name类型

`BiocManager::install("org.Hs.eg.db")`   
`library(clusterProfiler)`   
`library(AnnotationDbi)`    
`library(stats4)`   
`library(BiocGenerics)`   
`library(Biobase)`   
`library(IRanges)`   
`library(S4Vectors)`    
`library(org.Hs.eg.db)`   
`keytypes(org.Hs.eg.db)`


 - **Q**: 失灵了？？？不能在数据框后面加入entrezID???,第一次成功  
**A**: 修改参数：`drop = F` 即可  
`filteredSF0C.HC_EntrezID <- bitr(filteredSF0C.HC$Protein.accession, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")`    
`filteredSF0C.HC$EntrezID <- bitr(filteredSF0C.HC$Protein.accession, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = F)`  

 - 尝试纠正：  
`help(rm)`  
`rm(filteredSF0C.HC_EntrezID)`  
`rm(filteredSF0C.HC)`  
`write.csv(filteredSF0C.HC, "test_entrezid.csv", row.names = F)`  

 - 转换成功：`enrichGO` 进行下一步分析  

`help("enrichGO")`  

`GO_result_SF0C.HC <- enrichGO(filteredSF0C.HC$EntrezID,OrgDb = "org.Hs.eg.db", 
                              keyType = "ENTREZID", ont = "ALL", readable = T)`   

###作图？？是否导出GO结果
`barplot(GO_result_SF0C.HC,showCategory = 20)`  
`dotplot(GO_result_SF0C.HC, showCategory = 50)` 
### 输出结果为NULL 

- **富集到GO通路不显著**  
##尝试富集KEGG通路：  
`help(enrichKEGG) `   
`gene = filteredSF0C.HC$EntrezID`    
`KEGG_result_SF0C.HC <- enrichKEGG(filteredSF0C.HC$EntrezID,keyType = "kegg", organism = "hsa", pvalueCutoff = 0.1)`  

##尝试昨天学习的代码，无参分析，×,需要提前注释通路、描述等信息
`kegg_result_SF0C.HC <- enricher(filteredSF0C.HC$EntrezID, 
                      TERM2GENE = gene_KEGG[c('Pathway', 'gene_id')],
                      TERM2NAME = gene_KEGG[c('Pathway', 'Description')],
                      pAdjustMethod = 'BH', 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)`  

###GO&KEGG的区别 √
# 4.7整理配置Mfuzz数据
`setwd("D:/DATA_proteomics/100_mfuzz_data")`  

`mfuzz_t0_100 <- read.table("HC_SF0M_SF0S_SF0C_Mfuzz.txt", header = T, sep = "\t")`   

`mfuzz_t0_100_cluster1 <- subset(mfuzz_t0_100, mfuzz_t0_100$Regulated._Type == "Cluster1")
write.csv(mfuzz_t0_100_cluster1, file = "mfuzz_100_t0_1.csv", row.names = F)`    
### Repeat
`mfuzz_t0_100 <- read.table("HC_SF0M_SF0S_SF0C_Mfuzz.txt", header = T, sep = "\t")`   

`mfuzz_t0_100_cluster2 <- subset(mfuzz_t0_100, mfuzz_t0_100$Regulated._Type == "Cluster2")
write.csv(mfuzz_t0_100_cluster2, file = "mfuzz_100_t0_2.csv", row.names = F)`   


`mfuzz_t0_100 <- read.table("HC_SF0M_SF0S_SF0C_Mfuzz.txt", header = T, sep = "\t")`   
`mfuzz_t0_100_cluster3 <- subset(mfuzz_t0_100, mfuzz_t0_100$Regulated._Type == "Cluster3")
write.csv(mfuzz_t0_100_cluster3, file = "mfuzz_100_t0_3.csv", row.names = F)`  

`mfuzz_t0_100 <- read.table("HC_SF0M_SF0S_SF0C_Mfuzz.txt", header = T, sep = "\t")`   

`mfuzz_t0_100_cluster4 <- subset(mfuzz_t0_100, mfuzz_t0_100$Regulated._Type == "Cluster4")
write.csv(mfuzz_t0_100_cluster4, file = "mfuzz_100_t0_4.csv", row.names = F)`

##
`mfuzz_t1_100 <- read.table("HC_SF1M_SF1S_SF1C_Mfuzz.txt", header = T, sep = "\t")`   
`mfuzz_t1_100_cluster1 <- subset(mfuzz_t1_100, mfuzz_t1_100$Regulated.Type == "Cluster1")
write.csv(mfuzz_t1_100_cluster1, file = "mfuzz_100_t1_1.csv", row.names = F)`  


`mfuzz_t1_100_cluster2 <- subset(mfuzz_t1_100, mfuzz_t1_100$Regulated.Type == "Cluster2")
write.csv(mfuzz_t1_100_cluster2, file = "mfuzz_100_t1_2.csv", row.names = F)`   

`mfuzz_t1_100_cluster3 <- subset(mfuzz_t1_100, mfuzz_t1_100$Regulated.Type == "Cluster3")
write.csv(mfuzz_t1_100_cluster3, file = "mfuzz_100_t1_3.csv", row.names = F)`   


`mfuzz_t1_100_cluster4 <- subset(mfuzz_t1_100, mfuzz_t1_100$Regulated.Type == "Cluster4")
write.csv(mfuzz_t1_100_cluster4, file = "mfuzz_100_t1_4.csv", row.names = F)`  

