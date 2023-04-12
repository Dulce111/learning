# 3.28 Mfuzz分析学习
## 1 install packages
`if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")`

`BiocManager::install("Mfuzz")`
`browseVignettes("Mfuzz")`
`library(Mfuzz)`

## 2 下载示例数据进行学习分析
 - 内容来自：[CSDN-生信宝典 ](https://blog.csdn.net/qazplm12_3/article/details/117341589#:~:text=Mfuzz%E5%8C%85%E6%9C%80%E5%88%9D%E6%98%AF%E4%B8%BA%E5%A4%84%E7%90%86%E5%9F%BA%E5%9B%A0%E8%A1%A8%E8%BE%BE%E6%88%96%E8%9B%8B%E7%99%BD%E8%A1%A8%E8%BE%BE%E8%B0%B1%E6%95%B0%E6%8D%AE%E8%80%8C%E5%BC%80%E5%8F%91%E7%9A%84%E4%B8%80%E7%A7%8D%E8%81%9A%E7%B1%BB%E6%96%B9%E6%B3%95%EF%BC%8C%E6%A0%B8%E5%BF%83%E7%AE%97%E6%B3%95%E5%9F%BA%E4%BA%8E%20%E6%A8%A1%E7%B3%8Ac%E5%9D%87%E5%80%BC%E8%81%9A%E7%B1%BB%EF%BC%88Fuzzy,C-Means%20Clustering%EF%BC%8CFCM%EF%BC%89%20%EF%BC%8C%E7%94%A8%E4%BA%8E%E5%9C%A8%E5%85%B7%E6%9C%89%E6%97%B6%E9%97%B4%E5%BA%8F%E5%88%97%E7%89%B9%E5%BE%81%E7%9A%84%E8%BD%AC%E5%BD%95%E7%BB%84%E3%80%81%E8%9B%8B%E7%99%BD%E8%B4%A8%E7%BB%84%E6%95%B0%E6%8D%AE%E4%B8%AD%E5%88%86%E6%9E%90%E5%9F%BA%E5%9B%A0%E6%88%96%E8%9B%8B%E7%99%BD%E8%A1%A8%E8%BE%BE%E7%9A%84%E6%97%B6%E9%97%B4%E8%B6%8B%E5%8A%BF%EF%BC%8C%E5%B9%B6%E5%B0%86%E5%85%B7%E6%9C%89%E7%9B%B8%E4%BC%BC%E8%A1%A8%E8%BE%BE%E6%A8%A1%E5%BC%8F%E7%9A%84%E5%9F%BA%E5%9B%A0%E6%88%96%E8%9B%8B%E7%99%BD%E5%88%92%E5%88%86%E8%81%9A%E7%B1%BB%EF%BC%8C%E5%B8%AE%E5%8A%A9%E4%BA%86%E8%A7%A3%E8%BF%99%E4%BA%9B%E7%94%9F%E7%89%A9%E5%AD%A6%E5%88%86%E5%AD%90%E7%9A%84%E5%8A%A8%E6%80%81%E6%A8%A1%E5%BC%8F%E4%BB%A5%E5%8F%8A%E4%B8%8E%E5%8A%9F%E8%83%BD%E7%9A%84%E8%81%94%E7%B3%BB%E3%80%82) 


`setwd("D:/test_mfuzz_clustering")`  
`embroy_pro <- read.delim("mmc2.union_all_protein_exp.txt", row.names = 1, check.names = F) `
`embroy_pro <- as.matrix(embroy_pro)`

### step1 构建对象
`mfuzz_class <- new('ExpressionSet',exprs = embroy_pro)`
### step2 处理缺失or异常值
`mfuzz_class <- filter.NA(mfuzz_class,thres = 0.25)`   
`mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')`   
`mfuzz_class <- filter.std(mfuzz_class, min.std = 0)`

### step3 标准化数据
`mfuzz_class <- standardise(mfuzz_class)`  
>  手动定义目标聚类的个数，根据文献设为10 ，即，期望获得10个聚类
 
 - 设定随机数种子  
`set.seed(123)`  
`cluster_num <- 10`   
`mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))`

 - 作图参考 ？`mfuzz.plot2`    

`mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = colnames(embroy_pro))`
![mfuzz_result](./pictures/mfuzz_result.png)

**Q:如何获得各聚类群中包含的蛋白质？**  
**A: 查看每个聚类群中各自包含的蛋白数量** 

`cluster_size <- mfuzz_cluster$size`    
`names(cluster_size) <- 1:cluster_num`    
`cluster_size`

### step4 查看每个蛋白所属的聚类群
`head(mfuzz_cluster$cluster)`

> Mfuzz通过计算一个叫***membership***的统计量判断蛋白质所属的聚类群，以最大的membership值为准
 - 查看membership值   
`head(mfuzz_cluster$membership)`

### step5 提取所有蛋白质所属的聚类群，并和原始表达值整合在一起

`protein_cluster <- mfuzz_cluster$cluster`    
`protein_cluster <- cbind(embroy_pro[names(protein_cluster), ], protein_cluster)`   
`head(protein_cluster)`      
`write.table(protein_cluster, 'output_protein_cluster.txt', sep ='\t', col.names = NA, quote = F)`

###  若想提取数据分析过程中标准化后的表达值（作图值而非原始值）
`protein_cluster <- mfuzz_cluster$cluster`   
`protein_standard <- mfuzz_class@assayData$exprs`   
`protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)`   
 - 报错：*Error in cbind(...) : number of rows of matrices must match (see arg 2)*   
 - 查看数据行数和列数，看哪一步出问题

### step6 final analysis 挑选感兴趣的时间表达特征的聚类群，即可在表中进一步分析一个群中蛋白质的功能等
### · other questions
- 划分多少个群合适？   
 - 手动尝试，若出现折线存在某种时间趋势，则表明存在冗余的聚类群，可能需要减少群数量，若一个群内波动过大，则可能出现了多个密度中心

**Reference** ：
*Gao Yawei, Liu Xiaoyu, Tang Bin et al. Protein Expression Landscape of Mouse Embryos during Pre-implantation Development. Cell Rep, 2017, 21: 3957-3969.*
###


# 3.28 下午
# learning ***cluster packages*** from [CSDN](https://blog.csdn.net/woodcorpse/article/details/79274461)  
    
    
 - 数据来源：R内置数据集 *USArrests*
`install.packages("factoexra")`    
`library('factoextra')`   
`library('ggplot2')`   
`data("USArrests")`   
`library("cluster")`   
### 1 去掉缺失值
`USArrests = na.omit(USArrests)`    
`head(USArrests, n = 6)`   
 - 一些必要的数据检查，均值 标准差等  
`desc_stats = data.frame(Min=apply(USArrests, 2, min),
                        Med = apply(USArrests, 2, median),         
                        Mean = apply(USArrests, 2,mean),
                        SD = apply(USArrests, 2, sd),
                        Max = apply(USArrests, 2, max))`   
 - 保留小数点后一位  
`desc_stats = round(desc_stats, 1)`   
`desc_stats`   
### 2 数据标准化及其评估  
##方差和均值变化很大时进行标准化  
`df = scale(USArrests)`  
##数据集群性评估？？？  
`res = get_clust_tendency(df, 40, graph = T)`
`res$hopkins_stat`  
 - ##result>0.5 不聚合？  
`res$plot`

### 3 (Plus) 用函数clusGap估计最优聚合数，函数`fviz_gap_stat()`用于可视化
`set.seed(123)`   
`gap_stat = clusGap(df, FUN = kmeans, nstart = 25, K.max = 10, B = 500)`   
 - 报错：*Error in clusGap(df, FUN = kmeans, nstart = 25, K.max = 10, B = 500) : could not find function "clusGap"*  
 - reason: 没有加载cluster包  
##`library("cluster")`后解决

`fviz_gap_stat(gap_stat)`#可视化
![figure2](./pictures/gap_stat.png)

 - 使用kmeans进行聚类，选择25个随机集（怎么确定25个集？）  
`km.res = kmeans(df, 4, nstart = 25)`    
`fviz_cluster(km.res, USArrests)`   
 - 提取聚类轮廓图
`sil = silhouette(km.res$cluster, dist(df))`    
`rownames(sil) = rownames(USArrests)`   
`head(sil[, 1:3])`   
 ![聚类轮廓图]
- visualize  
`fviz_silhouette(sil)`
![figure3](./pictures/sil_cluster.png)
- 图中出现负值，用函数 *silhouette* 确定是哪个观测值
`neg_sil_index = which(sil[, "sil_width"] < 0)`   
`sil[neg_sil_index, , drop = F]`


### 4 **eclust()*增强的聚类分析--简化聚类分析的工作流程，可计算层次聚类和分区聚类，自动计算最佳聚类簇数
 - calculate k-means
`res.km = eclust(df, "kmeans")`   
- gap statistic plot 
`fviz_gap_stat(res.km$gap_stat)`  
![figure4](./pictures/gap_stat.png)
 - enhanced hierarchical clustering  
`res.hc = eclust(df, "hclust")`   
![figure5](./pictures/eclust_final.png)
`fviz_dend(res.hc, rect = T)`  
- 层级聚类  
`fviz_silhouette(res.hc)`  
![figure6](./pictures/eclust.png) 
`fviz_cluster(res.hc)`
![figure6](./pictures/res.hc_scatterplot.png)
##分析环境相关信息？？？
