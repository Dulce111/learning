# learning bioinformation-  3.24 bioconduct_MSstats 



## 1 下载安装biocontruct
 - Install/Update sandpaper

>options(repos = c(carpentries = "https://carpentries.r-universe.dev/", 
                  CRAN = "https://cloud.r-project.org"))

> install.packages("sandpaper")



 - update the workflows 
> library("sandpaper")
update_github_workflows()

 - update my R and Rstudio

> install.packages('installr')    
> library('installr')

> updateR()

*报错: cannot find “string”*

> install.packages("stringi")
library(stringi)

 - install bioconduct
> if (!require("BiocManager", quietly = TRUE))
  
  >install.packages("BiocManager")
BiocManager::install(version = "3.16")
 - 测试bioconduct安装情况
> BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::available()


## 2 MS analysis
 - install MSstats
 > if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")

> browseVignettes("MSstats")



 - Attention! 
 
 'MSstatsInput.csv' is the MSstats report from [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view). 
 
 PS: skyline can analysis .raw data and output .csv 

> input <- read.csv(file="MSstatsInput.csv")

> raw <- SkylinetoMSstatsFormat(input)

 -  MSstats can also Read in MaxQuant files

eg:
> proteinGroups <- read.table("proteinGroups.txt", sep     ="\t", header=TRUE)

>  infile <- read.table("evidence.txt", sep="\t", header=TRUE)

>/#  Read in annotation including condition and biological replicates per run.                                       
>/# Users should make this annotation file. It is not the output from MaxQuant.  ??? **don't understand**

>annot <- read.csv("annotation.csv", header=TRUE)
   
> raw <- MaxQtoMSstatsFormat(evidence=infile,                          annotation=annot,
 proteinGroups=proteinGroups)
 
 - practice MSstat>>dataProcess with .raw data
> library("MSstats")  
dataProcess()


***error**:The following object is masked from ‘package:grDevices’:savePlot

### 2.1 try dataProcess

>dataProcess(
  "2003151_WHJ_CSF_TMT_11.raw", logTrans = 2,  
  normalization = "equalizeMedians",  
  nameStandards = NULL,  
  featureSubset = "all",  
  remove_uninformative_feature_outlier = FALSE,  
  min_feature_count = 2,  
  n_top_feature = 3,  
  summaryMethod = "TMP",  
  equalFeatureVar = TRUE,  
  censoredInt = "NA",  
  MBimpute = TRUE,  
  remove50missing = FALSE,  
  fix_missing = NULL,  
  maxQuantileforCensored = 0.999,  
  use_log_file = TRUE,  
  append = FALSE,  
  verbose = TRUE,  
  log_file_path = NULL  
  )

 - error: 
 >Error in (function (classes, fdef, mtable)  :   
函数‘.checkDataValidity’标签‘"character"’找不到继承方法
 - reason: 函数无法识别输入参数的类型，即 dataprocess无法识别.raw data

 - **try .txt**
 
> getwd()  
setwd("D:/test_example/combined/txt")    
proteinGroups <- read.table("proteinGroups.txt", sep="\t", header=TRUE)

>dataProcess("proteinGroups.txt",logTrans = 2)

 - error:函数‘.checkDataValidity’标签‘"character"’找不到继承方法



###Mfuzz
BiocManager::install("Mfuzz")
browseVignettes("Mfuzz")
n

# 3.24 继续bioconduct
 - 聚类分析直接用cluster包
> library("cluster")

> library("ggplot2")

update package
>old.packages()  
> update.packages()


















































