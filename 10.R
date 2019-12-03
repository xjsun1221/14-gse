rm(list = ls())
ni = 10 #GSE70529,配对4组,太难不看系列

{options(stringsAsFactors = F)
  library(GEOquery)
  load("bg.Rdata")
  eSet <- bg[[ni]]
  exp <- exprs(eSet[[1]]) 
  exp[1:4,1:4]
  pd <- pData(eSet[[1]])
  anno = eSet[[1]]@annotation
  if(!require(stringr))install.packages("stringr")
  library(stringr)
  gse = names(bg)[ni]
  print(exp[1:4,1:4])
  View(pd)
}

identical(rownames(pd),colnames(exp))

group_list = case_when(str_detect(pd$title,"A")~"A",
                       str_detect(pd$title,"B")~"B",
                       str_detect(pd$title,"C")~"C",
                       str_detect(pd$title,"D")~"D")

#exp = log2(exp+1)
#PCA
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)#画主成分分析图需要加载这两个包
  library(factoextra) 
  # pca的统一操作走起
  dat.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  dir.create(names(bg)[ni])
  ggsave(paste0(names(bg)[ni],"/PCA.png"))
}


#热图 
cg=names(tail(sort(apply(exp,1,sd)),1000))
if(!require(pheatmap))install.packages("pheatmap")
library(pheatmap)
n=exp[cg,]

#加注释分组

annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row")

exp[1:4,1:4]
boxplot(exp[1,]~group_list) 


