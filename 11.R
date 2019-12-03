rm(list = ls())
ni = 11 #GSE70529,配对4组,太难不看系列

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

group_list = case_when(str_detect(pd$title,"non-del11q")~"non-del11q",
                       TRUE~"del11q")

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

#差异分析，用limma包来做
#需要表达矩阵和group_list，其他都不要动
library(limma)
design=model.matrix(~factor(group_list,levels = c("non-del11q","del11q")))
fit=lmFit(exp,design)
fit=eBayes(fit)

#差异基因排名
deg=topTable(fit,coef=2,number = Inf)
head(deg)
boxplot(exp[rownames(deg)[1],]~group_list)

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
#tibble::rownames_to_column()
head(deg)


#2.加symbol列，火山图要用
#id转换，查找芯片平台对应的包
eSet[[1]]@annotation
#http://www.bio-info-trainee.com/1399.html
#hgu133a
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)

#merge
deg <- inner_join(deg,ids,by="probe_id")
head(deg)
#3.加change列：上调或下调，火山图要用

logFC_t=1 #不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
change=ifelse(deg$P.Value>0.01,'stable', 
              ifelse( deg$logFC >logFC_t,'up', 
                      ifelse( deg$logFC < -logFC_t,'down','stable') )
)
deg <- mutate(deg,change)
head(deg)
table(deg$change)
deg <- mutate(deg,v = -log10(P.Value))

#4.加ENTREZID列，后面富集分析要用
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(unique(deg$symbol), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

head(deg)

library(dplyr)
dat <- mutate(deg,v=-log10(P.Value))
head(dat)
p <- ggplot(data = deg, 
            aes(x = logFC, 
                y = v)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()
for_label <- deg %>% 
  filter(abs(logFC) >2& P.Value< 0.01)
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
ggsave(paste0(names(bg)[ni],"/volcano.png"))


x=deg$logFC 
names(x)=deg$probe_id 
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))

library(pheatmap)
n=exp[cg,]

annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 

pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         #cluster_cols = F, 
         annotation_col=annotation_col) 
dev.off()
#保存

pdf(file = "heatmap.pdf")
pheatmap(n,show_colnames =F,
         show_rownames = F,
         scale = "row",
         #cluster_cols = F, 
         annotation_col=annotation_col) 
dev.off()
file.copy("heatmap.pdf",paste0(names(bg)[ni],"/heatmap.pdf"))
file.remove("heatmap.pdf")

