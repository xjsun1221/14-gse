rm(list = ls())
ni = 1 #GSE474，3组,完美

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

group_list=ifelse(str_detect(pd$title,"MObese"),"MObese",
                  ifelse(str_detect(pd$title,"NonObese"),"NonObese","Obese"))
group_list

exp = log2(exp+1)
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


#两两组合
group_list = as.character(group_list)
colnames(exp)
x = unique(group_list)
x
cp = list()
nor = 1
for(i in 1:length(x)){
  if(i == nor) next
  cp[[i]] = cbind(exp[,group_list ==x[[nor]]],
                  exp[,group_list ==x[[i]]])
  names(cp)[i] = x[[i]]
}
cp = cp[-nor]
i = 1
deg = list()
grp =list()
#批量差异分析
for(i in 1:length(cp)){
  exps = cp[[i]]
  grps = group_list[(group_list == x[nor])|(group_list== names(cp)[i])]
  grps = factor(grps,levels = c(x[nor],names(cp)[i]))
  grp[[i]] = grps
  library(limma)
  design=model.matrix(~grps)
  fit=lmFit(exps,design)
  fit=eBayes(fit)
  deg[[i]]=topTable(fit,coef=2,number = Inf)
  print(deg[[i]][1,1])
  boxplot(exps[rownames(deg[[i]])[1],]~grps)
}
#2.加symbol列，火山图要用
#id转换，查找芯片平台对应的包
anno
#http://www.bio-info-trainee.com/1399.html
#hgu133a
if(!require(hgu133a.db)) BiocManager::install("hgu133a.db")
library(hgu133a.db)
ids <- toTable(hgu133aSYMBOL)
head(ids)
library(dplyr)
i = 1
for(i in 1:length(deg)){
  #1.加probe_id列
  deg[[i]] <- mutate(deg[[i]],probe_id=rownames(deg[[i]]))
  #2.id转换
  deg[[i]] <- inner_join(deg[[i]],ids,by="probe_id")
  print(head(deg[[i]])) 
  logFC_t=1 
  #3.change
  {
    change=ifelse(deg[[i]] $P.Value>0.01,'stable', 
                  ifelse( deg[[i]]$logFC >logFC_t,'up', 
                          ifelse( deg[[i]]$logFC < -logFC_t,'down','stable') )
    )
    deg[[i]] <- mutate(deg[[i]],change)
    head(deg[[i]])
    print(table(deg[[i]]$change))
    deg[[i]] <- mutate(deg[[i]],v = -log10(P.Value))
  }
  #4.加ENTREZID列，后面富集分析要用
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  s2e <- bitr(unique(deg[[i]]$symbol), fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Hs.eg.db)
  head(s2e)
  head(deg[[i]])
  deg[[i]] <- inner_join(deg[[i]],s2e,by=c("symbol"="SYMBOL"))
  
  head(deg[[i]])
}
#批量火山图,批量热图
vo = function(x){
  p <- ggplot(data = x, 
              aes(x = logFC, 
                  y = v)) +
    geom_point(alpha=0.4, size=3.5, 
               aes(color=change)) +
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
    theme_bw()
  for_label <- x %>% 
    filter(abs(logFC) >4& P.Value< 0.00001)
  p +
    geom_point(size = 3, shape = 1, data = for_label) +
    ggrepel::geom_label_repel(
      aes(label = symbol),
      data = for_label,
      color="black"
    )
  ggsave(paste0(names(bg)[ni],"/volcano",k,".png"))
  
}

y = deg[[1]]
k=1
hp = function(y) {
  x = y$logFC
  names(x) = y$probe_id
  cg = c(names(head(sort(x), 100)),
         names(tail(sort(x), 100)))
  
  library(pheatmap)
  n = cp[[k]][cg, ]
  annotation_col = data.frame(group = grp[[k]])
  rownames(annotation_col) = colnames(n)
  pdf(file = paste0("heatmap", k, ".pdf"))
  test = pheatmap(
    n,
    show_colnames = F,
    show_rownames = F,
    scale = "row",
    #cluster_cols = F,
    annotation_col = annotation_col
  )
  print(test)
  dev.off()
  file.copy(paste0("heatmap", k, ".pdf"), paste0(names(bg)[ni], "/", paste0("heatmap", k, ".pdf")))
  file.remove(paste0("heatmap", k, ".pdf"))
}

for(k in 1:length(deg)){
  vo(deg[[k]])
}

for(k in 1:length(deg)){
  hp(deg[[k]])
}

