# E-MEXP-1422（找不到探针id）/  GSE474/
# E-MTAB-3017（找不到探针id）/  GSE58979-NASH（从临床信息中找不到有用的分组信息/
# GSE1462/                      GSE59045/
# GSE18732（探针id转换问题）/   GSE60291/
# GSE20950/                     GSE62832（不显著，可能需要重新分组）/
# GSE21785/                     GSE70529（可能需要重新分组）/
# GSE26526/                     GSE72158（FC太小）/
# GSE32575（阈值需要调整）/     illuminaHGv4/
# GSE43837/                     MsigDB/

#0写个循环下数据
x = "E-MEXP-1422（找不到探针id）/  GSE474/
E-MTAB-3017（找不到探针id）/  GSE58979-NASH（从临床信息中找不到有用的分组信息/
GSE1462/                      GSE59045/
GSE18732（探针id转换问题）/   GSE60291/
GSE20950/                     GSE62832（不显著，可能需要重新分组）/
GSE21785/                     GSE70529（可能需要重新分组）/
GSE26526/                     GSE72158（FC太小）/
GSE32575（阈值需要调整）/     illuminaHGv4/
GSE43837/                     MsigDB/"
library(stringr)
x2 = str_split(x,"/")%>% unlist()
x2
x2 = str_replace_all(x2,c("\n"," "),"")
x2
x2 = x2[str_starts(x2,"GSE")]
x2
x3 = str_split(x2,"（",simplify = T)
x4 = as.character(x3[,1])
x4
x4 = str_replace(x4,"-NASH","")
x4
library(GEOquery)
bg = list()
for (i in 1:length(x4)){
  bg[[i]] <- getGEO(x4[[i]], #系列编号
               destdir = '.', #当前目录
               getGPL = F)
}
names(bg) = x4
bg
dir(pattern = ".gz")
save(bg,file = "bg.Rdata")

