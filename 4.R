rm(list = ls())
ni = 4 #GSE59045,1组
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
