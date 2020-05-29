#!/export/personal/lijj/0.temp/2.software/R-patched/bin/R  
#findGSE histofile sizek outdir
args <- commandArgs(trailingOnly=TRUE)
cw=getwd()
his <-  args[1]
k<- as.numeric(args[2])
out<-  args[3]
library("findGSE")
findGSE(histo=his,sizek=k,outdir=out)
