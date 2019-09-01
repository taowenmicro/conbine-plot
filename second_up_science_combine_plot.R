library("phyloseq")
library(microbiomeSeq)
library("vegan")
library("grid")
library("gridExtra")
library("ggplot2")
ps = readRDS("./a3_DADA2_table//ps.rds")
ps1 = ps
ps1 = filter_taxa(ps1, function(x) sum(x ) > 200 , TRUE);ps1



path = "./phyloseq_5_RDA_CCA_cor/"
dir.create(path)

vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}
otu = as.data.frame(t(vegan_otu(ps1)))

mapping = as.data.frame( sample_data(ps1))
env.dat = mapping[,3:ncol(sample_data(ps1))]

env.st = decostand(env.dat, method="standardize", MARGIN=2)#



env_dat = env.st

otu2 = otu
env_dat2 = env_dat

##OTU distance
BC.beta = vegdist(t(otu2), method="bray")
JC.beta = vegdist(t(otu2), method="jaccard",binary=T)

##env distance
env.std = decostand(env_dat2, method = "standardize", MARGIN=2) #按列标准化到均值等于0，标准偏差等于1.


#calculate the relationship between community and each env_factor============

report =c()
for(i in 1:ncol(env.std)){     
  envdis =vegdist(env.std[,i],method = "euclidean", na.rm=T)
  mantel.BC = mantel(envdis, BC.beta, na.rm=T)
  mantel.JC = mantel(envdis, JC.beta, na.rm=T)
  report = rbind(report,c(colnames(env_dat2)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif))
}

colnames(report)<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
report


#使用偏相关
#Partial Mantel test分析，可以将某些环境因子作为控制矩阵，计算微生物群落矩阵与剩余目标环境因子矩阵的相关性。即排除不关注的环境因子，研究关注的环境因子与微生物之间的关系。
#
report = c()
#i = 4
for(i in 1:ncol(env.std)){   ##
  envdis = dist(env.std[,i])
  envdis2 = dist(env.std[,-i])
  pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
  pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)
  report = rbind(report,c(colnames(env_dat2)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif))
}
colnames(report)<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
report

report = as.matrix(report)
report
report = as.data.frame(report)

report
report[,2:dim(report)[2]]<-lapply(report[,2:dim(report)[2]],as.character)
report[,2:dim(report)[2]]<-lapply(report[,2:dim(report)[2]],as.numeric)

head(report)


### 下面进行筛选显著的关系,如果只有一组对象，很好做，直接筛选表格就够了，但是如果有多个组就很难办了

#首先我们使用BARY距离计算的mantel结果进行计算
report1 = report[1:3]
head(report1)
##设置显著性阈值
alpha = 0.3
library(tidyverse)
report1<- filter(report1, p.BC < alpha)
sub = report1$Envs
report1 = report[1:2]

p = plot_mantel_cor_3_more(env = env.st,report = report1,title = c("A") )
p

## 那么多个可怎么办？



source("./plot_mantel_cor_4_more.R")

report = read.csv("./report.csv",row.names = 1)
sub1 = c("env7","env8","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A"),sub1 = sub1 )
p


report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
# report$bb = report$R
head(report)
##下面我们指定每个矩阵需要和那些env相连接
# 注意ggplot不可以使用使用不同的数据框，但是相同的数据狂名称
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B"),sub1 = sub1,sub2 = sub2 )
p


report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
report$bb = report$R
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
sub3 = c("env11","env6","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B","C"),sub1 = sub1,sub2 = sub2 ,sub3 = sub3 )
p

report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
report$bb = report$R
report$cc = report$R

report
##下面我们指定每个矩阵需要和那些env相连接
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
sub3 = c("env11","env6","env4")
sub4 = c("env7","env6","env10","env12")

p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B","C","D") ,sub1 = sub1,sub2 = sub2 ,sub3 = sub3,sub4 = sub4)
p








source("./plot_mantel_cor_4_2_more.R")

report = read.csv("./report.csv",row.names = 1)
sub1 = c("env7","env8","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A"),sub1 = sub1 )
p


report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
# report$bb = report$R
head(report)
##下面我们指定每个矩阵需要和那些env相连接
# 注意ggplot不可以使用使用不同的数据框，但是相同的数据狂名称
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B"),sub1 = sub1,sub2 = sub2 )
p


report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
report$bb = report$R
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
sub3 = c("env11","env6","env4")
p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B","C"),sub1 = sub1,sub2 = sub2 ,sub3 = sub3 )
p

report = read.csv("./report.csv",row.names = 1)
report$aa = report$R
report$bb = report$R
report$cc = report$R

report
##下面我们指定每个矩阵需要和那些env相连接
sub1 = c("env7","env8","env4")
sub2 = c("env5","env10","env4")
sub3 = c("env11","env6","env4")
sub4 = c("env7","env6","env10","env12")

p = plot_mantel_cor_3_more(env = env.st,report = report,title = c("A","B","C","D") ,sub1 = sub1,sub2 = sub2 ,sub3 = sub3,sub4 = sub4)
p































