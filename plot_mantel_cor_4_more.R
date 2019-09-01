##更新内容，可以设置是否添加连线，去除某些不相关的连线


######环境变量之间的关系··································································


# title = c("A","B","c")
plot_mantel_cor_3_more = function(env = env.st,report = report,title = "16S_microbiology",sub1 = NULL,sub2 = NULL,sub3 = NULL,sub4 = NULL){
  colnames(report) = c(colnames(report)[1],paste("R",1:(dim(report) [2] - 1),sep = ""))
  report
  library(ggcorrplot)
  library(igraph)
  library(psych)
  
  occor = corr.test(env.st,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
  occor.r = occor$r # 取相关性矩阵R值
  
  occor.p = occor$p # 取相关性矩阵p值
  occor.r[occor.p>0.05] = 0
  
  # p = ggcorrplot(occor.r,p.mat = occor.p, method = "circle",lab = TRUE,insig = "blank",outline.color	= "white",
  #                ggtheme = ggplot2::theme_bw())
  # 
  # p
  # 
  ###############
  
  
  
  # write.csv(occor.r,"./occor.r.csv",quote = F)
  library(reshape2)
  
  occor.r2 = occor.r[lower.tri(occor.r, diag = TRUE)]
  # occor.r2 = occor.r[lower.tri(occor.r, diag = TRUE)]
  length(colnames(occor.r))
  # a = rep(0,length(occor.r2))
  # i = 1
  ### 构造环境变量三角矩阵图
  
  ## dd：构建三交矩阵横坐标
  cc = rep(length(colnames(occor.r)),(length(colnames(occor.r))))
  for (i in 1:(length(colnames(occor.r))-1)) {
    ee = rep(length(colnames(occor.r))-i,(length(colnames(occor.r))-i))
    dd = c(cc,ee)
    
    cc = dd
  }
  dd = (length(colnames(occor.r))+1) - dd
  
  
  # cc = rep(1:length(colnames(occor.r)))
  # for (i in 1:(length(colnames(occor.r))-1)) {
  #   ee = c(1:(length(colnames(occor.r))-i))
  #   rr = c(cc,ee)
  #   
  #   cc = rr
  # }
  # 
  # rr
  ##ww:构建三角矩阵纵坐标
  gg =  rep(length(colnames(occor.r)):1)
  for (i in 1:(length(colnames(occor.r))-1)) {
    ee = c((length(colnames(occor.r))-i):1)
    ww = c(gg,ee)
    
    gg = ww
  }
  
  ww
  ##构造变量对应关系
  wwt =data.frame(x = dd,y = ww)
  ##提取相关值大小，这是下三角矩阵转化
  # xa = melt(occor.r)
  wwt$mantelR = occor.r2
  wwt
  
  
  ### 下面添加群落和环境因子的相关结果
  
  ##构造每个环境因子的坐标
  
  x = c(1:(length(colnames(occor.r))) )
  y = c((length(colnames(occor.r))+1) :2)
  
  data2 = data.frame(x = x,y = y)
  # 
  data2 = cbind(data2,as.data.frame(report[2:dim(report)[2]]^2*500))
  data2

  axa = as.matrix(report[2:dim(report)[2]])
  axa[axa>0] <- 1
  axa[axa<0] <- -1
  axa[axa==0] <- 0
  colnames(axa) = paste(colnames(axa),"d",sep = "")
  data2 = cbind(data2,axa)
  data2


  data2$label = colnames(env.st)
  aaa = data2

  library(tidyverse)
  p = ggplot() + 
    geom_tile(aes(x = x, y = y),wwt,fill=NA,color='gray',size=0.5)+
    geom_point(aes(x = x, y = y,size=mantelR,fill=mantelR),wwt, shape=22,color='white')+
    scale_size(range = c(1, 8))+
    scale_fill_distiller(palette="RdYlBu")
  p
  
 
  ##如果三角之外只有一个点
  
  if (dim(report) [2] == 2) {
    
    
    data3<- filter(data2, label %in% sub)
    data3
    p = p+
      geom_curve(aes(x = max(wwt$x)*3/4, y = max(wwt$y)*3/4, xend = x, yend = y,group = as.factor(data3$R1d),color = as.factor(data3$R1d)),curvature = 0.2,data3,size = data3$R1) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,aaa,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*3/4, y = max(wwt$y)*3/4),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    
    
    
    
  }
  ## 有两个群落
  if (dim(report) [2] == 3) {
    
    data3<- filter(data2, label %in% sub1)
    data3
    
    p = p+
      geom_curve(aes(x = max(wwt$x)*2/3, y = max(wwt$y), xend = x, yend = y,group = as.factor(data3$R1d),color = as.factor(data3$R1d)),curvature = 0.2,data3,size = data3$R1) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,data3,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*2/3, y = max(wwt$y)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    data4<- filter(data2, label %in% sub2)
    data4
    p = p+
      geom_curve(aes(x = max(wwt$x)*4.2/5, y = max(wwt$y)*3/4, xend = x, yend = y,group = as.factor(data4$R2d),color = as.factor(data4$R2d)),curvature = 0.2,data4,size = data4$R2) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,aaa,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*4.2/5, y = max(wwt$y)*3/4),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    
  }
  p
  ## 有三个群落
  if (dim(report) [2] == 4) {
    data3<- filter(data2, label %in% sub1)
    data3
    p = p+
      geom_curve(aes(x = max(wwt$x)*1.5/3, y = max(wwt$y*1.5), xend = x, yend = y,group = as.factor(data3$R1d),color = as.factor(data3$R1d)),curvature = 0.2,data3,size = data3$R1) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,data3,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*1.5/3, y = max(wwt$y*1.5)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    p
    data4<- filter(data2, label %in% sub2)
    data4
    p = p+
      geom_curve(aes(x = max(wwt$x)*3.5/4, y = max(wwt$y), xend = x, yend = y,group = as.factor(data4$R2d),color = as.factor(data4$R2d)),curvature = 0.2,data4,size = data4$R2) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,data4,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*3.5/4, y = max(wwt$y)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    data5<- filter(data2, label %in% sub3)
    data5
    
    p = p+
      geom_curve(aes(x = max(wwt$x)*4.8/5, y = max(wwt$y)*2/4, xend = x, yend = y,group = as.factor(data5$R3d),color = as.factor(data5$R3d)),curvature = 0.2,data5,size = data5$R3) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,aaa,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*4.8/5, y = max(wwt$y)*2/4),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    
    
  }
  ## 有4个群落
  if (dim(report) [2] == 5) {
    data3<- filter(data2, label %in% sub1)
    data3
    p = p+
      geom_curve(aes(x = 2, y = max(wwt$y*1.5), xend = x, yend = y,group = as.factor(data3$R1d),color = as.factor(data3$R1d)),curvature = 0.2,data3,size = data3$R1) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,aaa,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = 2, y = max(wwt$y*1.5)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    p
    data4<- filter(data2, label %in% sub2)
    data4
    p = p+
      geom_curve(aes(x = 5, y = max(wwt$y*1.4), xend = x, yend = y,group = as.factor(data4$R2d),color = as.factor(data4$R2d)),curvature = 0.2,data4,size = data4$R2) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,data4,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = 5, y = max(wwt$y*1.4)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    p
    data5<- filter(data2, label %in% sub3)
    data5
    p = p+
      geom_curve(aes(x = max(wwt$x)*3.5/4, y = max(wwt$y), xend = x, yend = y,group = as.factor(data5$R3d),color = as.factor(data5$R3d)),curvature = 0.2,data5,size = data5$R3) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,data5,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*3.5/4, y = max(wwt$y)),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    
    p
    data6<- filter(data2, label %in% sub3)
    data6
    p = p+
      geom_curve(aes(x = max(wwt$x)*4.8/5, y = max(wwt$y)*2/4, xend = x, yend = y,group = as.factor(data6$R4d),color = as.factor(data6$R4d)),curvature = 0.2,data6,size = data6$R4) +
      geom_point(aes(x = x, y = y),pch = 21,size =4,aaa,color = "black",fill = "#FFF5EB")+
      geom_point(aes(x = max(wwt$x)*4.8/5, y = max(wwt$y)*2/4),pch = 21,size = 6,color = "black",fill = "#FEE6CE")
    
    p
  }
  
  
  
  ### 添加环境因子添加标签
  p = p + geom_text(aes(x = x, y = 0,label= colnames(env.st)),size=4,data2) +
    geom_text(aes(x = 0, y = y-1,label= colnames(env.st)),size=4,data2) 
  p
  p = p + geom_text(aes(x = x+0.5, y = y+0.5,label=data2$label),size=4,data2)
  p
  ##添加群落矩阵标签
  if (dim(report) [2] == 2) {
    p = p +
      geom_text(aes(x = max(wwt$x)*3/4+1, y = max(wwt$y)*3/4+1,label= title),size = 6) 
    p
    
  }
  # title = c("A","B","c")
  if (dim(report) [2] == 3) {
    p = p +
      geom_text(aes(x = max(wwt$x)*2/3 +1, y = max(wwt$y),label= title[1]),size = 6) +
      geom_text(aes(x = max(wwt$x)*4.2/5 +1, y = max(wwt$y)*3/4,label= title[2]),size = 6) 
    
    p
    
  }
  if (dim(report) [2] == 4) {
    p = p +
      geom_text(aes(x = max(wwt$x)*1.5/3+1, y = max(wwt$y*1.5)+1,label= title[1]),size = 6) +
      geom_text(aes(x = max(wwt$x)*3.5/4+1, y = max(wwt$y)+1,label= title[2]),size = 6) +
      geom_text(aes(x = max(wwt$x)*4.8/5+1, y = max(wwt$y)*2/4+1,label= title[3]),size = 6) 
    
    p
    
  }
  #四个群落
  if (dim(report) [2] == 5) {
    p = p +
      geom_text(aes(x = 2+1, y = max(wwt$y*1.5)+1,label= title[1]),size = 6) +
      geom_text(aes(x = 5+1, y = max(wwt$y*1.4)+1,label= title[2]),size = 6) +
      geom_text(aes(x = max(wwt$x)*3.5/4+1, y = max(wwt$y)+1,label= title[3]),size = 6) +
      geom_text(aes(x = max(wwt$x)*4.8/5+1, y = max(wwt$y)*2/4+1,label= title[4]),size = 6) 
    
    p
    
  }
  
  
  
  
  
  ### 修改主题
  p = p +theme_void()#横纵坐标群去掉
  p
  
  
  
}


# 
# plotnamea = paste(path,"New_matel_cor_plot.pdf",sep = "")
# ggsave(plotnamea, p, width = 8, height = 6)
# 
# 
# 
# data2
# # ggplot() + 
# #   geom_tile(aes(x = x, y = y),wwt,fill=NA,color='gray',size=1)+
# #   geom_point(aes(x = x, y = y,size=mantelR,fill=mantelR),wwt, shape=22,color='white')+
# #   scale_size(range = c(1, 8))+
# #   scale_fill_distiller(palette="RdYlBu")+
# #   geom_curve(aes(x = 10, y = 10, xend = x, yend = y),curvature = 0.2,data2) +
# #   geom_point(aes(x = x, y = y),pch = 21) +
# #   geom_point(aes(x = 10, y = 10),pch = 21,size = 4)
# 
# ###添加线的粗细和颜色映射
# 
# # write.csv(report,"./report.csv",quote = F)
# 
# data2$count = report$R
# data2$count = as.character(data2$count)
# data2$count = as.numeric(data2$count)
# data2
# 
# 
# as = rep("a",length(data2$count))
# 
# 
# for (i in 1:length(data2$count)) {
#   if (data2$count[i] > 0) {
#     as[i] = "+"
#   }
#   
#   if (data2$count[i] < 0) {
#     as[i] = "-"
#   }
#   if (data2$count[i] == 0) {
#     as[i] = "-"
#   }
# }
# as
# 
# data2$group = as
# data2$label = colnames(env.st)
# wwt
# 
# 
# data2$count1 = data2$count^2*500
