library(ggplot2)
library(cowplot)
library(ggpubr)
library(mise)
library(mgcv)

mise()
ReadCsv = function(df.result){
  df.result = t(df.result)
  df.result = df.result[c(2:nrow(df.result)),]
  return(df.result)
}

ComputeSD = function(df.result){
  df.result = t(df.result)
  df.result = df.result[c(2:nrow(df.result)),]
  SDarray = rep(0,ncol(df.result))
  for (colidx in 1:ncol(df.result)){
    SDarray[colidx] = sd(df.result[,colidx],na.rm=T)
  }
  return(SDarray)
}

computeColMeans = function(df.result){
  tmp = t(df.result)
  return(colMeans(tmp[c(2:nrow(tmp)),],na.rm=T))
}

MakeFileName = function(instancename,mismode,lossname,mymodel){
  instancename = paste(instancename,"mis",mismode,sep="")
  instancename = paste(instancename,'-',mymodel,"_",lossname,sep="")
  filename = paste(instancename,".csv",sep="")
  return(filename)
}

ConstructDFPlot = function(instancename,mismode,lossname,mean_median){
  # instancename = "planid-0816-1900"
  # mismode = 0
  # lossname = "weight"
  # instancename = paste(instancename,"mis",mismode,sep="")
  # instancename = paste(instancename,'-summary',"_",lossname,sep="")
  # filename = paste(instancename,".csv",sep="")
  
  df.result.summary = read.csv(MakeFileName(instancename,mismode,lossname,mymodel='summary'))
  df.result.DR = read.csv(MakeFileName(instancename,mismode,lossname,mymodel='DR'))
  df.result.PlugIn = read.csv(MakeFileName(instancename,mismode,lossname,mymodel='PlugIn'))
  df.result.WERM = read.csv(MakeFileName(instancename,mismode,lossname,mymodel='WERM'))
  
  confidence_coef = 1/2
  
  if (mean_median == "median"){
    DRCenter = df.result.summary$DR.50
    PlugInCenter = df.result.summary$PlugIn.50
    WERMCenter = df.result.summary$WERM.50
  }else{
    DRCenter = df.result.summary$DR.mean
    PlugInCenter = df.result.summary$PlugIn.mean
    WERMCenter = df.result.summary$WERM.mean
  }
  DRSD = ComputeSD(df.result.DR)
  DRLow = DRCenter - confidence_coef*DRSD
  DRHigh = DRCenter + confidence_coef*DRSD
  
  PlugInSD = ComputeSD(df.result.PlugIn)
  PlugInLow = PlugInCenter - confidence_coef*PlugInSD
  PlugInHigh = PlugInCenter + confidence_coef*PlugInSD
  
  WERMSD = ComputeSD(df.result.WERM)
  WERMLow = WERMCenter - confidence_coef*WERMSD
  WERMHigh = WERMCenter + confidence_coef*WERMSD
  
  df.Plot = data.frame(Nlist = df.result.summary$Nlist, 
                       DR.50 = DRCenter,
                       DR.25 = DRLow,
                       DR.75 = DRHigh,
                       PlugIn.50 = PlugInCenter,
                       PlugIn.25 = PlugInLow,
                       PlugIn.75 = PlugInHigh,
                       WERM.50 = WERMCenter,
                       WERM.25 = WERMLow,
                       WERM.75 = WERMHigh
  )
  return(df.Plot)
}

ConstructBoxPlot = function(instancename,lossname){
  for (mismode in c(0,1,2)){
    df.result.DR = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='DR')))
    vector.DR = df.result.DR[,ncol(df.result.DR)]
    
    df.result.PlugIn = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='PlugIn')))
    vector.PlugIn = df.result.PlugIn[,ncol(df.result.PlugIn)]
    
    df.result.WERM = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='WERM')))
    vector.WERM = df.result.WERM[,ncol(df.result.WERM)]
    
    assign(paste('vector.DR.mis',mismode,sep=""),vector.DR)
    assign(paste('vector.PlugIn.mis',mismode,sep=""),vector.PlugIn)
    assign(paste('vector.WERM.mis',mismode,sep=""),vector.WERM)
  }
  
  myN = length(vector.PlugIn.mis0)
  AAE = c(vector.PlugIn.mis0,vector.PlugIn.mis1,vector.PlugIn.mis2,
          vector.WERM.mis0,vector.WERM.mis1,vector.WERM.mis2,
          vector.DR.mis0,vector.DR.mis1,vector.DR.mis2)
  label = c(rep("PI",myN),rep("PI",myN),rep("PI",myN),
            rep("WERM",myN),rep("WERM",myN),rep("WERM",myN),
            rep("DR",myN),rep("DR",myN),rep("DR",myN))
  Type = c(rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN))
  plotData = data.frame(AAE,label,Type)
  plotData$label = factor(plotData$label,c("PI","WERM","DR"))
  plotData$mycolor = c(rep("red",3*myN),rep("green",3*myN),rep("blue",3*myN))
  return(plotData)
}

ConstructDrawGG = function(instancename,mismode,lossname,ylimits,mean_median){
  df.result = ConstructDFPlot(instancename,mismode,lossname,mean_median)
  
  # General 
  regmethod = 'auto'
  # xlimits = c(0,5000)
  xlimits = c(0,max(df.result$Nlist))
  spanval = 1
  point_size = 3
  alpha_point = 1
  
  twoD = T
  medianTF = T
  
  gg = ggplot(data = df.result, aes(x=Nlist))
  
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=DR.50,colour="DR.50"),size=1,method=regmethod,se=F, span=spanval)  
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1.5,method=regmethod,se=F, span=spanval,linetype='dashed')
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=WERM.50,colour="WERM.50"),size=3,method=regmethod,se=F, span=spanval,linetype='dotted')
    
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=DR.25,ymax=DR.75),alpha=0.2, fill="dodgerblue")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=WERM.25,ymax=WERM.75),alpha=0.2, fill="orange")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=PlugIn.25,ymax=PlugIn.75),alpha=0.2, fill="firebrick1")

  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=point_size,alpha=alpha_point,shape=4)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=point_size,alpha=alpha_point,shape=16)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=WERM.50,colour='WERM.50'),size=point_size,alpha=alpha_point,shape=9)

  gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","WERM.50"),values = c("blue", "firebrick2","orange"),labels=c("DR-R","Plug-In","WERM"))
  
  # gg = gg + ggtitle(Dtitle) 
  gg = gg + coord_cartesian(ylim=ylimits)
  gg = gg + theme_bw()
  gg = gg + scale_x_continuous(name = "m", limits=xlimits) + scale_y_continuous(name = "MAAE")
  gg = gg + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                  axis.line.y = element_line(size = 0.5, colour = "black"),
                  axis.line = element_line(size=1, colour = "black"),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.text = element_text(size=15),
                  plot.title=element_text(size = 40),
                  axis.text.x = element_text(size = 20),
                  axis.text.y = element_text(size = 20),
                  axis.title.y = element_text(size=25),
                  axis.title.x = element_text(size=30)
  )
  return(gg)
}

ConstructBoxGG = function(instancename,lossname,ylimits){
  box.result = ConstructBoxPlot(instancename,lossname)
  ggbox = ggplot(box.result,aes(x=Type, y=AAE))
  ggbox = ggbox + geom_boxplot(aes(fill=label))
  # gg = gg + stat_summary(fun=mean, geom="point", aes(group=label), position=position_dodge(.9), color="red", size=3)
  # gg = gg + geom_errorbar(aes(x=Type,ymax = Means + SDs, ymin = Means - SDs),position = "dodge")
  ggbox = ggbox + scale_fill_manual(values=c("#FF6666","#FFFF66","#3399FF"))
  
  ggbox = ggbox + coord_cartesian(ylim=ylimits)
  # gg = gg + geom_hline(yintercept=0,color='coral',size=1)
  ggbox = ggbox + theme_bw()
  # gg = gg + scale_fill_discrete(guide = guide_legend())
  ggbox = ggbox + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                        axis.line.y = element_line(size = 0.5, colour = "black"),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        legend.text = element_text(size=25),
                        legend.title = element_blank(),
                        plot.title=element_text(size = 40),
                        axis.text.x = element_text(size = 20),
                        axis.text.y = element_text(size = 20),
                        axis.title.y = element_text(size=25),
                        axis.title.x = element_text(size=20))
}

ConstructMultipleGG = function(instancename,lossname,ylimits,box_ylimit,mean_median){
  for (mismode in c(0,1,2)){
    mygg = ConstructDrawGG(instancename,mismode,lossname,ylimits,mean_median)  
    assign(paste("gg",mismode,sep=""),mygg)
  }
  mybox = ConstructBoxGG(instancename,lossname,box_ylimit)
  mygg = ggarrange(gg0, gg1, gg2, mybox + rremove("x.text"), 
                   labels = c("Mis0", "Mis1", "Mis2","Box"),
                   ncol = 4, nrow = 1)
  return(mygg)
}
# instancename = 'Result/napkin-mismode-0-0815-0100'; ylimits = c(0.0,0.15)
# instancename = 'Result/napkin-mismode-1-0815-1100'

# instancename = 'Result/planid-mismode-0-0815-2000'; ylimits = c(0.0,0.15)
# instancename = 'Result/planid-mismode-1-0815-1100'

# instancename = 'Result/napkin-mismode-0-0816-0200'; ylimits = c(0.0,0.15)
# instancename = 'Result/napkin-mismode-1-0816-0200'; ylimits = c(0.0,0.15)
# instancename = 'Result/napkin-mismode-2-0816-0200'; ylimits = c(0.0,0.15)

# instancename = 'Result/planid-mismode-0-0816-0200'; ylimits = c(0.0,0.1)
# instancename = 'Result/planid-mismode-1-0816-0200'; ylimits = c(0.0,0.3)
# instancename = 'Result/planid-mismode-2-0816-0200'; ylimits = c(0.0,0.1)

# instancename = "planid-0816-1900"
# mismode = 0
# lossname = "weight"
# instancename = paste(instancename,"mis",mismode,sep="")
# instancename = paste(instancename,'-summary',"_",lossname,sep="")
# paste(instancename,".csv",sep="")

# instancename = "Result/napkin-0816-2300"; mismode = 0; lossname = "abs"; ylimits=c(0,0.25); boxlimits = c(0,0.5)
# instancename = "Result/planid-0816-2300"; lossname = "abs"; ylimits=c(0,0.25); boxlimits = c(0,0.5)
instancename = "Result/napkin-1817tmp"; lossname = "abs"; ylimits=c(0,0.25); boxlimits = c(0,0.5)
# df.result = ConstructDFPlot(instancename,mismode,lossname,'mean')
# gg = ConstructDrawGG(instancename,mismode,lossname,ylimits,'mean')
gg = ConstructMultipleGG(instancename,lossname,ylimits,boxlimits,'mean')
# box.result = ConstructBoxPlot(instancename,lossname)



# 
# # General 
# regmethod = 'auto'
# # xlimits = c(0,5000)
# xlimits = c(0,max(df.result$Nlist))
# spanval = 1
# point_size = 3
# alpha_point = 1
# 
# twoD = T
# medianTF = T
# 
# gg = ggplot(data = df.result, aes(x=Nlist))
# 
# 
# if(medianTF == T){
#   gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=DR.50,colour="DR.50"),size=1,method=regmethod,se=F, span=spanval)  
#   gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1.5,method=regmethod,se=F, span=spanval,linetype='dashed')
#   # gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.heuristic.50,colour="Y.heuristic.50"),size=2,method=regmethod,se=F, span=spanval,linetype='dotdash')
#   gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=WERM.50,colour="WERM.50"),size=3,method=regmethod,se=F, span=spanval,linetype='dotted')
#   gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=DR.25,ymax=DR.75),alpha=0.2, fill="dodgerblue")
#   # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.1, fill="blue")
#   gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=WERM.25,ymax=WERM.75),alpha=0.2, fill="orange")
#   gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=PlugIn.25,ymax=PlugIn.75),alpha=0.2, fill="firebrick1")
#   # gg = gg + geom_smooth(data=df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1.5,method='gam',se=F,formula=y~s(x,k=10),)
#   gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=point_size,alpha=alpha_point,shape=4)
#   gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=point_size,alpha=alpha_point,shape=16)
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.heuristic.50,colour='Y.heuristic.50'),size=point_size,alpha=alpha_point,shape=8)
#   gg = gg + geom_point(data=df.result,aes(x=Nlist,y=WERM.50,colour='WERM.50'),size=point_size,alpha=alpha_point,shape=9)
#   if(twoD == F){
#     gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.scale.50,colour="Y.scale.50"),size=1.5,method=regmethod,se=F,span=spanval)
#     gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","Y.scale.50"),values = c("gold", "red","seagreen"),labels=c("CWO","Naive","Weight-HD"))
#     # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.scale.50,colour='Y.scale.50'),size=1.5,alpha=0.2)
#   }else{
#     gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","WERM.50"),values = c("blue", "firebrick2","orange"),labels=c("DR-R","Plug-In","WERM"))
#     # gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","Y.heuristic.50","WERM.50"),values = c("blue", "firebrick2","orange","darkgreen"),labels=c("WERM-ID-R-Global","Plug-in","WERM-ID-R-Heuristic","WERM-ID"))   
#   }
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=1.5,alpha=0.2)
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=1.5,alpha=0.2)
# }else{
#   gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.multi.mean,colour="Y.multi.mean"),size=1.5,method=regmethod,se=F, span=spanval)
#   gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.naive.mean,colour="Y.naive.mean"),size=1.5,method=regmethod,se=F, span=spanval)
#   gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.multi.25,ymax=Y.multi.75),alpha=, fill="blue")
#   gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.05, fill="red")
#   # gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="Y.naive.mean"),size=1.5,method='lm',se=F,formula=y~splines::bs(x,5),span=spanval)
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.multi.mean,colour='Y.multi.mean'),size=1.5,alpha=0.2)
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.naive.mean,colour='Y.naive.mean'),size=1.5,alpha=0.2)
#   if(twoD == F){
#     gg = gg + scale_color_manual("",breaks=c("Y.multi.mean","Y.naive.mean","Y.scale.mean"),values = c("blue", "red","green"),labels=c("CWO","Naive","Weight-HD"))   
#     gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.scale.mean,colour="Y.scale.mean"),size=1.5,method=regmethod,se=F,span=2)
#     # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.scale.mean,colour='Y.scale.mean'),size=1.5,alpha=0.2)
#   }else{
#     gg = gg + scale_color_manual("",breaks=c("Y.multi.mean","Y.naive.mean"),values = c("blue", "red"),labels=c("CWO","Naive"))   
#   }
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.multi.mean,colour='Y.multi.mean'),size=1.5,alpha=0.2)
#   # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.naive.mean,colour='Y.naive.mean'),size=1.5,alpha=0.2)
# }
# 
# # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.multi.25,ymax=Y.multi.75),alpha=0.25, fill="blue")
# # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.25, fill="red")
# # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.truth.x0.25,ymax=Y.truth.x0.75),alpha=0.25, fill="green")
# 
# # gg = gg + ggtitle(Dtitle) 
# gg = gg + coord_cartesian(ylim=ylimits)
# gg = gg + theme_bw()
# gg = gg + scale_x_continuous(name = "m", limits=xlimits) + scale_y_continuous(name = "MAAE")
# gg = gg + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                 axis.line.y = element_line(size = 0.5, colour = "black"),
#                 axis.line = element_line(size=1, colour = "black"),
#                 panel.border = element_blank(),
#                 panel.background = element_blank(),
#                 legend.text = element_text(size=15),
#                 plot.title=element_text(size = 40),
#                 axis.text.x = element_text(size = 20),
#                 axis.text.y = element_text(size = 20),
#                 axis.title.y = element_text(size=25),
#                 axis.title.x = element_text(size=30)
# )



# ggbox = ggplot(box.result,aes(x=Type, y=AAE))
# ggbox = ggbox + geom_boxplot(aes(fill=label))
# # gg = gg + stat_summary(fun=mean, geom="point", aes(group=label), position=position_dodge(.9), color="red", size=3)
# # gg = gg + geom_errorbar(aes(x=Type,ymax = Means + SDs, ymin = Means - SDs),position = "dodge")
# ggbox = ggbox + scale_fill_manual(values=c("#FF6666","#FFFF66","#3399FF"))
# 
# ggbox = ggbox + coord_cartesian(ylim=ylimits)
# # gg = gg + geom_hline(yintercept=0,color='coral',size=1)
# ggbox = ggbox + theme_bw()
# # gg = gg + scale_fill_discrete(guide = guide_legend())
# ggbox = ggbox + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                 axis.line.y = element_line(size = 0.5, colour = "black"),
#                 axis.line = element_line(size=1, colour = "black"),
#                 panel.border = element_blank(),
#                 panel.background = element_blank(),
#                 legend.text = element_text(size=25),
#                 legend.title = element_blank(),
#                 plot.title=element_text(size = 40),
#                 axis.text.x = element_text(size = 20),
#                 axis.text.y = element_text(size = 20),
#                 axis.title.y = element_text(size=25),
#                 axis.title.x = element_text(size=20))


