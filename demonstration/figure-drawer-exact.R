library(ggplot2)
library(cowplot)
library(mise)
library(mgcv)

# mise()
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

ConstructDFPlot = function(instancename,mean_median){
  df.result.summary = read.csv(paste(instancename,'-summary.csv',sep=""))
  df.result.DR = read.csv(paste(instancename,'-DR.csv',sep=""))
  df.result.PlugIn = read.csv(paste(instancename,'-PlugIn.csv',sep=""))
  df.result.WERM = read.csv(paste(instancename,'-WERM.csv',sep=""))
  
  confidence_coef = 1/2
  
  if (mean_median == "median"){
    DRCenter = df.result.summary$DR.50
    PlugInCenter = df.result.summary$PlugIn.50
    WERMCenter = df.result.summary$WERM.50
  }else{
    DRCenter = df.result.summary$DR.mean
    PlugInCenter = df.result.summary$PlugIn.mean
    WERMCenter = computeColMeans(df.result.PlugIn)
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


instancename = 'Result/napkin-mismode-0-0810-1300'
# instancename = 'Result/napkin-mismode-1-0811-1800'
# instancename = 'Result/napkin-mismode-2-0811-2030'
df.result = ConstructDFPlot(instancename,'mean')


# df.result.time = ConstructTimePlot(instancename)
# df.result.time = read.csv('Result/napkin-0802-1030-D20-time-global.csv')
# df.result.time = t(df.result.time)
# df.result.time = df.result.time[c(2:nrow(df.result.time)),]

# General 
regmethod = 'auto'
ylimits = c(0.0,0.1)
# xlimits = c(0,5000)
xlimits = c(0,max(df.result$Nlist))
spanval = 1
point_size = 3
alpha_point = 1

twoD = T
medianTF = T



gg = ggplot(data = df.result, aes(x=Nlist))


if(medianTF == T){
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=DR.50,colour="DR.50"),size=1,method=regmethod,se=F, span=spanval)  
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1.5,method=regmethod,se=F, span=spanval,linetype='dashed')
  # gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.heuristic.50,colour="Y.heuristic.50"),size=2,method=regmethod,se=F, span=spanval,linetype='dotdash')
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=WERM.50,colour="WERM.50"),size=3,method=regmethod,se=F, span=spanval,linetype='dotted')
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=DR.25,ymax=DR.75),alpha=0.2, fill="dodgerblue")
  # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.1, fill="blue")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=WERM.25,ymax=WERM.75),alpha=0.2, fill="orange")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=PlugIn.25,ymax=PlugIn.75),alpha=0.2, fill="firebrick1")
  # gg = gg + geom_smooth(data=df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1.5,method='gam',se=F,formula=y~s(x,k=10),)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=point_size,alpha=alpha_point,shape=4)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=point_size,alpha=alpha_point,shape=16)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.heuristic.50,colour='Y.heuristic.50'),size=point_size,alpha=alpha_point,shape=8)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=WERM.50,colour='WERM.50'),size=point_size,alpha=alpha_point,shape=9)
  if(twoD == F){
    gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.scale.50,colour="Y.scale.50"),size=1.5,method=regmethod,se=F,span=spanval)
    gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","Y.scale.50"),values = c("gold", "red","seagreen"),labels=c("CWO","Naive","Weight-HD"))
    # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.scale.50,colour='Y.scale.50'),size=1.5,alpha=0.2)
  }else{
    gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","WERM.50"),values = c("blue", "firebrick2","orange"),labels=c("DR-R","Plug-In","WERM"))
    # gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","Y.heuristic.50","WERM.50"),values = c("blue", "firebrick2","orange","darkgreen"),labels=c("WERM-ID-R-Global","Plug-in","WERM-ID-R-Heuristic","WERM-ID"))   
  }
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=1.5,alpha=0.2)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=1.5,alpha=0.2)
}else{
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.multi.mean,colour="Y.multi.mean"),size=1.5,method=regmethod,se=F, span=spanval)
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.naive.mean,colour="Y.naive.mean"),size=1.5,method=regmethod,se=F, span=spanval)
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.multi.25,ymax=Y.multi.75),alpha=, fill="blue")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.05, fill="red")
  # gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="Y.naive.mean"),size=1.5,method='lm',se=F,formula=y~splines::bs(x,5),span=spanval)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.multi.mean,colour='Y.multi.mean'),size=1.5,alpha=0.2)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.naive.mean,colour='Y.naive.mean'),size=1.5,alpha=0.2)
  if(twoD == F){
    gg = gg + scale_color_manual("",breaks=c("Y.multi.mean","Y.naive.mean","Y.scale.mean"),values = c("blue", "red","green"),labels=c("CWO","Naive","Weight-HD"))   
    gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=Y.scale.mean,colour="Y.scale.mean"),size=1.5,method=regmethod,se=F,span=2)
    # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.scale.mean,colour='Y.scale.mean'),size=1.5,alpha=0.2)
  }else{
    gg = gg + scale_color_manual("",breaks=c("Y.multi.mean","Y.naive.mean"),values = c("blue", "red"),labels=c("CWO","Naive"))   
  }
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.multi.mean,colour='Y.multi.mean'),size=1.5,alpha=0.2)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=Y.naive.mean,colour='Y.naive.mean'),size=1.5,alpha=0.2)
}

# gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.multi.25,ymax=Y.multi.75),alpha=0.25, fill="blue")
# gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.naive.25,ymax=Y.naive.75),alpha=0.25, fill="red")
# gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=Y.truth.x0.25,ymax=Y.truth.x0.75),alpha=0.25, fill="green")

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
gg