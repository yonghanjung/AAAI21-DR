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
  
  confidence_coef = 1/5
  
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
    
    # df.result.WERM = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='WERM')))
    # vector.WERM = df.result.WERM[,ncol(df.result.WERM)]
    
    assign(paste('vector.DR.mis',mismode,sep=""),vector.DR)
    assign(paste('vector.PlugIn.mis',mismode,sep=""),vector.PlugIn)
    # assign(paste('vector.WERM.mis',mismode,sep=""),vector.WERM)
  }
  
  myN = length(vector.PlugIn.mis0)
  WAAE = c(vector.PlugIn.mis0,vector.PlugIn.mis1,vector.PlugIn.mis2,
          # vector.WERM.mis0,vector.WERM.mis1,vector.WERM.mis2,
          vector.DR.mis0,vector.DR.mis1,vector.DR.mis2)
  label = c(rep("PI",myN),rep("PI",myN),rep("PI",myN),
            # rep("WERM",myN),rep("WERM",myN),rep("WERM",myN),
            rep("DML",myN),rep("DML",myN),rep("DML",myN))
  Type = c(rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           # rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN))
  plotData = data.frame(WAAE,label,Type)
  # plotData$label = factor(plotData$label,c("PI","WERM","DR"))
  plotData$label = factor(plotData$label,c("PI","DML"))
  plotData$mycolor = c(rep("red",3*myN),rep("blue",3*myN))
  return(plotData)
}



ConstructDrawGG = function(instancename,mismode,lossname,ylimits,xlimits,mean_median){
  df.result = ConstructDFPlot(instancename,mismode,lossname,mean_median)
  
  # General 
  regmethod = 'auto'
  # xlimits = c(0,5000)
  # xlimits = c(0,max(df.result$Nlist))
  spanval = 0.8
  point_size = 3
  alpha_point = 1
  
  twoD = T
  medianTF = T
  
  gg = ggplot(data = df.result, aes(x=Nlist))
  
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=DR.50,colour="DR.50"),size=1,method=regmethod,se=F, span=spanval)  
  gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=PlugIn.50,colour="PlugIn.50"),size=1,method=regmethod,se=F, span=spanval)
  # gg = gg + geom_smooth(data = df.result, aes(x=Nlist,y=WERM.50,colour="WERM.50"),size=1,method=regmethod,se=F, span=spanval)
    
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=DR.25,ymax=DR.75),alpha=0.2, fill="dodgerblue")
  gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=PlugIn.25,ymax=PlugIn.75),alpha=0.2, fill="firebrick1")
  # gg = gg + geom_ribbon(data=df.result, aes(x=Nlist, ymin=WERM.25,ymax=WERM.75),alpha=0.2, fill="orange")

  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=PlugIn.50,colour='PlugIn.50'),size=point_size,alpha=alpha_point,shape=4)
  gg = gg + geom_point(data=df.result,aes(x=Nlist,y=DR.50,colour='DR.50'),size=point_size,alpha=alpha_point,shape=16)
  # gg = gg + geom_point(data=df.result,aes(x=Nlist,y=WERM.50,colour='WERM.50'),size=point_size,alpha=alpha_point,shape=9)

  # gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50","WERM.50"),values = c("blue", "firebrick2","orange"),labels=c("DR","Plug-In","WERM"))
  gg = gg + scale_color_manual("",breaks=c("DR.50","PlugIn.50"),values = c("blue", "firebrick2"),labels=c("DML","Plug-In"))
  
  # gg = gg + ggtitle(Dtitle) 
  gg = gg + coord_cartesian(ylim=ylimits,xlim=xlimits)
  gg = gg + theme_bw()
  gg = gg + scale_x_continuous(name = "N", limits=xlimits) + scale_y_continuous(name = "WAAE")
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
  ggbox = ggplot(box.result,aes(x=Type, y=WAAE))
  ggbox = ggbox + geom_boxplot(aes(fill=label))
  # gg = gg + stat_summary(fun=mean, geom="point", aes(group=label), position=position_dodge(.9), color="red", size=3)
  # gg = gg + geom_errorbar(aes(x=Type,ymax = Means + SDs, ymin = Means - SDs),position = "dodge")
  # ggbox = ggbox + scale_fill_manual(values=c("#FF6666","#FFFF66","#3399FF"))
  ggbox = ggbox + scale_fill_manual(values=c("#FF6666","#FFFF66"))
  # ggbox = ggbox + scale_fill_manual(values=c("blue","firebrick2","orange"))
  
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

ConstructErrorGG = function(instancename,lossname,ylimits){
  # Box plot with mean and sd 
  box.result = ConstructBoxPlot(instancename,lossname)
  mySummary = c()
  for (mylabel in c('DML','PI')){
    for (myType in c("Correct","Mis-1","Mis-2")){
      filtered_box = subset(box.result,label==mylabel&Type==myType)    
      
      myMean = mean(filtered_box[,'WAAE'])
      mySD = sd(filtered_box[,'WAAE'])
      myLow = max(myMean-mySD,0)
      myHigh = min(myMean+mySD,1)
      
      mycolor = unique(filtered_box[,'mycolor'])
      mySummary = rbind(mySummary,c(myMean,mySD,mylabel,myType,mycolor,myLow,myHigh))
    }
  }
  # for (mylabel in c('DR','PI','WERM')){
  #   for (myType in c("Correct","Mis-1","Mis-2")){
  #     filtered_box = subset(box.result,label==mylabel&Type==myType)    
  #     
  #     myMean = mean(filtered_box[,'WAAE'])
  #     mySD = sd(filtered_box[,'WAAE'])
  #     myLow = max(myMean-mySD,0)
  #     myHigh = min(myMean+mySD,1)
  #     
  #     mycolor = unique(filtered_box[,'mycolor'])
  #     mySummary = rbind(mySummary,c(myMean,mySD,mylabel,myType,mycolor,myLow,myHigh))
  #   }
  # }
  colnames(mySummary) = c('WAAE','SD','Label','Scenario','Color','Low','High')
  mySummary = as.data.frame(mySummary)
  mySummary[,'WAAE'] = as.numeric(mySummary[,'WAAE'])
  mySummary[,'SD'] = as.numeric(mySummary[,'SD'])
  mySummary[,'Low'] = as.numeric(mySummary[,'Low'])
  mySummary[,'High'] = as.numeric(mySummary[,'High'])
  # mySummary[,'order'] = c(1,2,3,7,8,9,4,5,6)
  mySummary[,'order'] = c(1,2,3,4,5,6)
  mySummary = mySummary[order(mySummary[,'order']),]
  # mySummary[,'group'] = c(rep('A',3),rep('B',3),rep('C',3))
  mySummary[,'group'] = c(rep('A',3),rep('B',3))
  
  # boxlimits = c(0,0.3)
  myPosition = position_dodge(.5)
  gg = ggplot(mySummary, aes(x=Scenario,y=WAAE,color=group,ymin=Low,ymax=High))
  gg = gg + geom_point(shape=15,size=5,position=myPosition)
  gg = gg + geom_errorbar(position = myPosition,width=0.3,size=1.5)
  # gg = gg + scale_color_manual("Label",values=c("#3399FF","orange","firebrick1"))
  gg = gg + scale_color_manual("Label",values=c("#3399FF","firebrick1"))
  gg = gg + coord_cartesian(ylim=boxlimits)
  gg = gg = gg + theme_bw()
  gg = gg + guides(color=guide_legend())
  # gg = gg + scale_fill_discrete(guide = guide_legend())
  gg = gg + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
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
  # c("#FF6666","#FFFF66","#3399FF")
  return(gg)
}

ConstructMultipleGG = function(instancename,lossname,ylimits0,ylimits1,ylimits2,box_ylimit,xlimits,mean_median){
  for (mismode in c(0,1,2)){
    myylimit = get(paste('ylimits',mismode,sep=""))
    mygg = ConstructDrawGG(instancename,mismode,lossname,myylimit,xlimits,mean_median)  
    assign(paste("gg",mismode,sep=""),mygg)
  }
  # mybox = ConstructBoxGG(instancename,lossname,box_ylimit)
  mybox = ConstructErrorGG(instancename,lossname,box_ylimit)
  mygg = ggarrange(gg0, gg1, gg2, mybox + rremove("x.text"), 
                   labels = c("Mis0", "Mis1", "Mis2","Box"),
                   ncol = 4, nrow = 1)
  # mygg = mygg + ggtitle("HAO")
  return(mygg)
}

ConstructGGOutput = function(instancename,lossname,ylimits0,ylimits1,ylimits2,box_ylimit,xlimits,mean_median){
  for (mismode in c(0,1,2)){
    myylimit = get(paste('ylimits',mismode,sep=""))
    mygg = ConstructDrawGG(instancename,mismode,lossname,myylimit,xlimits,mean_median)  
    assign(paste("gg",mismode,sep=""),mygg)
  }
  # mybox = ConstructBoxGG(instancename,lossname,box_ylimit)
  mybox = ConstructErrorGG(instancename,lossname,box_ylimit)
  myresult = list(gg0,gg1,gg2,mybox)
  return(myresult)
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
# instancename = "Result/napkin-1817tmp"; lossname = "abs"; ylimits=c(0,0.25); boxlimits = c(0,0.5)

# instancename = "Result/planid-0817-1830"; lossname = "abs"; ylimits0=c(0,0.1); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); boxlimits = c(0,0.5)
# instancename = "Result/napkin-0817-2200"; lossname = "abs"; ylimits0=c(0,0.1); ylimits1 = c(0,0.3); ylimits2 = c(0,0.1); boxlimits = c(0,0.5)

# instancename = "Result/napkin-0818-2200-tmp"; lossname = "weight"; ylimits0=c(0,0.1); ylimits1 = c(0,0.25); ylimits2 = c(0,0.1); boxlimits = c(0,0.5)
# instancename = "Result/planid-0818-2200-tmp"; lossname = "weight"; ylimits0=c(0,0.1); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); boxlimits = c(0,0.5)




# Reweight 
# instancename = "Result/planid-0819-1330"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); xlimits = c(0,8000); boxlimits = c(0,0.5)

# Noreweight / No regularization ==> FINAL 
# instancename = "Result/planid-0821-0130-noreweight"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); xlimits = c(0,5000); boxlimits = c(0,0.3)

# Noreweight-ver2 / No regularization 
# instancename = "Result/planid-0821-0130-noreweight-ver2"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); xlimits = c(0,5000); boxlimits = c(0,0.3)

# Final 
# instancename = "Result/napkin-0819-0230"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.25); ylimits2 = c(0,0.15); xlimits = c(0,10000); boxlimits = c(0,0.25) # Final
# instancename = "Result/planid-0821-0130-noreweight"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.3); ylimits2 = c(0,0.15); xlimits = c(0,5000); boxlimits = c(0,0.3)

# For DML presenting 
# instancename = "Result/napkin-0908-1630"; lossname = "weight"; ylimits0=c(0,0.15); ylimits1 = c(0,0.25); ylimits2 = c(0,0.15); xlimits = c(0,10000); boxlimits = c(0,0.25) # Final
# instancename = "Result/napkin-0908-2000"; lossname = "weight"; ylimits0=c(0,0.25); ylimits1 = c(0,0.25); ylimits2 = c(0,0.25); xlimits = c(0,10000); boxlimits = c(0,0.25) # Final
# instancename = "Result/planid-0908-2000"; lossname = "abs"; ylimits0=c(0,0.25); ylimits1 = c(0,0.25); ylimits2 = c(0,0.25); xlimits = c(0,10000); boxlimits = c(0,0.25) # Final

# DML Final 
instancename = "Result/napkin-0909-0000"; lossname = "weight"; ylimits0=c(0,0.25); ylimits1 = c(0,0.25); ylimits2 = c(0,0.25); xlimits = c(0,10000); boxlimits = c(0,0.25); mean_median = 'mean' # Final


mismode = 2
df.result.summary = read.csv(MakeFileName(instancename,mismode,lossname,mymodel='summary'))
df.result.DR = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='DR')))
df.result.PlugIn = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='PlugIn')))
df.result.WERM = ReadCsv(read.csv(MakeFileName(instancename,mismode,lossname,mymodel='WERM')))

colLength = 20
varDR = rep(0,colLength)
varPI = rep(0,colLength)
varWERM = rep(0,colLength)
for (colidx in 1:colLength){
  varDR[colidx] = sd(df.result.DR[,colidx])
  varPI[colidx] = sd(df.result.PlugIn[,colidx])
  varWERM[colidx] = sd(df.result.WERM[,colidx])  
}
myVar = data.frame(X=c(1:colLength),varDR,varPI,varWERM)
varGG = ggplot() + geom_smooth(data=myVar, aes(x=X,y=varDR),color="blue") + geom_point(data=myVar, aes(x=X,y=varDR),color="blue")
varGG = varGG + geom_smooth(data=myVar, aes(x=X,y=varPI),color="red") + geom_point(data=myVar, aes(x=X,y=varPI),color="red")
varGG = varGG + geom_smooth(data=myVar, aes(x=X,y=varWERM),color="orange") + geom_point(data=myVar, aes(x=X,y=varWERM),color="orange")



gg = ConstructMultipleGG(instancename,lossname,ylimits0,ylimits1,ylimits2,boxlimits,xlimits,mean_median)
myresult = ConstructGGOutput(instancename,lossname,ylimits0,ylimits1,ylimits2,boxlimits,xlimits,mean_median)

# File Save
FileSave = F 
if (FileSave){
  print("FileSave")
  mismode = 2
  # filename = paste("Result/Plot/napkin/napkin-mis",mismode,".pdf",sep="")
  filename = paste("Result/Plot/planid/planid-mis",mismode,".pdf",sep="")
  # filename = paste("Result/Plot/planid/planid-mis",mismode,".pdf",sep="")
  pdf(file=filename, width=8,heigh=6)
  myresult[[(mismode+1)]]  
  dev.off()
  
  # filename = "Result/Plot/napkin/napkin-box.pdf"
  filename = "Result/Plot/planid/planid-box.pdf"
  pdf(file=filename, width=8,heigh=6)
  myresult[[4]]
  dev.off()
}








