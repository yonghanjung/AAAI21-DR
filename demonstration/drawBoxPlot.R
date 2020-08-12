library(reshape)
library(ggplot2)
mise()
ylimits = c(0,0.4)
probleminstance = 'napkin'



simResult$X = c(1:nrow(simResult))
colnames(simResult)[1] = "simIdx"
myN = length(simResult$PI.Mis0)

if (probleminstance == 'fd' || probleminstance == 'genfd' || probleminstance == 'napkin' || probleminstance == 'planid'){
  AAE = c(simResult$PI.Mis0,simResult$PI.Mis1,simResult$PI.Mis2,
         simResult$WERM.Mis0,simResult$WERM.Mis1,simResult$WERM.Mis2,
       simResult$DR.Mis0,simResult$DR.Mis1,simResult$DR.Mis2)
  label = c(rep("PI",myN),rep("PI",myN),rep("PI",myN),
            rep("WERM",myN),rep("WERM",myN),rep("WERM",myN),
            rep("DR",myN),rep("DR",myN),rep("DR",myN))
  Type = c(rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN))
  plotData = data.frame(AAE,label,Type)
  plotData$label = factor(plotData$label,c("PI","WERM","DR"))
  plotData$mycolor = c(rep("red",3*myN),rep("green",3*myN),rep("blue",3*myN))
}else{
  AAE = c(simResult$param.Mis0,simResult$param.Mis1,simResult$param.Mis2,
    simResult$dl.Mis0,simResult$dl.Mis1,simResult$dl.Mis2)
  label = c(rep("PI",myN),rep("PI",myN),rep("PI",myN),
            rep("DR",myN),rep("DR",myN),rep("DR",myN))
  Type = c(rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
           rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN))
  plotData = data.frame(AAE,label,Type)
  plotData$label = factor(plotData$label,c("PI","DR"))
}

# plotData$AE = factor(plotData$AE)
# plotData$label = factor(plotData$label)
# plotData$Type = factor(plotData$Type)

gg = ggplot(plotData,aes(x=Type,y=AAE,fill=label))
gg = gg + geom_boxplot()
if (probleminstance == 'fd' || probleminstance == 'genfd' || probleminstance == 'napkin' || probleminstance == 'planid'){
  gg = gg + scale_fill_manual(values=c("#FF6666","#FFFF66","#3399FF"))
}else{
  gg = gg + scale_fill_manual(values=c("#FF6666","#3399FF"))
}
gg = gg + coord_cartesian(ylim=ylimits)
# gg = gg + geom_hline(yintercept=0,color='coral',size=1)
gg = gg + theme_bw()
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
