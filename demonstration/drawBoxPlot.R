library(reshape)
library(ggplot2)
mise()
ylimits = c(0,0.4)
# probleminstance = "verma"
# probleminstance = "genfd"
# probleminstance = "fd"
probleminstance = 'napkin'
# probleminstance = 'planid'

# simResult = read.csv("verma-0112-2002-99-500-5-3.csv")
# simResult = read.csv("FD-0111-0155-15-500-5-2.csv")
# simResult = read.csv("FD-0116-0050-15-500-5-2.csv")
# simResult = read.csv("genFD-0116-0050-12-500-5-2.csv")
# simResult = read.csv("genFD-0110-2243-15-500-5-2.csv")
# simResult = read.csv('tmp_FD_0126_1033.csv')
# simResult = read.csv('tmp_verma_0126_1307.csv')
# simResult = read.csv('tmp_verma_0126_1307_035.csv')
# simResult = read.csv('tmp_genFD_0202_1109.csv')
# simResult = read.csv('tmp_FD_0202_1352.csv')
# simResult = read.csv('tmp_genFD_0203_0214.csv')
# simResult = read.csv('tmp_verma_0203_0216.csv')
# simResult = read.csv('tmp_genFD_0203_0219.csv')

# simResult = read.csv('tmp_newverma_0205_500_1211_010_030.csv')

# simResult = read.csv("tmp_newverma_0207_100_0522_020_035")

# simResult = read.csv("tmp_newverma_0207_500_0522_020_035.csv")

# simResult = read.csv("verma_0207_0544_500_020_035.csv")

# simResult = read.csv('fd1d-0216-1531-500-020-030.csv')

# simResult = read.csv('tmp_napkin_0323_1615_100_010_030.csv')
# simResult = read.csv('planid-0705-1050-100-010-030.csv')
# simResult = read.csv('napkin-0714-1626-100-03-03.csv')
# simResult = read.csv('tmp_planid.csv')
# simResult = read.csv('napkin-0718-2100-100-025-025.csv')
# simResult = read.csv('napkin-0718-2100-100-010-030.csv')
simResult = read.csv('napkin-0720-1200-100-025-025.csv')


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
