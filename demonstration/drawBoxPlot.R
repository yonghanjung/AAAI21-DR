library(mise)
library(reshape)
library(ggplot2)
mise()
ylimits = c(0,0.5)

df.result = read.csv('Result/napkin-box-0812-0230.csv')
df.result$X = c(1:nrow(df.result))
colnames(df.result)[1] = "simIdx"
myN = nrow(df.result)

AAE = c(df.result$PI.Mis0,df.result$PI.Mis1,df.result$PI.Mis2,
        df.result$WERM.Mis0,df.result$WERM.Mis1,df.result$WERM.Mis2,
        df.result$DR.Mis0,df.result$DR.Mis1,df.result$DR.Mis2)
Means = c(mean(df.result$PI.Mis0), mean(df.result$PI.Mis1), mean(df.result$PI.Mis2),
          mean(df.result$WERM.Mis0), mean(df.result$WERM.Mis1), mean(df.result$WERM.Mis2),
          mean(df.result$DR.Mis0), mean(df.result$DR.Mis1), mean(df.result$DR.Mis2))
SDs = c(sd(df.result$PI.Mis0), sd(df.result$PI.Mis1), sd(df.result$PI.Mis2),
        sd(df.result$WERM.Mis0), sd(df.result$WERM.Mis1), sd(df.result$WERM.Mis2),
        sd(df.result$DR.Mis0), sd(df.result$DR.Mis1), sd(df.result$DR.Mis2))
label = c(rep("PI",myN),rep("PI",myN),rep("PI",myN),
          rep("WERM",myN),rep("WERM",myN),rep("WERM",myN),
          rep("DR",myN),rep("DR",myN),rep("DR",myN))
Type = c(rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
         rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN),
         rep("Correct",myN),rep("Mis-1",myN),rep("Mis-2",myN))
plotData = data.frame(AAE,Means,SDs,label,Type)
plotData$label = factor(plotData$label,c("PI","WERM","DR"))
plotData$mycolor = c(rep("red",3*myN),rep("green",3*myN),rep("blue",3*myN))



# plotData$AE = factor(plotData$AE)
# plotData$label = factor(plotData$label)
# plotData$Type = factor(plotData$Type)

gg = ggplot(plotData,aes(x=Type, y=AAE))
gg = gg + geom_boxplot(aes(fill=label))
# gg = gg + stat_summary(fun=mean, geom="point", aes(group=label), position=position_dodge(.9), color="red", size=3)
# gg = gg + geom_errorbar(aes(x=Type,ymax = Means + SDs, ymin = Means - SDs),position = "dodge")
gg = gg + scale_fill_manual(values=c("#FF6666","#FFFF66","#3399FF"))

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
