

library(ggplot2)
library(dplyr)
library(RColorBrewer)

TrianglePlot = function(TadRt,ChrMax,OutPlotName)
{
  TadRt = TadRt[-1,]
  # Triangle Plot use TADlevel To High.Hign is inner and low is outer.
  TadRt = mutate(TadRt,
                 high = TadRt$TADlevel)
  # generate plot data
  TriaPlotData = NULL
  for (i in 1:nrow(TadRt)) {
    lin = TadRt[i,]
    TriaPlotData = rbind(TriaPlotData,
                         c(lin$start,0,lin$midd,lin$high,lin$TADlevel),
                         c(lin$midd,lin$high,lin$end,0,lin$TADlevel)
    )
  }
  TriaPlotData = as.data.frame(TriaPlotData)
  colnames(TriaPlotData) = c("xs","ys","xe","ye","type")
  TriaPlotData$type = as.factor(TriaPlotData$type)
  # plot use segment!
  TriPlot = ggplot(data = TriaPlotData,
                   aes(x = xs,y = ys,xend = xe,yend = ye,
                       color = type)) + 
    geom_segment(size = 1) + 
    xlab("Chr Length(*50k)") + xlim(0,ChrMax)
  # theme
  TriPlot + theme_bw() + guides(color = "none") + theme(panel.grid = element_blank(),
                                                        panel.border = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.title.y = element_blank())
  ggsave(OutPlotName,width = 40,height = 5)
}

RectPlot = function(TadRt,ChrMax,OutPlotName)
{
  # RectPlot
  RectPlotData = NULL
  for (i in 1:nrow(TadRt)) {
    lin = TadRt[i,]
    RectPlotData = rbind(RectPlotData,
                         # 1+lin$TADlevel*5,100-lin$TADlevel*5
                         c(lin$start,lin$end,0,100,lin$TADlevel)
    )
  }
  RectPlotData = as.data.frame(RectPlotData)
  colnames(RectPlotData) = c("xs",'xe','ys','ye','type')
  RectPlotData$type = as.factor(RectPlotData$type) 
  
  RectPlot = ggplot(data = RectPlotData,
                    aes(xmin = xs,xmax = xe,
                        ymin = ys,ymax = ye,
                        fill = type)) + xlim(0,ChrMax) +
    xlab("TAD[Chr Length(*50k)]") + 
    geom_rect(color = 'black',size = 0.5) + 
    # coord_polar() +
    scale_fill_brewer(palette = "Greens")
  
  RectPlot + theme_bw() + guides(fill = "none") + theme(panel.grid = element_blank(),
                                                        panel.border = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.title.y = element_blank())
  ggsave(OutPlotName,width = 40,height = 5)
}

# Read Tads
TadRt = read.table(snakemake@input[[1]],sep = "\t",header = F)
colnames(TadRt) = c("start","end","TADlevel","mean",'score')
ChrMax = max(TadRt$end)

# Basic Data Dispose
TadRt = mutate(TadRt,midd = (start+end)/2)
TadRt = arrange(TadRt,TADlevel)
# Filter trust tads
TadRt = filter(TadRt,score >= max(TadRt$score[-1])*0.1)

# Plots
TrianglePlot(TadRt,ChrMax,snakemake@output[[1]])
RectPlot(TadRt,ChrMax,snakemake@output[[2]])

