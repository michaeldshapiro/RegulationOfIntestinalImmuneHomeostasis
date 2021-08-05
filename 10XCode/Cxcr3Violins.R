
library(HandyPack)
library(Seurat)
library(ggplot2)

rm(list=ls())
graphics.off()

source('SeuratWare02.R')

## ###################################################
saveFigureDF = function(df,fileName)
{
    fileName = str_replace(fileName,'violinPlots/','figureTables/violinPlot_')
    fileName = str_replace(fileName,'pdf','txt')

    Write.Table(df,fileName)
}


## ###################################################
f = getSeuratObject()
f@active.assay = 'RNA'

cells = getShortNames()[16:30]
print(cells %in% getShortNames())


gene = 'Cxcr3'
print(gene %in% rownames(f))

idx = f@meta.data$shortName %in% cells
fPrime = f[gene,idx]
expression = fPrime@assays$RNA@data
expression = as.matrix(expression)
expression = as.numeric(expression)

df = data.frame(cell=fPrime@meta.data$shortName,expression)
df$cell = factor(df$cell,levels=getShortNames())

g = ggplot(df,aes(x=cell,y=expression,fill=cell)) +
    geom_violin(draw_quantiles= c(0.25, 0.5, 0.75),scale='width') +
    scale_fill_manual(values=getThePalette(),
                      breaks=getShortNames()) +
    geom_jitter(width=0.2) +
    theme(axis.text.x = element_text(angle = 80,hjust=1,size=10)) +
    theme(legend.position='none') +
    ggtitle(gene)

print(g)

ggsave(plot=g,
       filename='violinPlots/Cxcr3.pdf',
       height=8,width=15,units='in')

saveFigureDF(df,'violinPlots/Cxcr3.pdf')
