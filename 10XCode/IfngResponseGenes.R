
library(Seurat)
library(HandyPack)
library(ggplot2)
library(ggsignif)



source('SeuratWare02.R')
## ###################################################
f = getSeuratObject()
f@active.assay = 'RNA'

responseGenes = Read.Table('IfngResponseGenes.txt')
responseGenes = responseGenes$genes
idx = responseGenes %in% rownames(f)
responseGenes = responseGenes[idx]

df = data.frame(f@reductions$umap@cell.embeddings)
df$quadrant = f@meta.data$quadrant
df$quadrant = factor(df$quadrant,levels=getQuadrants())
df$treatment = f@meta.data$treatment


M = f@assays$RNA@data
M = as.matrix(M)
M = M[responseGenes,]
expression = colMeans(M)
df$expression = expression

df = df[df$treatment=='infected',]



gViolin = ggplot(df,aes(x=quadrant,y=expression,fill=quadrant)) +
    geom_violin(draw_quantiles= c(0.25, 0.5, 0.75)) +
    scale_fill_manual(values=getQuadrantPalette(),
                      breaks=getQuadrants()) +
    ggtitle('Ifng response genes') +
    geom_signif(comparisons=list(c('WT_infected','KO_infected')),
                map_signif_level=TRUE,textsize=6)

print(gViolin)
ggsave(plot=gViolin,
       filename='violinPlots/IfngResponseGenesViolin.pdf')

Write.Table(df,
            'figureTables/violinPlot_IfngResponseGenes.txt')


idx = df$quadrant == 'WT_infected'
x = df[idx,'expression']
y = df[!idx,'expression']
test = wilcox.test(x,y)

print(test$p.value)
