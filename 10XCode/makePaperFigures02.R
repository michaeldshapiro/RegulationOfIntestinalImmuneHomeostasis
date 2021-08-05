

library(HandyPack)
library(ggplot2)

rm(list=ls())
graphics.off()

source('SeuratWare02.R')


## ###################################################
## ###################################################

getFigureGenes = function()
{
    df = Read.Table('paperFigures/GenesForDotPlotFigure.csv')
    return(unique(df$gene))
}

## ###################################################
denude = function(p)
{
    pPrime = p +
        theme(legend.position = "none",
              axis.text=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              title=element_blank())

    return(pPrime)
}



## ###################################################
## ###################################################

genes = getFigureGenes()
f = getSeuratObject()
f@active.assay = 'RNA'

groups = getGroups()
groupNames = unlist(lapply(groups,getShortName))

N = length(genes) * length(groups)
df = data.frame(group=character(N),
                gene=character(N),
                expressedIn=numeric(N),
                meanExpression=numeric(N),
                stringsAsFactors=FALSE)

M = matrix(0,nrow=length(groups),ncol=length(genes))
rownames(M) = groupNames
colnames(M) = genes

cutoff = 0.02



finger = 1
for(i in 1:length(groups))
{
    group = groups[i]
    ## ###################################################
    ## Subset to group:
    idx = f@meta.data$seurat_cluster == group
    fGroup = f[,idx]
    numInGroup = sum(idx)
    expression = as.matrix(fGroup[['RNA']]@data)
    
    for(gene in genes)
    {
        writeLines(paste(group,gene))
        
        df$group[finger] = getShortName(group) 
        df$gene[finger] = gene

        idx = expression[gene,] >= cutoff
        df$expressedIn[finger] = sum(idx) / numInGroup
        df$meanExpression[finger] = mean(expression[gene,idx])
        M[groupNames[i],gene] = mean(expression[gene,idx])
        finger = finger + 1
    }
}

h = heatmap(M)
df$group = factor(df$group,levels=rev(groupNames))
df$gene = factor(df$gene,levels=genes)

dfFiltered = df
cutoff = .25
idx = dfFiltered$expressedIn <= cutoff
dfFiltered = dfFiltered[!idx,]

textSize = 15

g = ggplot(dfFiltered,aes(x=gene,y=group,size=expressedIn,color=meanExpression)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45,hjust=1,size=textSize),
          axis.text.y = element_text(size=textSize)) +
    scale_color_viridis_c()

print(g)

if(! dir.exists('paperFigures'))
    dir.create('paperFigures')

ggsave(plot=g,
       filename='paperFigures/genesVCellTypesDotPlotFiltered.pdf',
       height=10,width=20,units='in')

ggsave(plot=denude(g),
       filename='paperFigures/genesVCellTypesDotPlotFilteredNoFrills.pdf',
       height=10,width=20,units='in')

Write.Table(dfFiltered,'figureTables/genesVCellTypesDotPlot.txt')

