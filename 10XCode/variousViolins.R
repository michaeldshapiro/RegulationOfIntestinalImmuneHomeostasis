
library(Seurat)
library(ggplot2)
library(HandyPack)
library(stringr)
library(ggsignif)

rm(list=ls())
graphics.off()

source('SeuratWare02.R')

## ###################################################
getExpressionDF = function(f,cellGroups,genes,quadrants)
{
    ## Get expression into a data frame:
    df = data.frame(quadrant=character(0),
                    cells=character(0),
                    gene=character(0),
                    expression=numeric(0))
                    

    ## Iterate over quadrants:
    for(q in quadrants)
    {
        ## Iterate through the groups by name:
        for(n in names(cellGroups))
        {
            ## Iterate over the cells in this group:
            for(cell in cellGroups[[n]])
            {
                idx = (f@meta.data$quadrant == q &
                       f@meta.data$shortName == cell)

                if(sum(idx) == 0)
                {
                    print(paste('skipping',q,cell))
                    next
                }
                
                fThis = f[,idx]
                M = fThis@assays$RNA@data
                M = as.matrix(M)

                ## Iterate over genes:
                for(gene in genes)
                {
                    expression = M[gene,]
                    a = data.frame(quadrant=q,
                                   cells=n,
                                   gene=gene,
                                   expression=expression)

                    df = rbind(df,a)
                }
            }
        }
    }

    df = exifyDF(df,genes)

    return(df)
}

## ###################################################
exifyDF = function(df,genes)
{
    df$gene = factor(df$gene,levels=genes)
    df$quadrant = factor(df$quadrant,levels=getQuadrants())

    df$x = paste(df$gene,df$cells,df$quadrant)

    levels = c()
    for(gene in genes)
        for(cell in unique(df$cells))
            for(q in unique(df$quadrant))
            {
                a = paste(gene,cell,q)
                levels = c(levels,a)
            }

    df$x = factor(df$x,levels=levels)

    return(df)
}

## ###################################################
makeComparisons = function(df)
{
    genes = unique(df$gene)
    cells = unique(df$cells)
    quadrants = unique(df$quadrant)

    shortX = c()
    for(gene in genes)
        for(cell in cells)
            shortX = c(shortX,paste(gene,cell))

    comparisons = list()
    for(i in 1:length(shortX))
        comparisons[[i]] = c(paste(shortX[i],quadrants[1]),
                             paste(shortX[i],quadrants[2]))

    return(comparisons)
}

## ###################################################
makeThePlot = function(df,title)
{
    g = ggplot(df,aes(x=x,y=expression,fill=quadrant)) +
        geom_violin(draw_quantiles= c(0.25, 0.5, 0.75),scale='width') +
        theme(axis.text.x = element_text(angle = 80,hjust=1,size=10)) +
        scale_fill_manual(values=getQuadrantPalette(),
                          breaks=getQuadrants()) +
        geom_signif(comparisons=makeComparisons(df),
                    map_signif_level=TRUE,textsize=6) +
        xlab('') +
        ggtitle(title)

    g = stretchY(g)

    return(g)
}

## ###################################################
makeIfngPlot = function(df,title)
{
    g = ggplot(df,aes(x=x,y=expression,fill=cells)) +
        geom_violin(draw_quantiles= c(0.25, 0.5, 0.75),scale='width') +
        geom_jitter(width=0.2) +
        theme(axis.text.x = element_text(angle = 80,hjust=1,size=10)) +
        scale_fill_manual(values=getThePalette(),
                          breaks=getShortNames()) +
        geom_signif(comparisons=makeComparisons(df),
                    map_signif_level=TRUE,textsize=6) +
        xlab('') +
        ggtitle(title)
    
    g = stretchY(g)
    
    return(g)
}

## ###################################################
stretchY = function(g,bump=.2)
{
    a = ggplot_build(g)
    y = a$layout$panel_scales_y[[1]]$range$range[2]
    g = g + ylim(0,y+bump)

    return(g)
}

## ###################################################
saveFigureDF = function(df,fileName)
{
    fileName = str_replace(fileName,'violinPlots/','figureTables/violinPlot_')
    fileName = str_replace(fileName,'pdf','txt')

    Write.Table(df,fileName)
}

## ###################################################
## ###################################################
f = getSeuratObject()
f@active.assay = 'RNA'

## ###################################################
## Ifng:
cellGroups = list(OtherLymphoidCells='OtherLymphoidCells',
                  NaturalKillerCells='NaturalKillerCells',
                  TCells1='TCells1',
                  TCells2='TCells2',
                  TCells3='TCells3')
genes = 'Ifng'
quadrants = c('WT_infected','KO_infected')

ifngDF = getExpressionDF(f,cellGroups,genes,quadrants)
gIfng = makeIfngPlot(ifngDF,'Ifng')

print(gIfng)
ggsave(plot=gIfng,
       filename='violinPlots/Ifng.pdf',
       height=8,width=15,units='in')

saveFigureDF(ifngDF,'violinPlots/Ifng.pdf')


## ###################################################
## Mesothelial cells, Lcn2, Saa3, Msln, Upk3b 
cellGroups = list(Mesothelial=c('MesothelialCells1',
                                'MesothelialCells2',
                                'MesothelialCells3',
                                'MesothelialCells4',
                                'MesothelialCells5'))
genes = c('Lcn2', 'Saa3', 'Msln', 'Upk3b') 
quadrants = c('WT_uninfected','KO_uninfected')

mesothelialDF = getExpressionDF(f,cellGroups,genes,quadrants)
gMesothelial = makeThePlot(mesothelialDF,'Mesothelial Cells')
dev.new()
print(gMesothelial)
ggsave(plot=gMesothelial,
       filename='violinPlots/Mesothelial.pdf',
       height=8,width=15,units='in')

saveFigureDF(mesothelialDF,'violinPlots/Mesothelial.pdf')

## ###################################################
## Il33, Cxcl1, Il6 in Fibroblasts
cellGroups = list(Fibroblasts=c('Fibroblasts1',
                                'Fibroblasts2'))
genes = c('Il33','Cxcl1','Il6')
quadrants = c('WT_uninfected','KO_uninfected')

FibroblastsDF = getExpressionDF(f,cellGroups,genes,quadrants)
gFibroblasts = makeThePlot(FibroblastsDF,'Fibroblasts')
dev.new()
print(gFibroblasts)
ggsave(plot=gFibroblasts,
       filename='violinPlots/Fibroblasts.pdf',
       height=8,width=15,units='in')

saveFigureDF(FibroblastsDF,'violinPlots/Fibroblasts.pdf')

## ###################################################
## Ccl8, Ccl7, Cxcl2, Ccl12, Ccl2, Ccl4, Il1b in Muscularis Macrophages
cellGroups = list(MuscularisMacrophages=c('MuscularisMacrophages1',
                                          'MuscularisMacrophages2',
                                          'MuscularisMacrophages3',
                                          'MuscularisMacrophages4'))
genes = c('Ccl8','Ccl7','Cxcl2','Ccl12','Ccl2','Ccl4','Il1b')
quadrants = c('WT_uninfected','KO_uninfected')

MuscularisMacrophagesDF = getExpressionDF(f,cellGroups,genes,quadrants)
gMuscularis = makeThePlot(MuscularisMacrophagesDF,'Muscularis Macrophages')
dev.new()
print(gMuscularis)
ggsave(plot=gMuscularis,
       filename='violinPlots/MuscularisMacrophages.pdf',
       height=8,width=15,units='in')

saveFigureDF(MuscularisMacrophagesDF,'violinPlots/MuscularisMacrophages.pdf')

## ###################################################
## Arg1, Retnla, Chil3 in Muscularis Macrophages
cellGroups = list(MuscularisMacrophages=c('MuscularisMacrophages1',
                                          'MuscularisMacrophages2',
                                          'MuscularisMacrophages3',
                                          'MuscularisMacrophages4'))
genes = c('Arg1','Retnla','Chil3')
quadrants = c('WT_uninfected','KO_uninfected')

MuscularisMacrophagesDFTwo = getExpressionDF(f,cellGroups,genes,quadrants)
gMuscularisTwo = makeThePlot(MuscularisMacrophagesDFTwo,'Muscularis Macrophages')
dev.new()
print(gMuscularisTwo)
ggsave(plot=gMuscularisTwo,
       filename='violinPlots/MuscularisMacrophagesTwo.pdf',
       height=8,width=15,units='in')

saveFigureDF(MuscularisMacrophagesDFTwo,'violinPlots/MuscularisMacrophagesTwo.pdf')

## ###################################################
## Stat1 in all cells
cellGroups = list(All=getShortNames())
genes = c('Stat1')
quadrants = c('WT_infected','KO_infected')

AllDF = getExpressionDF(f,cellGroups,genes,quadrants)
gAll = makeThePlot(AllDF,'All Cells')
dev.new()
print(gAll)
ggsave(plot=gAll,
       filename='violinPlots/AllCells.pdf',
       height=8,width=15,units='in')

saveFigureDF(AllDF,'violinPlots/AllCells.pdf')


## ###################################################
## Arg1, Retnla, Chil3 in Muscularis Macrophages
cellGroups = list(MuscularisMacrophages=c('MuscularisMacrophages1',
                                          'MuscularisMacrophages2',
                                          'MuscularisMacrophages3',
                                          'MuscularisMacrophages4'))
genes = c('Arg1','Retnla','Chil3')
quadrants = c('WT_infected','KO_infected')

MuscularisMacrophagesDFThree = getExpressionDF(f,cellGroups,genes,quadrants)
gMuscularisThree = makeThePlot(MuscularisMacrophagesDFThree,'Muscularis Macrophages')
dev.new()
print(gMuscularisThree)
ggsave(plot=gMuscularisThree,
       filename='violinPlots/MuscularisMacrophagesThree.pdf',
       height=8,width=15,units='in')

saveFigureDF(MuscularisMacrophagesDFThree,'violinPlots/MuscularisMacrophagesThree.pdf')

## ###################################################
## ###################################################
## Get the p-values:
pValues = function(df,figure)
{
    pValueDF = data.frame(figure=character(0),
                          group1=character(0),
                          group2=character(0),
                          pValue=numeric(0))

    comparisons = makeComparisons(df)
    for(i in 1:length(comparisons))
    {
        group1 = comparisons[[i]][1]
        group2 = comparisons[[i]][2]
        x = df[df$x==group1,]$expression
        y = df[df$x==group2,]$expression

        test = wilcox.test(x,y)

        a = data.frame(figure,
                       group1,
                       group2,
                       pValue=test$p.value)

        pValueDF = rbind(pValueDF,a)
    }

    return(pValueDF)
}

pValueDF = pValues(ifngDF,'Ifng')
pValueDF = rbind(pValueDF,
                 pValues(mesothelialDF,'Mesothelial'))

pValueDF = rbind(pValueDF,
                 pValues(FibroblastsDF,'Fibroblasts'))
                 
pValueDF = rbind(pValueDF,
                 pValues(MuscularisMacrophagesDF,'MuscularisMacrophages'))

pValueDF = rbind(pValueDF,
                 pValues(MuscularisMacrophagesDFTwo,'MuscularisMacrophagesTwo'))

pValueDF = rbind(pValueDF,
                 pValues(AllDF,'All'))

pValueDF = rbind(pValueDF,
                 pValues(MuscularisMacrophagesDFThree,'MuscularisMacrophagesThree'))

Write.Table(pValueDF,
            'violinPlots/violinPlotPValues.txt')

Write.Table(pValueDF,
            'figureTables/violinPlotPValues.txt')



                 
