
library(HandyPack)
library(ggplot2)
library(cowplot)

rm(list=ls())

graphics.off()

source('SeuratWare02.R')

## ###################################################
## ###################################################
denude = function(p)
{
    pPrime = p +
        theme(legend.position = "none",
              axis.text=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks=element_blank(),
              title=element_blank())

    return(pPrime)
}



## ###################################################
## Counts per cell type per quadrant:

basicCellCountDF = Read.Table('basicCounts/countsOfGroupVsQuadrant.txt')

countPlotDF = data.frame(group=numeric(0),
                         groupName=character(0),
                         count=numeric(0),
                         quadrant=character(0),
                         stringsAsFactors=FALSE)

## 3:6 are the numeric columns of the table:
for(i in 3:6)
{
    a = data.frame(group=basicCellCountDF$group,
                   shortName=basicCellCountDF$shortName,
                   count=basicCellCountDF[,i],
                   quadrant=names(basicCellCountDF)[i],
                   stringsAsFactors=FALSE)
    countPlotDF = rbind(countPlotDF,a)
}

quadrantOrder = getQuadrants()
countPlotDF$quadrant = factor(countPlotDF$quadrant,
                              levels=quadrantOrder)
shortNames = getShortNames()
countPlotDF$shortName = factor(countPlotDF$shortName,
                               levels=shortNames)

if(! dir.exists('paperFigures'))
    dir.create('paperFigures')

textSize = 15
textAngle = 60

palette = getThePalette()

Write.Table(countPlotDF,'figureTables/cellCountOfCellTypesPerQuadrant.txt')

## ###################################################
## We now do each plot by quadrant:
for(q in unique(countPlotDF$quadrant))
{
    g = ggplot(countPlotDF[countPlotDF$quadrant==q,]
              ,aes(x=shortName,y=count,fill=shortName)) +
        geom_col() +
        facet_wrap(~quadrant,nrow=2) +
        theme(axis.text.x = element_text(angle = textAngle,hjust=1,size=textSize),
              axis.text.y = element_text(size=textSize)) +
        theme(legend.position = "none") +
        xlab('Cell type') +
        ggtitle(paste('Counts per cell type',q)) +
        scale_fill_manual(values=palette,breaks=getShortNames())
    
    dev.new()
    print(g)
    ggsave(plot=g,
           filename=paste0('paperFigures/countsPerCellType_',q,'.pdf'),
           height=6,width=12,units='in')

    ggsave(plot=denude(g),
       filename=paste0('paperFigures/countsPerCellType_',q,'NoFrills.pdf'),
       height=6,width=12,units='in')
}

## ###################################################
## ###################################################
## Counts of DE genes wrt treatment and genotype per by
## quadrant:

## Uninfected WT v KO is in dayZero
## Infected   WT v KO is in differentialResponseToInfection
## WT uninfected v infected is in DE_ResponseToInfection
## KO uninfected v infected is in DE_ResponseToInfection

groups = getGroups()
numGroups = length(groups)

comparisons = c('Uninfected WT v KO',
                'Infected WT v KO',
                'WT Uninfected v Infected',
                'KO Uninfected v Infected')
deCountDF = data.frame(shortName=character(0),
                       count=numeric(0),
                       comparison=character(0),
                       stringsAsFactors=FALSE)
for(i in 1:4)
{
    a = data.frame(shortName=character(numGroups),
                       count=numeric(numGroups),
                       comparison=comparisons[i],
                       stringsAsFactors=FALSE)
    for(j in 1:length(groups))
    {
        group = groups[j]
        shortName = getShortName(group)
        a$shortName[j] = shortName
        if(i == 1)
        {
            fileName = paste0('dayZero/markerGenesKOoverWTUninfected_',
                              group,'_',shortName,'.txt')
        } else if(i == 2) {
            fileName = paste0('differentialResponseToInfection/markerGenesKOoverWTInfected_',
                              group,'_',shortName,'.txt')
        } else if(i == 3) {
            fileName = paste0('DE_ResponseToInfection/deGenesInfOverUninfWTgroup',
                              group,'_',shortName,'.txt.')
        } else if(i == 4) {
            fileName = paste0('DE_ResponseToInfection/deGenesInfOverUninfKOgroup',
                              group,'_',shortName,'.txt.')
        }
        if(! file.exists(fileName))
            next

        deDF = Read.Table(fileName)
        a$count[j] = nrow(deDF)
    }
    deCountDF = rbind(deCountDF,a)
}

deCountDF$shortName = factor(deCountDF$shortName,
                             levels=shortNames)



## ###################################################
## We're going to break this out by comparisons:
for(comp in comparisons)
{
    idx = deCountDF$comparison == comp 
    subDF = deCountDF[idx,]
    g = ggplot(subDF,aes(x=shortName,y=count,fill=shortName)) +
        geom_col() +
        ggtitle(comp) +
        theme(axis.text.x = element_text(angle = textAngle,hjust=1,size=textSize),
              axis.text.y = element_text(size=textSize)) +
        theme(legend.position = "none") +
        facet_wrap(~comparison) +
        xlab('Cell type') +
        scale_fill_manual(values=palette,breaks=getShortNames())

    dev.new()
    print(g)

    shortComp = str_replace(comp,' ','_')
    fileName = paste0('paperFigures/countsOfDEGenes',shortComp,'.pdf')
    ggsave(plot=g,
           filename=fileName,
           height=6,width=12,units='in')

    fileName = paste0('paperFigures/countsOfDEGenes',shortComp,'_NoFrills.pdf')
    ggsave(plot=denude(g),
           filename=fileName,
           height=6,width=12,units='in')
}

    

