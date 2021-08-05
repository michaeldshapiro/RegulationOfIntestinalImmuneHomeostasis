
library(HandyPack)
library(ggplot2)
library(stringr)

rm(list=ls())
graphics.off()

source('SeuratWare02.R')

## ###################################################
## ###################################################


quadrants = getQuadrants()
shortNames = getShortNames()

comparisons = list(WT=c(1,2),
                   KO=c(3,4),
                   uninfected=c(1,3),
                   infected=c(2,4))

subDir = nameAndMakeDir('subsampling')


theShortNames = c()
first = c()
second = c()
count = c()

finger = 1

for(i in 1:4)
{
    theseQuadrants = quadrants[comparisons[[i]]]
    
    for(shortName in shortNames)
    {
        fileName = paste0('subsampling/subsampledDEGenes_nCount_',
                          shortName,
                          '_',
                          theseQuadrants[2],
                          '_over_',
                          theseQuadrants[1],
                          '.txt')
        if(! file.exists(fileName))
            next
        
        deGenes = Read.Table(fileName)

        theShortNames[finger] = shortName
        first[finger] = theseQuadrants[1]
        second[finger] = theseQuadrants[2]
        count[finger] = nrow(deGenes)

        finger = finger + 1
    }
}

deCounts = data.frame(shortName=theShortNames,
                      quadrant1=first,
                      quadrant2=second,
                      DEGenes=count,
                      stringsAsFactors=FALSE)
## Write.Table(deCounts,
##             'subsampling/countControlledDEGenes_nCount.txt')

deCounts$shortName = factor(deCounts$shortName,levels=shortNames)
deCounts$comparison = paste(deCounts$quadrant2,
                            'vs',
                            deCounts$quadrant1)

g = ggplot(deCounts,aes(x=shortName,y=DEGenes,fill=shortName)) +
    geom_col() +
    facet_wrap(~comparison) +
    ggtitle('DE genes controlled for cell counts') +
    theme(legend.position='none',
          axis.text.x = element_text(angle = 60,hjust=1,size=15)) +
    scale_fill_manual(values=getThePalette(),
                      breaks=getShortNames())

## dev.new()
## print(g)
## ggsave(plot=g,
##        'subsampling/countControlledDEGenes_nCount.pdf')


for(comp in unique(deCounts$comparison))
{
    thisDF = deCounts[deCounts$comparison==comp,]

    a = data.frame(shortName=getShortNames(),
                   quadrant1=thisDF$quadrant1[1],
                   quadrant2=thisDF$quadrant2[1],
                   DEGenes=0,
                   comparison=comp)
    idx = a$shortName %in% thisDF$shortName
    thisDF = rbind(thisDF,a[!idx,])

    gSub = ggplot(thisDF,aes(x=shortName,y=DEGenes,fill=shortName)) +
    geom_col() +
    ggtitle(comp) +
    theme(legend.position='none',
          axis.text.x = element_text(angle = 60,hjust=1,size=15)) +
    scale_fill_manual(values=getThePalette(),
                      breaks=getShortNames())

    dev.new()
    print(gSub)

    fileName = paste0('subsampling/countControlledGenes_',
                      comp,
                      '.pdf')

    fileName = str_replace_all(fileName,' ','_')
    ggsave(plot=gSub,
           filename=fileName,
           height=8,width=12,units='in')

    Write.Table(thisDF,
                paste0('figureTables/subsampledGeneCount_',
                       comp,
                       '.txt'))
}


    
