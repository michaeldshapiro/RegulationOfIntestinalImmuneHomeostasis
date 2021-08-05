
library(HandyPack)
library(Seurat)

rm(list=ls())

source('SeuratWare02.R')

## ###################################################
## ###################################################

subsample = function(f,shortName,theseQuadrants,upper=100,lower=50)
{
    ## Restrict to these cells and quadrants:
    idx = (f@meta.data$shortName == shortName &
           f@meta.data$quadrant %in% theseQuadrants)
    ff = f[,idx]

    ## Subsample:
    idx = ff@meta.data$quadrant == theseQuadrants[1]
    ff1 = ff[,idx]
    ff2 = ff[,!idx]
    N1 = sum(idx)
    N2 = sum(!idx)

    if(N1 >= N2)
    {
        ff1 = ff1[,sample(N1,upper)]
        ff2 = ff2[,sample(N2,lower)]
    } else {
        ff1 = ff1[,sample(N1,lower)]
        ff2 = ff2[,sample(N2,upper)]
    }
    ff = merge(ff1,ff2)

    return(ff)
}

## ###################################################
getSubsampledDEGenes = function(f,shortName,quadrants,cutoff=0.05)
{
    ## Subsample:
    ff = subsample(f,shortName,quadrants)

    ## Get the DE genes:
    markerDF = FindMarkers(ff,
                           assay='RNA',
                           group.by=ff@meta.data$quadrant,
                           ident.1=quadrants[2],
                           only.pos=FALSE,
                           test='MAST',
                           latent.vars='nCount_RNA')

    ## Prepend the genes:
    ## markerDF = cbind(data.frame(gene=rownames(ff),
    ##                             stringsAsFactors=FALSE),
    ##                  markerDF)

    ## Subset by p_val_adj:
    idx = markerDF$p_val_adj <= cutoff
    markerDF = markerDF[idx,]

    ## Sort by avg_logFC:
    markerDF = markerDF[order(-markerDF$avg_log2FC),]

    return(markerDF)
}

## ###################################################
## ###################################################

df = Read.Table('basicCounts/countsOfGroupVsQuadrant.txt')
f = getSeuratObject()

quadrants = getQuadrants()
shortNames = getShortNames()

comparisons = list(WT=c(1,2),
                   KO=c(3,4),
                   uninfected=c(1,3),
                   infected=c(2,4))

lower = 50
upper = 100

theShortNames = c()
first = c()
second = c()
finger = 1

subDir = nameAndMakeDir('subsampling')

for(i in 1:4)
{
    theseQuadrants = quadrants[comparisons[[i]]]
    

    for(shortName in shortNames)
    {
        idx = df$shortName == shortName
        these = df[idx,theseQuadrants]

        ## ###################################################
        ## Can we sub-sample:
        if(min(these) > lower &
           max(these) > upper)
        {
            ## Record the comparisons:
            theShortNames[finger] = shortName
            first[finger] = theseQuadrants[1]
            second[finger] = theseQuadrants[2]
            finger = finger + 1

            ## Get the DE genes:
            deGenes = getSubsampledDEGenes(f,shortName,theseQuadrants)

            ## Save these:
            fileName = paste0('subsampling/subsampledDEGenes_nCount_',
                              shortName,
                              '_',
                              theseQuadrants[2],
                              '_over_',
                              theseQuadrants[1],
                              '.txt')
            Write.Table(deGenes,
                        fileName)
        }
    }
}

canSubsample = data.frame(shortName=theShortNames,
                          quadrant1=first,
                          quadrant2=second,
                          stringsAsFactors=FALSE)

Write.Table(canSubsample,
            'subsampling/forSubsampling_nCount.txt')


