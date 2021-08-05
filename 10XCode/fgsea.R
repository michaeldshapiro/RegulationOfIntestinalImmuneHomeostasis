
library(fgsea)
library(HandyPack)
library(BioHandy)
library(stringr)

rm(list=ls())

## ###################################################
trimPathways = function(pathways)
{
    for(i in 1:length(pathways))
        pathways[[i]] = pathways[[i]][pathways[[i]] != '']

    return(pathways)
}

## ###################################################
toMouse = function(human)
{
    n = nchar(human)
    if(n == 1)
        return(human)

    mouse = paste0(substr(human,1,1),
                   str_to_lower(substr(human,2,n)))

    return(mouse)
}

## ###################################################
allToMouse = function(pathways)
{
    pathwaysOut = list()
    for(n in names(pathways))
    {
        pathwaysOut[[n]] = unlist(lapply(pathways[[n]],toMouse))
    }

    return(pathwaysOut)
}

## ###################################################
deListLE = function(ell)
{
    N = length(ell)
    b = character(N)

    for(i in 1:N)
        b[i] = paste(unlist(ell[[i]]),collapse=', ')

    return(b)
}


## ###################################################
getPathways = function()
{
    file = '~/geneSets/mouse/h.all.v7.2.symbols.gmt'
    pathways = gmtPathways(file)

    return(pathways)
}


## ###################################################
## ###################################################

pathwayTags = c('HallmarkGeneSets')

pathways = getPathways()

rnkFiles = Sys.glob('GSEA/*/*.rnk')

padjCutoff = .05

for(rnk in rnkFiles)
{
    ranksDF = Read.Table(rnk)
    ranks = ranksDF$t
    names(ranks) = ranksDF$ID
    
    
    fgseaRes = fgsea(pathways, ranks, minSize=15, maxSize=500) 
    fgseaRes = fgseaRes[order(fgseaRes$padj),]
    
    fgseaRes$leadingEdge = deListLE(fgseaRes$leadingEdge)
    
    tag = paste0('_',pathwayTags,'_All.txt')
    fOut = str_replace(rnk,'LFC_','')
    fOut = str_replace(fOut,'\\.rnk',tag)
    
    print(fOut)
    
    Write.Table(fgseaRes,
                fOut)   
    
    idx = fgseaRes$padj < padjCutoff
    fgseaRes = fgseaRes[idx,]
    if(nrow(fgseaRes) == 0)
        next
    
    fOut = str_replace(fOut,'_All','')
    Write.Table(fgseaRes,
                fOut)   
}

    
