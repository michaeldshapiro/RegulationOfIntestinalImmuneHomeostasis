library(Seurat)
library(tictoc)
library(stringr)
library(ggplot2)
library(HandyPack)

source('SeuratWare02.R')

## ###################################################
## ###################################################

ff = getSeuratObject()


## ###################################################
twoGenesScatterPlot = function(gene1,gene2,cellTypes,facet)
{
    ## ###################################################
    ## Subset to cellTypes:
    idx = ff@meta.data$shortName %in% cellTypes
    f = ff[,idx]
    
    ## ###################################################
    ## Get the gene expression:
    gene1Expression = as.matrix(f[['RNA']]@data)[gene1,]
    gene2Expression = as.matrix(f[['RNA']]@data)[gene2,]
    
    ## ###################################################
    ## Make the plot:    
    scatterDF = data.frame(gene1=gene1Expression,
                           gene2=gene2Expression,
                           cellType=f@meta.data$shortName,
                           quadrant=f@meta.data$quadrant,
                           stringsAsFactors = FALSE)
    
    cells = paste(cellTypes,collapse=', ')
    titleStr = paste(gene1,gene2,'expression in',cells)
    
    h = ggplot(scatterDF,aes(x=gene1,y=gene2,color=cellType)) +
        geom_point(size=2) +
        ggtitle(titleStr) +
        scale_color_manual(values=getThePalette(),
                           breaks=getShortNames()) +
        xlab(gene1) + ylab(gene2)

    if(facet)
        h = h + facet_wrap(~ quadrant) 
    
    dev.new()
    print(h)
    
    figDir = 'twoGeneFigures'
    if(! dir.exists(figDir))
        dir.create(figDir)
    
    cells = paste(cellTypes,collapse='_')
    if(! facet)
    {
        fileName = paste0(figDir,'/',
                          gene1,
                          '_',
                          gene2,
                          '_',
                          'ExpressionIn',
                          '_',
                          cells,
                          '.pdf')
    } else {
        fileName = paste0(figDir,'/',
                          gene1,
                          '_',
                          gene2,
                          '_',
                          'ExpressionIn',
                          '_',
                          cells,
                          '_faceted',
                          '.pdf')
    }
    
    ggsave(plot=h,
           filename=fileName)

    ## Now save the figure table:
    tableName = str_replace(fileName,figDir,'figureTables')
    tableName = str_replace(tableName,'pdf','txt')

    Write.Table(scatterDF,
                tableName)
    

    ## ###################################################
    ## Make the table of counts:
    countCases = function(df,cellTypes,gene1,gene2)
    {
        M = matrix(0,nrow=length(cellTypes),ncol=4)
        rownames(M) = cellTypes
        colnames(M) = c('both',gene1,gene2,'neither')

        for(cell in cellTypes)
        {
            idx = df$cellType == cell
            DF = df[idx,]
            M[cell,1] = sum(DF$gene1 > 0 &
                            DF$gene2 > 0)
            M[cell,2] = sum(DF$gene1 > 0 &
                            DF$gene2 == 0)
            M[cell,3] = sum(DF$gene1 == 0 &
                            DF$gene2 > 0)
            M[cell,4] = sum(DF$gene1 == 0 &
                            DF$gene2 == 0)            
        }
        countDF = data.frame(cellTypes,
                             M,
                             stringsAsFactors=FALSE)

        return(countDF)
    }

    if(!facet)
        countDF = countCases(scatterDF,cellTypes,gene1,gene2)
    
    if(facet)
    {
        quadrants = getQuadrants()
        
        for(i in 1:4)
        {
            q = quadrants[i]
            idx = scatterDF$quadrant == q
            subDF = scatterDF[idx,]
            subCount = countCases(subDF,cellTypes,gene1,gene2)
            subCount$cellTypes = paste(subCount$cellTypes,q)
            
            if(i == 1)
            {
                countDF = subCount
            } else {
                countDF = rbind(countDF,subCount)
            }
        }
    }
    
    fileName = str_replace(fileName,'ExpressionIn','ExpressionCount')
    fileName = str_replace(fileName,'\\.pdf','.txt')
    Write.Table(countDF,
                fileName)
}

## ###################################################
## ###################################################
