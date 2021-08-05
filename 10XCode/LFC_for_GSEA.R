
library(Seurat)
library(HandyPack)
library(tictoc)

rm(list=ls())

source('SeuratWare02.R')

## ###################################################
## ###################################################
    if(! dir.exists('GSEA'))
        dir.create('GSEA')

## We retrieve the desired cell types and comparisons
## from GSEA_cellTypes.csv
cellTypes = Read.Table('GSEA_cellTypes.csv')

## ###################################################
## Get the object:
f = getSeuratObject()

treatments = c('uninfected','infected')

for(i in 1:nrow(cellTypes))
{
    treatment = cellTypes$treatment[i]
    gseaDir = nameAndMakeDir(paste0('GSEA/',treatment))
    
    ## ###################################################
    ## Subset to the treatment:
    idx = f@meta.data$treatment == treatment
    ff = f[,idx]
    
    shortName = cellTypes$shortName[i]
    ## ###################################################
    fileName = paste0(gseaDir,'/LFC_',
                      treatment,'_',
                      shortName,
                      '.rnk')
    if(file.exists(fileName))
        next
    
    Tic(paste(treatment,shortName))
    ## ###################################################
    ## Subset to group:
    idx = ff@meta.data$shortName == shortName
    fff = ff[,idx]
    
    ## ###################################################
    ## Compute LFC:
    markerDF = FindMarkers(fff,
                           assay='RNA',
                           group.by=fff@meta.data$genotype,
                           ident.1='KO',
                           ident.2='WT',
                           test='MAST',
                           logfc.threshold=0)
    
    ## ###################################################
    ## For sanity's sake:
    markerDF = markerDF[!is.na(markerDF$avg_log2FC),]
    
    ## ###################################################
    ## Order by LFC:
    markerDF = markerDF[order(-markerDF$avg_log2FC),]
    
    
    asRank = data.frame(ID=rownames(markerDF),
                        t=markerDF$avg_log2FC,
                        stringsAsFactors=FALSE)
    
    Write.Table(asRank,
                fileName)
    toc()
}

