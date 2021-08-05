
library(Seurat)
library(HandyPack)
library(tictoc)
library(ggplot2)
library(ggrepel)
library(stringr)
library(cowplot)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)
library(ggsci)


rm(list=ls())

source('SeuratWare02.R')

## ###################################################
## ###################################################
runThePipeline = function(assay='RNA')
{
    f = getSeuratObject()
    f@active.assay = assay
    
    countGroupsByQuadrant(f)
    basicPlots(f)
    plotQuadrantsSeparately(f)
    findGroupMarkerGenes(f)
    findDayZeroMarkers(f)
    DEResponseToInfection(f)
    KO_vs_WT_comparisons(f)
    differentialResponseToInfection(f)

    return('done!')
}

## ###################################################
prependGenes = function(df)
{
    df = cbind(data.frame(gene=rownames(df),
                          stringsAsFactors=FALSE),
               df)

    return(df)
}

## ###################################################
countGroupsByQuadrant = function(f)
{
    Tic('countGroupsByQuadrant')

    groups = getGroups()

    shortName = c()
    WT_uninfected = c()
    WT_infected = c()
    KO_uninfected = c()
    KO_infected = c()


    for(i in 1:length(groups))
    {
        ## The group name:
        shortName[i] = getShortName(groups[i])

        ## WT_uninfected
        idx = (f@meta.data$Seurat_Cluster == groups[i] &
               f@meta.data$genotype == 'WT' &
               f@meta.data$treatment == 'uninfected')
        WT_uninfected[i] = sum(idx)

        ## WT_infected
        idx = (f@meta.data$Seurat_Cluster == groups[i] &
               f@meta.data$genotype == 'WT' &
               f@meta.data$treatment == 'infected')
        WT_infected[i] = sum(idx)

        
        ## KO_uninfected
        idx = (f@meta.data$Seurat_Cluster == groups[i] &
               f@meta.data$genotype == 'KO' &
               f@meta.data$treatment == 'uninfected')
        KO_uninfected[i] = sum(idx)

        ## KO_infected
        idx = (f@meta.data$Seurat_Cluster == groups[i] &
               f@meta.data$genotype == 'KO' &
               f@meta.data$treatment == 'infected')
        KO_infected[i] = sum(idx)
    }

    countsDF = data.frame(group=groups,
                          shortName,
                          WT_uninfected,
                          WT_infected,
                          KO_uninfected,
                          KO_infected)

    tableDir = nameAndMakeDir('basicCounts')
    Write.Table(countsDF,
                paste(tableDir,'countsOfGroupVsQuadrant.txt',sep='/'))
                      
   toc() 
}

## ###################################################
denude = function(p)
{
    pPrime = p +
        theme(legend.position = "none",
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              title=element_blank())

    return(pPrime)
}

## ###################################################
makePlotDF = function(f)
{
    plotDF = data.frame(f@reductions$umap@cell.embeddings)
    plotDF$quadrant = f@meta.data$quadrant
    plotDF$Seurat_Cluster = f@meta.data$Seurat_Cluster
    plotDF$Seurat_Cluster = factor(plotDF$Seurat_Cluster)
    plotDF$shortName = f@meta.data$shortName

    tsne = data.frame(f@reductions$tsne@cell.embeddings)
    plotDF = cbind(plotDF,tsne)

    return(plotDF)
}

## ###################################################
basicPlots = function(f)
{
    Tic('basicPlots')
    plotDir = nameAndMakeDir('basicPlots')

    plotDF = makePlotDF(f)

    legendDotSize = 5


    n = length(unique(plotDF$Seurat_Cluster))
    palette = getThePalette()

    p1 = ggplot(plotDF,aes(x=UMAP_1,y=UMAP_2,color=shortName)) +
        geom_point(size=.3) +
        labs(color='') +
        ggtitle('UMAP colored by seurat cluster') +
        scale_color_manual(values=palette,breaks=getShortNames()) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize)))

    ggsave(plot=p1,
           filename=paste(plotDir,'umapCellType.jpg',sep='/'),
           dpi=600,
           height=7,width=12,units='in')
    ggsave(plot=p1,
           filename=paste(plotDir,'umapCellType.pdf',sep='/'),
           height=7,width=12,units='in')
    ggsave(plot=denude(p1),
           filename=paste(plotDir,'umapCellTypeNoFrills.pdf',sep='/'))
    

    p2 = p1 +
        facet_wrap(~ quadrant)
    ggsave(plot=p2,
           filename=paste(plotDir,'umapByCellTypeAndQuadrant.jpg',sep='/'),
           dpi=600,
           height=7,width=12,units='in')
    ggsave(plot=p2,
           filename=paste(plotDir,'umapByCellTypeAndQuadrant.pdf',sep='/'),
           height=7,width=12,units='in')
    ggsave(plot=denude(p2),
           filename=paste(plotDir,'umapByCellTypeAndQuadrantNoFrills.pdf',sep='/'))        

    
   p3 = ggplot(plotDF,aes(x=tSNE_1,y=tSNE_2,color=shortName)) +
        geom_point(size=.3) +
        labs(color='') +
        ggtitle('TSNE colored by Seurat_Cluster') +
        scale_color_manual(values=palette,breaks=getShortNames()) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize)))
    
    ggsave(plot=p3,
           filename=paste(plotDir,'tsneCellType.jpg',sep='/'),
           dpi=600,
           height=7,width=12,units='in')
    ggsave(plot=p3,
           filename=paste(plotDir,'tsneCellType.pdf',sep='/'),
           height=7,width=12,units='in')
    ggsave(plot=denude(p3),
           filename=paste(plotDir,'tsneCellTypeNoFrills.pdf',sep='/'))        

    p4 = p3 +
        facet_wrap(~ quadrant)
    ggsave(plot=p4,
           filename=paste(plotDir,'tsneByCellTypeAndQuadrant.jpg',sep='/'),
           dpi=600,
           height=7,width=12,units='in')
    ggsave(plot=p4,
           filename=paste(plotDir,'tsneByCellTypeAndQuadrant.pdf',sep='/'),
           height=7,width=12,units='in')
    ggsave(plot=denude(p4),
           filename=paste(plotDir,'tsneByCellTypeAndQuadrantNoFrills.pdf',sep='/'))        

    ## ###################################################
    p5 = ElbowPlot(f) +
        ggtitle('PCs')

    fileName = paste(plotDir,'/PCAElbowPlot.jpg',sep='')
    ggsave(plot=p5,
           filename=fileName,
           height=8,width=8,units='in')

    toc()
}

## ###################################################
plotQuadrantsSeparately = function(f)
{
    plotDF = makePlotDF(f)
    plotDir = nameAndMakeDir('quadrantPlots')

    n = length(unique(plotDF$Seurat_Cluster))
    palette = getThePalette()
    legendDotSize = 5

    for(q in unique(f@meta.data$quadrant))
    {
        subPlot = ggplot(plotDF[plotDF$quadrant==q,],aes(x=UMAP_1,y=UMAP_2,color=shortName)) +
            geom_point(size=.8) +
            labs(color='') +
            ggtitle('UMAP colored by seurat cluster') +
            scale_color_manual(values=palette,breaks=getShortNames()) +
            guides(color = guide_legend(override.aes = list(size=legendDotSize)))
        
        subPlot = denude(subPlot)
        
        fileName = paste0(plotDir,'/',
                          q,'_umap.pdf')
        ggsave(plot=subPlot,
               filename=fileName,
               height=8,width=8,units='in')
    }
}

## ###################################################
findGroupMarkerGenes = function(f)
{
    markerDir = nameAndMakeDir('groupMarkerGenes')
    groups = getGroups()

    ## ###################################################
    for(g in groups)
    {
        shortName = getShortName(g)
        Tic(paste('finding markers, group',g,shortName,Sys.time()))
        groupMarkerDF = FindMarkers(f,
                                    assay=f@active.assay,
                                    only.pos=TRUE,
                                    group.by=f@meta.data$Seurat_Cluster,
                                    ident.1=g,
                                    test='MAST')
        ## ###################################################
        ## Prepend the genes:
        groupMarkerDF = prependGenes(groupMarkerDF)

        ## ###################################################
        ## Subset and save:
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]

        if(nrow(groupMarkerDF) == 0)
            next
    
        fileName = paste(markerDir,
                         '/markerGenes_',
                         g,
                         '_',
                         shortName,
                         '.txt',sep='')

        Write.Table(groupMarkerDF,
                    fileName)
        toc()
    }
}

## ###################################################
findDayZeroMarkers = function(f)
{
    ## ###################################################
    ## Get the uninfected:
    idx = f@meta.data$treatment == 'uninfected'
    fPrime = f[,idx]

    groups = getGroups()
    markerDir = nameAndMakeDir('dayZero')

    for(g in groups)
    {
        shortName = getShortName(g)
        Tic(paste('finding day zero markers, group',g,shortName,Sys.time()))

        ## ###################################################
        ## Get this group:
        idx = fPrime@meta.data$Seurat_Cluster == g

        if(sum(idx) == 0)
            next
        
        fGroup = fPrime[,idx]

        idx = fGroup@meta.data$genotype == 'WT'
        if(sum(idx) < 3 |
           sum(!idx) < 3)
            next

        groupMarkerDF = FindMarkers(fGroup,
                                    assay=f@active.assay,
                                    group.by=fGroup@meta.data$genotype,
                                    ident.1='KO',
                                    ident.2='WT',
                                    test='MAST')

        ## ###################################################
        ## Prepend the genes:
        groupMarkerDF = prependGenes(groupMarkerDF)

        ## ###################################################
        ## Subset and save:
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]
    
        if(nrow(groupMarkerDF) == 0)
            next

        fileName = paste(markerDir,
                         '/markerGenesKOoverWTUninfected_',
                         g,
                         '_',
                         shortName,
                         '.txt',sep='')
        
        Write.Table(groupMarkerDF,
                    fileName)
        toc()
    }
}

## ###################################################
differentialResponseToInfection = function(f)
{
    ## ###################################################
    ## Get the infected:
    idx = f@meta.data$treatment == 'infected'
    fPrime = f[,idx]

    groups = getGroups()
    markerDir = nameAndMakeDir('differentialResponseToInfection')

    for(g in groups)
    {
        shortName = getShortName(g)
        Tic(paste('finding markers, differentialResponseToInfection group',g,shortName,Sys.time()))
        
        ## ###################################################
        ## Get this group:
        idx = fPrime@meta.data$Seurat_Cluster == g
        fGroup = fPrime[,idx]
        
        idx = fGroup@meta.data$genotype == 'WT'
        if(sum(idx) < 3 |
           sum(!idx) < 3)
            next
        
        groupMarkerDF = FindMarkers(fGroup,
                                    assay=f@active.assay,
                                    group.by=fGroup@meta.data$genotype,
                                    ident.1='KO',
                                    ident.2='WT',
                                    test='MAST')
        
        ## ###################################################
        ## Prepend the genes:
        groupMarkerDF = prependGenes(groupMarkerDF)
        
        ## ###################################################
        ## Subset and save:
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]

        if(nrow(groupMarkerDF) == 0)
            next
        
        fileName = paste(markerDir,
                         '/markerGenesKOoverWTInfected_',
                         g,
                         '_',
                         shortName,
                         '.txt',sep='')

        Write.Table(groupMarkerDF,
                    fileName)
        toc()
    }
}

## ###################################################
## Somehow this went missing!
## For each cell type it needs to find the differentially
## expressed genes in both the WT and KO cases.
## It needs to save the full table before subseting.
DEResponseToInfection = function(f)
{
    deDir = nameAndMakeDir('DE_ResponseToInfection')
    fullTableDir = nameAndMakeDir(c('DE_ResponseToInfection','fullTables'))

    groups = getGroups()
    genotypes = c('WT','KO')

    for(g in groups)
    {
        for(type in genotypes)
        {
            Tic(paste('finding marker genes DEResponseToInfection',
                      type,getShortName(g)))
            idx = (f@meta.data$genotype == type &
                   f@meta.data$Seurat_Cluster == g)
            
            fThis = f[,idx]

            idx = fThis@meta.data$treatment == 'infected'
            if(sum(idx) < 3 |
               sum(!idx) < 3)
                next
            
            markerDF = FindMarkers(fThis,
                                   assay=f@active.assay,
                                   group.by=fThis@meta.data$treatment,
                                   ident.1='infected',
                                   only.pos=FALSE,
                                   test='MAST')
            markerDF = prependGenes(markerDF)
            ## ###################################################
            ## Save the full table:
            shortName = getShortName(g)
            fileName = paste(fullTableDir,
                             '/deGenesInfOverUninf',
                             type,
                             'group',
                             g,
                             '_',
                             shortName,
                             '.txt.',sep='')
            Write.Table(markerDF,
                        fileName)
            ## ###################################################
            ## Now subset and save:
            cutoff = 0.05
            idx = markerDF$p_val_adj <= cutoff
            markerDF = markerDF[idx,]
            fileName = paste(deDir,
                             '/deGenesInfOverUninf',
                             type,
                             'group',
                             g,
                             '_',
                             shortName,
                             '.txt.',sep='')
            Write.Table(markerDF,
                        fileName)
            toc()  
        }
    }
}


## ###################################################
KO_vs_WT_comparisons = function(f)
{
    Tic('KO_vs_WT_comparisons')
    deDir = nameAndMakeDir('DE_ResponseToInfection')
    fullTableDir = nameAndMakeDir(c('DE_ResponseToInfection','fullTables'))
    plotDir = nameAndMakeDir(c(deDir,'lfcPlots'))
    comparisonTableDir = nameAndMakeDir(c(deDir,'comparisonTables'))
    groups = getGroups()
    vennDF = data.frame(group=groups,
                        shortName='',
                        WT_Only=NA,
                        KO_Only=NA,
                        WT_and_KO=NA,
                        stringsAsFactors=FALSE)

    
    for(i in 1:length(groups))
    {
        g = groups[i]
        shortName = getShortName(g)
        vennDF$shortName[i] = shortName
        longName = shortName
        ## ###################################################
        ## Get the two full tables:
        wtName = paste(fullTableDir,
                       '/deGenesInfOverUninf',
                       'WT',
                       'group',
                       g,
                       '_',
                       shortName,
                       '.txt',sep='')

        koName = paste(fullTableDir,
                       '/deGenesInfOverUninf',
                       'KO',
                       'group',
                       g,
                       '_',
                       shortName,
                       '.txt',sep='')

        ## ###################################################
        ## Do we have both files?
        if((! file.exists(wtName)) |
           (! file.exists(koName)))
            next
        
        wtDF = Read.Table(wtName)
        koDF = Read.Table(koName)

        ## ###################################################
        ## Relabel columns and merge:
        for(i in 2:ncol(wtDF))
        {
            names(wtDF)[i] = paste('wt',names(wtDF)[i],sep='_')
            names(koDF)[i] = paste('ko',names(koDF)[i],sep='_')
        }
        merged = merge(wtDF,koDF,by='gene')

        ## ###################################################
        ## Get the subset significant in at least one:
        cutoff = 0.05
        wtIdx = merged$wt_p_val_adj <= cutoff
        koIdx = merged$ko_p_val_adj <= cutoff
        merged$significantIn = ''
        merged$significantIn[wtIdx] = 'WT only'
        merged$significantIn[koIdx] = 'KO only'
        merged$significantIn[wtIdx & koIdx] = 'WT and KO'

        wtGenes = merged$gene[wtIdx]
        koGenes = merged$gene[koIdx]
        
        merged = merged[wtIdx | koIdx,]

        ## ###################################################
        ## Make log fold change plots:
        p = ggplot(merged,aes(x=wt_avg_logFC,y=ko_avg_logFC,color=significantIn)) +
            geom_point(size=.75) +
            ggtitle(paste('Group',g,longName))

        fileName = paste(plotDir,'/lfcKOvWT_group',g,'_',shortName,'.jpg',sep='')
        ggsave(plot=p,
               filename=fileName)

        ## ###################################################
        ## Get the venn diagram lists:
        vennDF$WT_Only[g] = length(setdiff(wtGenes,koGenes))
        vennDF$KO_Only[g] = length(setdiff(koGenes,wtGenes))
        vennDF$WT_and_KO[g] = length(intersect(wtGenes,koGenes))

        geneList = merged[,c('gene','significantIn','wt_avg_logFC','ko_avg_logFC')]
        fileName = paste(comparisonTableDir,'/DE_WT_vs_KO_Comparison_group',g,shortName,'.txt',sep='')
        Write.Table(geneList,
                    fileName)

    }
    fileName = paste(comparisonTableDir,'/DE_WT_vs_KO_VennCounts.txt',sep='')
    Write.Table(vennDF,
                fileName)
    toc()
}


