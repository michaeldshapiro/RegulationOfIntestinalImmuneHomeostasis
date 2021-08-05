
library(Seurat)
library(ggplot2)

rm(list=ls())
graphics.off()

source('SeuratWare02.R')

## ###################################################
analyseCluster = function(f,npcs=30,dims=1:20,res)
{
    f = CreateSeuratObject(f@assays$RNA@data)
    f = FindVariableFeatures(f)
    f = ScaleData(f)
    f = RunPCA(f,npcs=npcs)
    f = RunUMAP(f,reduction ="pca",dims=dims)
    ## f = RunTSNE(f,reduction ="pca",dims=dims)
    f = FindNeighbors(f,reduction='pca',dims=dims)

    f = FindClusters(f,res=.2)

    return(f)
}



## ###################################################
plotCellType = function(f,cells,n)
{
    ## Subset to uninfected:
    idx = f@meta.data$treatment == 'uninfected'
    ff = f[,idx]

    ## Subset to the cells:
    idx = ff@meta.data$shortName %in% cells
    fff = ff[,idx]
    fff@active.assay = 'RNA'

    F = analyseCluster(fff)

    df = data.frame(F@reductions[['umap']]@cell.embeddings)
    df$genotype = fff@meta.data$genotype

    g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=genotype)) +
        geom_point() +
        scale_color_manual(values=c('black','orange'),
                           breaks=c('WT','KO')) +
        theme(legend.position='none')

    Write.Table(df,
                paste0('figureTables/subclusteredUMAP_',
                       n,
                       '.txt'))

    return(g)
}

## ###################################################
plotSubClustersAndViolins = function(f,cells,cellName,gene)
{
    ## Subset to WT uninfected:
    idx = f@meta.data$quadrant == 'WT_uninfected'
    ff = f[,idx]
    
    ## Subset to the cells:
    idx = ff@meta.data$shortName %in% cells
    fff = ff[,idx]
    fff@active.assay = 'RNA'
    
    F = analyseCluster(fff)
    saveRDS(F,
            paste0('SeuratObject/subclustered',
                   cellName,
                   '.rds'))
            
    
    df = data.frame(F@reductions[['umap']]@cell.embeddings)
    df$genotype = fff@meta.data$genotype
    df$cluster = F@meta.data$seurat_clusters

     gClusters = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster)) +
        geom_point() 
    
    fileName = paste0('FranzeFigures/subclustered',
                      cellName,'.pdf')
    ggsave(plot=gClusters,
           filename=fileName)
    
    geneData = FetchData(F,gene)
    df$expression = geneData[,1]

    Write.Table(df,
                paste0('figureTables/subclustered',
                       cellName,
                       '_',
                       gene,
                       'Violin.txt'))
    

    gViolin = ggplot(df,aes(x=cluster,y=expression,color=cluster)) +
        geom_violin(draw_quantiles= c(0.25, 0.5, 0.75)) +
        geom_jitter(width=0.2)

    fileName = paste0('FranzeFigures/subclustered',
                      cellName,
                      '_',gene,
                      '.pdf')
    ggsave(plot=gViolin,
           filename=fileName)

}

## ###################################################
f = getSeuratObject()
    
cells = list(Fibroblasts=c('Fibroblasts1',
                           'Fibroblasts2'),
             MesothelialCells=c('MesothelialCells1',
                                'MesothelialCells2',
                                'MesothelialCells3',
                                'MesothelialCells4',
                                'MesothelialCells5'),
             MuscularisMacrophages=c('MuscularisMacrophages1',
                                     'MuscularisMacrophages2',
                                     'MuscularisMacrophages3',
                                     'MuscularisMacrophages4'))

figDir = 'FranzeFigures'
if(! dir.exists(figDir))
    dir.create(figDir)

gene = 'Pdgfra'

plotSubClustersAndViolins(f,cells$Fibroblasts,'Fibroblasts',gene)

for(n in names(cells))
{
    g = plotCellType(f,cells[[n]],n)
    fileName = paste0(figDir,'/',n,'_cluster.pdf')
    ggsave(plot=g,
           filename=fileName)
    
    dev.new()
    print(g + ggtitle(n))
}

