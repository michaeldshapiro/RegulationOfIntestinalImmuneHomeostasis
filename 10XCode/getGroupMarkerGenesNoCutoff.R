
source('SeuratWare02.R')

## ###################################################
prependGenes = function(df)
{
    df = cbind(data.frame(gene=rownames(df),
                          stringsAsFactors=FALSE),
               df)

    return(df)
}


findGroupMarkerGenesNoCutoff = function()
{
    f = getSeuratObject()
    markerDir = nameAndMakeDir('groupMarkerGenesForGSEA')
    groups = getGroups()

    ## ###################################################
    for(g in groups)
    {
        shortName = getShortName(g)
        Tic(paste('finding markers, group',g,shortName,Sys.time()))
        groupMarkerDF = FindMarkers(f,
                                    assay='RNA',
                                    group.by=f@meta.data$seurat_cluster,
                                    ident.1=g,
                                    test='MAST')
        ## ###################################################
        ## Prepend the genes:
        groupMarkerDF = prependGenes(groupMarkerDF)

        ## ###################################################
        ## Save:
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
findGroupMarkerGenesNoCutoff()

## ###################################################
## Oops, I forgot to order these by lfc:

files = Sys.glob('groupMarkerGenesForGSEA/*txt')

for(f in files)
{
    df = Read.Table(f)
    df = df[order(-df$avg_log2FC),]

    Write.Table(df,f)
}

