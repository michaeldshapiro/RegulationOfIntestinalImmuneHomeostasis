
library(stringr)
library(HandyPack)
library(uniformly)

## ###################################################
#' Returns the Seurat object
#'
#' This retrieves the Seurat object and does some work to tidy things
#' up.
#' @return This returns the Seurat object.
#' @import Seurat
#' @import HandyPack
#' @import stringr
#' @param which - whether to use object 5 or 7
#' @export
getSeuratObject = function(which=7)
{
    stopifnot(which %in% c(5,7))
    
    if(which == 5)
    {
        versionDir = '~/FranzeSingleCell05' 
        fileName = paste0(versionDir,'/SeuratObject/SeuratObject.rds')
        f = readRDS(fileName)
        f@active.assay = 'RNA'


        return(f)
    }

    if(which == 7)
    {
        versionDir = '~/FranzeSingleCell07'
        fileName = paste0(versionDir,'/SeuratObject/SeuratObject.rds')
        SeuratObject = readRDS(fileName)

        return(SeuratObject)        
    }
    
}




## ###################################################
#' A handy little function for making directories
#'
#' A handy little function for making directories
#'
#' @param theDirs - a vector of nested directories to make in the
#'     current directory.
#' @return This returns the path to the deepest directory
#' @export
nameAndMakeDir = function(theDirs)
{
    thisDir = theDirs[1]

    if(! dir.exists(thisDir))
        dir.create(thisDir)

    if(length(theDirs) == 1)
        return(thisDir)

    for(i in 2:length(theDirs))
    {
        thisDir = paste(thisDir,theDirs[i],sep='/')
        if(! dir.exists(thisDir))
            dir.create(thisDir)
    }
    return(thisDir)
}

## ###################################################
#' This gets the groups, i.e., seurat_clusters as numbers
#' 
#' This gets the groups, i.e., seurat_clusters as numbers
#'
#' @param f - the Seurat object
#' @return This cluster numbers in numerical order
#' @export
getGroups = function()
{
    df = getGroupTable()
    return(df$Seurat_Cluster)
}

## ###################################################
#' This gets the group names, i.e., the names of the
#' Seurat clusters
#' 
#' This gets the groups names.  So far this is just the
#' numbers turned into characters.
#'
#' @param f - the Seurat object
#' @return This cluster names
#' @export
getGroupNames = function(f)
{
    groups = unique(f@meta.data$seurat_cluster)
    groups = groups[order(groups)]

    return(as.character(groups))
}

## ###################################################
#' Random color generator
#'
#' This produces a random set of colors for use as
#' a pallette randomly chosen from color names that aren't grey or
#' grey. 
#'
#' @param n - number of colors to generate
#' @param seed - set the seed for repeatability
#' @return a vector of hex values
#' @export
randomColors = function(n,seed=1)
{
    set.seed(seed)

    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

    palette = sample(color,n)
}

## ###################################################
#' This gets the palette for coloring cell types, e.g. in
#' UMAP plots
#'
#' @return A vector of colors as hex values
#' export
getThePalette = function(ncolor=30)
{
    return(get_random_distinct_colors(ncolor=30,seed=1))
}
    

## ###################################################
#' Random color generator
#'
#' This produces a random set of colors for use as
#' a pallette randomly chosen from the color cube.
#'
#' @param ncolor - the number of colors
#' @param seed - set the seed for repeatability
#' @return A list of colors as hex values
#' @export
get_random_distinct_colors = function(ncolor,seed = 100) {
  require(uniformly)
  set.seed(seed)
  rgb_mat <- runif_in_cube(n=ncolor,d=3,O=rep(0.5,3),r=0.5)
  rgb(r=rgb_mat[,1],g=rgb_mat[,2],b=rgb_mat[,3])
}

## ###################################################
#' Looks for a set of old cells in the new object
#'
#' @param oldCells - a subset of the colnames
#' @param newObje - the new Seurat object
#' @return a logical of length ncol(newObj)
#' @export
findOldCellsInNew = function(oldCells,newObj)
{
    oldCells = stripCellNames(oldCells)
    newCells = stripCellNames(colnames(newObj))

    where = c()
    for(i in 1:length(oldCells))
    {
        for(j in 1:length(newCells))
        {
            if(oldCells[i] == newCells[j])
            {
                where[i] = j
                break
            }
        }
    }
    return(where)
}

## ###################################################
#' A helper function to uniformize cell names
#'
#' @param cellNames - These are colnames from Seurat objects.
#' @return the cell names with trailing numbers removed
#' @export
stripCellNames = function(cellNames)
{
    cellNames = str_replace(cellNames,'-1','')
    a = str_split(cellNames,'_')
    cellNames = unlist(lapply(a,function(x) return(x[1])))

    return(cellNames)
}

## ###################################################
#' This retrieves the names of the cell clusters
#'
#' @return A character vector of the cluster names
#' @export
getShortNames = function()
{
    df = getGroupTable()
    return(df$shortName)
}

## ###################################################
#' This converts the seurat_cluster (or Seurat_Cluster)
#' number to the corresponding name
#'
#' @param g The number of a cluster
#' @return the name of the cluster
#' @export
getShortName = function(g)
{
    df = getGroupTable()
    idx = df$Seurat_Cluster == g
    stopifnot(sum(idx) == 1)

    return(df$shortName[idx])
}

## ###################################################
#' This is used to make the ancillary table that relates
#' group number to group name.
#'
#' This is called once when the Seurat object has its
#' clusters named and renumbered.
#'
#' @param f The Seurat object whose data we record.
#' @return none
#' @export
makeGroupTable = function(f)
{
    df = f@meta.data[,c('Seurat_Clusters','shortName')]
    df = df[! duplicated(df$Seurat_Clusters),]
    
    df = df[order(df$Seurat_Clusters),]
    
    Write.Table(df,
                '~/FranzeSingleCell07/groupTable.txt')
}

## ###################################################
#' This retrieves the group table produced by makeGroupTable
#'
#' @return A data frame
#' @export
getGroupTable = function()
{
    return(Read.Table('~/FranzeSingleCell07/groupTable.txt'))
}

## ###################################################
#' This retrieves the quadrants in proper order for plotting
#'
#' @return a character vector of the quadrants
#' @export
getQuadrants = function()
{
    quadrants = c('WT_uninfected',
                  'KO_uninfected',
                  'WT_infected',
                  'KO_infected')

    return(quadrants)
}

## ###################################################
#' This returns the quadrant colors for plotting
#'
#' @return a character vector of colors
#' @export
getQuadrantPalette = function()
{
    return(c('white','orange','gray','brown'))
}

