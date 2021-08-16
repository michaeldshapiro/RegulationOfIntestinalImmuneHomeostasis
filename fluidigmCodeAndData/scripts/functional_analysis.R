library(fgsea)

setwd("../output")

## read differentially expressed genes
de = read.csv("glia_DE.csv")

## read hallmark gene sets
pathways = gmtPathways("../data/h.all.v7.1.symbols.gmt")

## read mouse to human mapping
mapping = read.csv("../data/mouse_human_homologues.csv")
mouse_mapping = as.character(mapping$X10090)
names(mouse_mapping) = mapping$X9606

toMouse = function(human){
  return(as.character(mouse_mapping[human]))
}

## function from Michael for mapping pathway to mouse
allToMouse = function(pathways)
{
  pathwaysOut = list()
  for(n in names(pathways))
  {
    pathwaysOut[[n]] = unlist(lapply(pathways[[n]],toMouse))
    pathwaysOut[[n]] = pathwaysOut[[n]][!(is.na(pathwaysOut[[n]]))]
  }
  
  return(pathwaysOut)
}

## convert hallmark gene sets from human to mouse gene names
pathways_mouse = allToMouse(pathways)
de = read.csv("glia_condition_markers_mast_nFeature.csv")

## set seed for reproducibility and run FGSEA
set.seed(42)
de$avg_logFC = -de$avg_logFC
de = de[order(de$avg_logFC, decreasing=T),]
ranks=setNames(de$avg_logFC,de$X)
fgseaRes = fgsea(pathways_mouse, ranks, minSize=15, maxSize=2000, nperm=10000)
results=as.data.frame(fgseaRes)
results=results[order(results$pval),]
results$leadingEdge = vapply(results$leadingEdge, paste, collapse = ", ", character(1L))
write.csv(results,"fgsea_results_glia_hpoly_nFeature.csv",row.names=F)