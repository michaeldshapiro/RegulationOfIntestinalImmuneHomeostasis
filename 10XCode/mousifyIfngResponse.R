
library(BioHandy)
library(stringr)
library(HandyPack)

genes = readLines('gseaIfngResponse.txt')
genes = str_split(genes,'\t')
genes = unlist(genes)

mouseGenes = humanToMouse(genes)

Write.Table(data.frame(genes=mouseGenes),
                      'IfngResponseGenes.txt')
