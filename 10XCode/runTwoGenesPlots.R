
rm(list=ls())
graphics.off()

source('twoGenesScatterPlot.R')

## ###################################################

cells = c('OtherLymphoidCells',
          'NaturalKillerCells',
          'TCells1',
          'TCells2',
          'TCells3')

genePairs = list()
genePairs[[1]] = c('Cd8a','Ifng')
genePairs[[2]] = c('Cd4','Ifng')
genePairs[[3]] = c('Cxcr3','Cd8a')

for(i in 1:length(genePairs))
{
    twoGenesScatterPlot(genePairs[[i]][1],
                        genePairs[[i]][2],
                        cells,
                        facet=FALSE)
}

