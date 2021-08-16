library(Seurat)
library(dplyr)
library(scales)
library(biomaRt)
library(plyr)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(Seurat)
library(gplots)
library(scales)
library(ggplot2)
library(fgsea)
library(viridis)
library(pheatmap)

## set working directory
setwd('../data')

## read data from infected glia
HP1 = read.csv('HP_1_counts.tsv', sep = '\t')
HP2 = read.csv('HP_2_counts.tsv', sep = '\t')
HP3 = read.csv('HP_3_counts.tsv', sep = '\t')

## read data from control glia
C = read.csv('C_counts.tsv', sep = '\t')

## merge infected and control data
all = merge(HP1,HP2, by = "gene")
all = merge(all,HP3, by = "gene")
all = merge(all,C, by = "gene")
tail(all$gene)

## read info on control samples 
sample_info = read.csv("project-PRO472-sample-sheet.csv")
sample_info$Sample.Name

## function to get info on cell number/shape from sample name
get_cell_num = function(x){
  return(strsplit(x,'_')[[1]][1])
}

## check annotations on cell number/shape
sample_info$n_cell = as.character(lapply(as.character(sample_info$Sample.Name),get_cell_num))
unique(sample_info$n_cell)

## retain only single cells 
failed = subset(sample_info,!(sample_info$n_cell %in% c("1","1S","1?","1Irr","1+D","1L","1B")))
failed_samples = substring(failed$Sample.LIMSID,7)
failed_samples = paste0("C_",failed_samples)

mat = as.matrix(all[,2:ncol(all)])
rownames(mat) = all$gene
colnames(mat) = colnames(all)[2:ncol(all)]
mat = mat[,!(colnames(mat) %in% failed_samples)]


## add batch information
batch_1 = subset(sample_info,sample_info$Container.Name %in% c("FP1_FL020817","FP2_FL020817"))
batch_2 = subset(sample_info,sample_info$Container.Name %in% c("FP1_FL140817","FP2_FL140817"))
batch_1_samples = substring(batch_1$Sample.LIMSID,7)
batch_2_samples = substring(batch_2$Sample.LIMSID,7)
batch_1_samples = paste0("C_",batch_1_samples)
batch_2_samples = paste0("C_",batch_2_samples)

add_batch = function(x,batch){
  info = strsplit(x,'_')[[1]]
  return(paste(info[1],batch,info[2],sep = '_'))
}

colnames(mat)[colnames(mat) %in% batch_1_samples] = 
  as.character(lapply(colnames(mat)[colnames(mat) %in% batch_1_samples], add_batch, "1"))

colnames(mat)[colnames(mat) %in% batch_2_samples] = 
  as.character(lapply(colnames(mat)[colnames(mat) %in% batch_2_samples], add_batch, "2"))


## remove spike ins
is.ercc = which(grepl("^ERCC", rownames(mat)))
mat = mat[!(grepl("^ERCC", rownames(mat))),]


## remove macrophages
sum(mat["ENSMUSG00000030786",] >10) #7 removed
mat = mat[,!(mat["ENSMUSG00000030786",] >10)]


## remove neurons
hist(mat["ENSMUSG00000028546",],xlim = c(1,200), breaks = 2000)
sum(mat["ENSMUSG00000028546",] >30) #95

#"ENSMUSG00000030110" Ret
mat = mat[,!(mat["ENSMUSG00000028546",] >30)]


mat = mat[rowSums(mat) > 0,]

## get gene names
ensembl = useMart("ensembl", dataset='mmusculus_gene_ensembl')
gene_names_ids.df=getBM(c("ensembl_gene_id", "external_gene_name"), mart= ensembl)
gene_id = as.data.frame(rownames(mat))
names(gene_id) = "ensembl_gene_id"
gene_id = join(gene_id,gene_names_ids.df, by = "ensembl_gene_id")
gene_id$external_gene_name[is.na(gene_id$external_gene_name)] = gene_id$ensembl_gene_id[is.na(gene_id$external_gene_name)]

## make gene names unique
gene_id$external_gene_name = make.unique(gene_id$external_gene_name)
gene_id$external_gene_name
rownames(mat) = gene_id$external_gene_name

## create seurat object
data = CreateSeuratObject(counts = mat)

## define mitochondrial genes
location = select(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=as.character(gene_id$ensembl_gene_id),
                  column="CDSCHROM", keytype="GENEID")
location = location[match(as.character(gene_id$ensembl_gene_id), location$GENEID),]
is.mito = location$CDSCHROM == "chrM" & !is.na(location$CDSCHROM)
data[["percent.mt"]] = PercentageFeatureSet(data, features = rownames(data)[is.mito])

## change to output directory
setwd("../output")

## plot quality metrics
pdf("quality_metrics.pdf")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2,nrow = 2)
dev.off()

plot1 = FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


## dim(data)
#[1] 31272   394

## filter for > 2000 genes and less than 10 % mitochondrial read counts
data = subset(data, subset = nFeature_RNA > 2000 & percent.mt < 10)
###dim(data)
#[1] 31272   217

## filter out genes with less than 10 counts
data = data[rowSums(data) >= 10,]
#dim(data)
#[1] 20603   217

## filter data for genes expressed in at least 3 cells
n_cells_exp_gene = rowSums(GetAssayData(data) > 0)
data = data[n_cells_exp_gene >= 3,]
#dim(data)
#[1] 18643   217

## normalise data
data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

## find variable features
data = FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
# 
## add condition and batch information
data$condition = substr(colnames(data),1,1)
data$batch = substr(colnames(data),1,4)

## scale data
all.genes = rownames(data)
data = ScaleData(data, features = all.genes)

## run PCA
data = RunPCA(data, features = VariableFeatures(object = data))

## evaluate number of PCs to use, diagnostic plots
data = JackStraw(data, num.replicate = 1000)
data = ScoreJackStraw(data, dims = 1:20)
pdf("choice_n_pcs.pdf")
JackStrawPlot(data, dims = 1:20)
ElbowPlot(data)
dev.off()

## plot UMAPs constructed with different numbers of PCs
dir.create("seed_1")
for (i in 4:20){
  data= FindNeighbors(data, dims = 1:i)
  data = FindClusters(data, resolution = 0.1)
  data = RunUMAP(data, dims = 1:i, seed.use = 1)
  pdf(paste0("./seed_1/umap_louvain_seed_1_", i, ".pdf"))
  p = DimPlot(data, reduction = "umap", 
              group.by = c("RNA_snn_res.0.1","condition","batch"))
  print(p)
  dev.off()
}


## plot umaps with different seeds
dir.create("10_pcs")
for (i in 1:20){
  data = FindNeighbors(data, dims = 1:10)
  data = FindClusters(data, resolution = 0.1)
  data = RunUMAP(data, dims = 1:10, seed.use = i)
  pdf(paste0("./10_pcs/umap_louvain_10pcs_seed_", i, ".pdf"))
  p = DimPlot(data, reduction = "umap", 
              group.by = c("RNA_snn_res.0.1","condition","batch"))
  print(p)
  dev.off()
}

## plot heatmaps for genes associated with each PC
pdf('dim_heatmaps.pdf')
for (i in 1:30){
  DimHeatmap(object = data, dims = i, reduction = "pca", cells = 200, balanced = TRUE)
}
dev.off()

## decision to use 4 PCs, calculate neighbour graph and UMAP using these.
data = FindNeighbors(data, dims = 1:4)
data = FindClusters(data)
data = RunUMAP(data, dims = 1:4, seed.use = 1)

## Plot UMAP
DimPlot(data, reduction = "umap", 
        group.by = c("RNA_snn_res.0.8","condition","batch"))
DimPlot(data, reduction = "umap")

## save seurat object
saveRDS(data, "seurat_object.RDS")
data = readRDS("seurat_object.RDS")

## function to perform Louvain clustering at a given resolution 
## and to create heatmaps of marker genes for each cluster
dir.create("cluster_resolutions")
cluster = function(data,resolution){
  data = FindClusters(data,resolution = resolution)
  resolution = paste0("RNA_snn_res.",resolution)
  pdf(paste0("./cluster_resolutions/umap_louvain_",resolution,'.pdf'))
  p = DimPlot(data, reduction = "umap",group.by = "seurat_clusters")
  print(p)
  dev.off()
  
  markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers_lFC1 = subset(markers,markers$avg_logFC > 1)
  top10 = markers_lFC1 %>% group_by(cluster) %>% top_n(-10, p_val)
  top10 = c(top10$gene)
  
  pdf(paste0("./cluster_resolutions/seurat_heatmap_",resolution,".pdf"))
  p = DoHeatmap(data, features = top10)
  print(p)
  dev.off()
  
  mat = GetAssayData(object = data)
  mat = as.matrix(mat[top10,])
  mat = mat[,order(data[["seurat_clusters"]])]
  
  # Create vector with levels of object@ident
  identities = levels(data@active.ident)
  
  # Create vector of default ggplot2 colors
  my_color_palette = hue_pal()(length(identities))
  
  mat_norm = mat - rowMeans(mat)
  
  color.palette = colorRampPalette(c("blue", "black", "yellow"))
  pdf(paste0("cluster_resolutions/seurat_heatmap2_",resolution,".pdf"))
  p = heatmap.2(mat_norm,Colv = F,Rowv = T, col=color.palette, symbreak=TRUE, trace='none', 
                ColSideColors=my_color_palette[data[["seurat_clusters"]][["seurat_clusters"]][order(data[["seurat_clusters"]])]],cexCol = 0.25)
  print(p)
  dev.off()
}


## iterate through different resolutions 
## and perform clustering
for (res in c(0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2)){
  cluster(data,res)
}

#decision to use resolution 0.1
data = FindClusters(data,resolution = 0.1)

data$cluster = ""
data$cluster[data$RNA_snn_res.0.1 == 0] = "Glia_1"
data$cluster[data$RNA_snn_res.0.1 == 1] = "Glia_2"

data@active.ident = factor(data$cluster,levels = c("Glia_1","Glia_2"))

## define colours for clusters
cluster_cols = c("#0072B2","#F0E442")

## plot UMAP coloured by cluster
#stored as data$cluster (active identity)
pdf("umap.pdf")
p = DimPlot(data, reduction = "umap", pt.size = 3,
            label = F,label.size=5, cols = cluster_cols)
p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) 
print(p)
dev.off()

#label conditions
data$condition_label = ""
data$condition_label[data$condition == "C"] = "Control"
data$condition_label[data$condition == "H"]= "H. poly"
data$condition_label = as.factor(data$condition_label)

#define colours for condition
condition_cols = c("#009E73", "#E69F00")

## plot UMAP coloured by condition
#stored as data$condition_label
pdf("condition_umap.pdf")
p = DimPlot(data, reduction = "umap", pt.size = 3, label = F, label.size=5,
            group.by = "condition_label", cols = condition_cols) 
p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) 
print(p)
dev.off()

#extract x and y limits for consistency in feature plots
ylim = ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
xlim = ggplot_build(p)$layout$panel_scales_x[[1]]$range$range

##make UMAP showing interferon response genes
terms = gmtPathways("../data/mmusculus.GO:BP.name.gmt")

#get those interferon response genes which are in the data
interferon_genes = rownames(data)[rownames(data) %in% terms[['GO:0034341']]]

#write selected genes to file
write.csv(interferon_genes, "interferon_genes_all.csv", row.names = F, quote = F)

#get mean scaled expression of interferon genes per cell
data$response_to_ifn_scaled =  colMeans(x = GetAssayData(object = data, slot = "scale.data")[interferon_genes, ])

pdf("interferon_response_all_genes_scaled.pdf")
p = FeaturePlot(object = data, features = "response_to_ifn_scaled", pt.size = 3,
                label = F)
p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_blank()) + scale_color_viridis() +
  xlim(xlim) + ylim(ylim)
print(p)
dev.off() 


## make the same plot but without scaling of genes

#get mean expression of interferon genes per cell
data$response_to_ifn =  colMeans(x = GetAssayData(object = data)[interferon_genes, ])

pdf("interferon_response_all_genes.pdf")
p = FeaturePlot(object = data, features = "response_to_ifn", pt.size = 3,
                label = F)
p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_blank()) + scale_color_viridis() +
  xlim(xlim) + ylim(ylim)
print(p)
dev.off() 

## Violin plot of known markers
markers = c("Gfap","Plp1","Sox10","S100b","Erbb3","Fabp7")
values = GetAssayData(data)
results = data.frame(expression = numeric(), cluster = character(), gene = character())

for (marker in markers){
  glia_1_vals = values[rownames(data) == marker, data$cluster == "Glia_1"]
  glia_1_res = data.frame(expression = as.numeric(glia_1_vals),cluster = rep("Glia_1",length(glia_1_vals)),gene = rep(marker,length(glia_1_vals)))
  glia_2_vals = values[rownames(data) == marker, data$cluster == "Glia_2"]
  glia_2_res = data.frame(expression = as.numeric(glia_2_vals),cluster = rep("Glia_2",length(glia_2_vals)),gene = rep(marker,length(glia_2_vals)))
  results = rbind(results,glia_1_res,glia_2_res)
}


pdf("violin_plot_glial_markers.pdf",height = 2, width = 6)
p = ggplot(results, aes(x=gene, y=expression,fill = cluster)) +  
  geom_violin(position=position_dodge(0.5)) + geom_boxplot(width=0.1,position=position_dodge(0.5)) + xlab('') + 
  ylab('normalised expression') + scale_fill_manual(values=cluster_cols) + theme(legend.title = element_blank())
print(p) 
dev.off()

## UMAP of known markers
for (marker in markers){
  pdf(paste0("umap_",marker,".pdf"))
  p = FeaturePlot(object = data, features = marker, pt.size = 3,
                  label = F)
  p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20),plot.title = element_blank()) + scale_color_viridis() +
    xlim(xlim) + ylim(ylim)
  print(p)
  dev.off()
}


pdf("dotplot_known_markers.pdf", height = 3.5, width =6)
p = DotPlot(data[,data$cluster %in% c("Glia_1","Glia_2")], features = markers) + RotatedAxis() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_text(size=10),legend.text=element_text(size=8)) +
  scale_color_viridis()
print(p)
dev.off()


## DE list Glia_1 vs Glia_2

#run DE using Wilcoxon test
glia_de = FindMarkers(data, ident.1 = "Glia_1", ident.2 = "Glia_2")
write.csv(glia_de,"glia_DE.csv")

#segregate up and down regulated genes
selected1 = subset(glia_de, glia_de$avg_logFC > 0)
selected2 = subset(glia_de, glia_de$avg_logFC < 0)

#order genes by logFC
selected1 = selected1[order(-selected1$avg_logFC),]
selected2 = selected2[order(selected2$avg_logFC),]

#select top 20 up and down regulated genes
genes = c(rownames(selected1)[1:20],rownames(selected2)[1:20])
genes = genes[!(is.na(genes))]

#create dotplot
pdf("dotplot_glia1_glia2_viridis_lfc_order.pdf", height = 3.5, width = 12)
p = DotPlot(data, features = genes) + RotatedAxis() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_text(size=10),legend.text=element_text(size=8)) +
  scale_color_viridis()
print(p)
dev.off()


## DE list for naive vs H. poly glia

#run DE using Mast with number of features as a latent variable
glia_ctrl_inf = FindMarkers(data, ident.1 = "Control", ident.2 = "H. poly", test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(glia_ctrl_inf,"glia_condition_markers_mast_nFeature.csv")

#segregate up and down regulated genes
selected1 = subset(glia_ctrl_inf, glia_ctrl_inf$avg_logFC > 0)
selected2 = subset(glia_ctrl_inf, glia_ctrl_inf$avg_logFC < 0)

#order genes by logFC
selected1 = selected1[order(-selected1$avg_logFC),]
selected2 = selected2[order(selected2$avg_logFC),]

#select top 20 up and down regulated genes
genes = c(rownames(selected1)[1:20],rownames(selected2)[1:20])
genes = genes[!(is.na(genes))]

#create dotplot
pdf("dotplot_glia_condition_mast_lfc_order_viridis_nFeature.pdf", height = 3.5, width = 12)
p = DotPlot(data, features = genes) + RotatedAxis() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_text(size=10),legend.text=element_text(size=8)) +
  scale_color_viridis()
print(p)
dev.off()


## Fisher exact test - null hyopothesis there is no relationship between cluster and condition

contin_tab = table(data$cluster,data$condition_label)
fisher.test(contin_tab)

# Fisher's Exact Test for Count Data
# 
# data:  contin_tab
# p-value = 2.283e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.677133 10.874061
# sample estimates:
# odds ratio 
#   5.317458 
table(data$condition)

## write genes to file so can be used as a background for gene set overrepresentation analysis

write.csv(rownames(data),"background.csv",row.names = F,quote = F)

## create umap with batch names for supplementary materials

data$batch_name = data$batch
data$batch_name[data$batch == "HP_1"] = "H. poly 1"
data$batch_name[data$batch == "HP_2"] = "H. poly 2"
data$batch_name[data$batch == "HP_3"] = "H. poly 3"
data$batch_name[data$batch == "C_1_"] = "Control 1"
data$batch_name[data$batch == "C_2_"] = "Control 2"

pdf("batch_umap2.pdf")
p = DimPlot(data, reduction = "umap",  group.by = "batch_name", pt.size = 3,
            label = F,label.size=5)
p = p + theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) 
print(p)
dev.off()

data@active.ident = data$condition_label

## create heatmap for response to reviewers

selected1 = subset(glia_de, glia_de$avg_logFC > 0)
selected2 = subset(glia_de, glia_de$avg_logFC < 0)

selected1 = selected1[order(-selected1$avg_logFC),]
selected2 = selected2[order(selected2$avg_logFC),]

genes = c(as.character(selected1$X[1:20]),as.character(selected2$X[1:20]))
genes = genes[!(is.na(genes))]

mat_norm = GetAssayData(object = data, slot = "scale.data")
mat_norm = as.matrix(mat_norm[genes,])
ordering = data[["seurat_clusters"]]

clusters = data@meta.data[,c("seurat_clusters","condition_label")][order(ordering),]
mat_norm = mat_norm[,order(ordering)]

#create vector with levels of object@ident
identities = levels(data@active.ident)

#cut at a zscore of 2.5
mat_norm[mat_norm > 2.5] = 2.5 
extreme = (max(abs(min(mat_norm)),abs(max(mat_norm))))
breaks = seq(from=-extreme, to=max(extreme), length.out=100)

ann_colours = list()
ann_colours$seurat_clusters = c("#0072B2","#F0E442")
ann_colours$condition_label =  c("#009E73", "#E69F00")

names(ann_colours$seurat_clusters) = levels(clusters$seurat_clusters)
names(ann_colours$condition_label) = levels(clusters$condition_label)
colour = colorRampPalette(c("blue", "black", "yellow"))(100)

png("heatmap_wilcox_markers_clusters.png", height = 5, width = 8, units = "in", res = 300)
pheatmap(mat_norm, cluster_rows = F, cluster_cols = F,annotation_col = clusters,show_rownames = T, show_colnames = F, color= colour,fontsize = 10,annotation_colors = ann_colours, breaks = breaks)
dev.off()

## save data

saveRDS(data,"seurat_object_final.RDS")
