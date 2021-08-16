library(Seurat)

## edit seurat object to remove redundant metadata,
# add colours to metadata, and rename columns for 
# interpretability

data = readRDS("../output/seurat_object_final.RDS")
data$condition = as.character(data$condition_label)

##rename metadata columns
data$condition[data$condition == "Control"] = "Naive"
data$cluster[data$cluster == "Glia_1"] = "EGC1"
data$cluster[data$cluster == "Glia_2"] = "EGC2"

##add colours for condition
data$condition_colors = data$condition
data$condition_colors[data$condition == "Naive"] = "#009E73"
data$condition_colors[data$condition == "H. poly"] = "#E69F00"

##add colours for cluster
data$cluster_colors = data$cluster
data$cluster_colors[data$cluster == "EGC1"] = "#0072B2"
data$cluster_colors[data$cluster == "EGC2"] = "#F0E442"

##rename interferon response gene columns 
data$inteferon_gamma_response_genes = data$response_to_ifn
data$inteferon_gamma_response_genes_scaled =  data$response_to_ifn_scaled

##rename batches
data$batch[data$batch == "C_1_"] = "Naive 1"
data$batch[data$batch == "C_2_"] = "Naive 2"
data$batch[data$batch == "HP_1"] = "H. poly 1"
data$batch[,data$batch == "HP_2"] = "H. poly 2"
data$batch[data$batch == "HP_3"] = "H. poly 3"

##only keep necessary metadata columms
data@meta.data = data@meta.data[,c("cluster", "cluster_colors","condition","condition_colors","batch","inteferon_gamma_response_genes","inteferon_gamma_response_genes_scaled")]

saveRDS(data,"../output/EGCsHPoly.rds")